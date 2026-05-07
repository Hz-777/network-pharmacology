"""Visualization module: Venn diagrams, PPI networks, bubble plots."""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx
import plotly.graph_objects as go
import plotly.express as px
from matplotlib_venn import venn2, venn3
from typing import Optional, Tuple
import io
import base64


import matplotlib.font_manager as fm
_available = {f.name for f in fm.fontManager.ttflist}
_cjk_fonts = [f for f in ["Arial Unicode MS", "PingFang SC", "Heiti SC", "Microsoft YaHei", "SimHei"] if f in _available]
plt.rcParams["font.family"] = _cjk_fonts + ["DejaVu Sans"] if _cjk_fonts else ["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


# --- Venn Diagram ---

def plot_venn(sets: dict, title: str = "靶点韦恩图") -> plt.Figure:
    """
    Draw Venn diagram for 2 or 3 sets.
    sets: dict of {label: set_of_genes}
    """
    fig, ax = plt.subplots(figsize=(7, 6), facecolor="white")
    labels = list(sets.keys())
    values = [sets[k] for k in labels]

    colors = ["#4C72B0", "#DD8452", "#55A868"]

    if len(sets) == 2:
        v = venn2(
            values,
            set_labels=labels,
            set_colors=colors[:2],
            alpha=0.6,
            ax=ax,
        )
    elif len(sets) == 3:
        v = venn3(
            values,
            set_labels=labels,
            set_colors=colors[:3],
            alpha=0.6,
            ax=ax,
        )

    ax.set_title(title, fontsize=16, fontweight="bold", pad=15)
    fig.tight_layout()
    return fig


# --- PPI Network ---

def plot_ppi_network(
    ppi_df: pd.DataFrame,
    centrality_df: pd.DataFrame,
    top_n: int = 20,
    title: str = "PPI 蛋白互作网络",
) -> plt.Figure:
    """Draw PPI network graph highlighting hub genes."""
    if ppi_df.empty:
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, "无PPI数据", ha="center", va="center", fontsize=16)
        return fig

    # Find edge columns
    cols = ppi_df.columns.tolist()
    src_col = next((c for c in cols if "preferredName_A" in c or c == cols[0]), cols[0])
    tgt_col = next((c for c in cols if "preferredName_B" in c or c == cols[1]), cols[1])
    score_col = next((c for c in cols if "score" in c.lower()), None)

    G = nx.from_pandas_edgelist(ppi_df, source=src_col, target=tgt_col)

    # Limit to connected subgraph around hub genes
    if not centrality_df.empty and len(G.nodes) > top_n * 2:
        hub_genes = centrality_df.head(top_n)["Gene"].tolist()
        neighbors = set(hub_genes)
        for g in hub_genes:
            if g in G:
                neighbors.update(list(G.neighbors(g))[:5])
        G = G.subgraph(list(neighbors)).copy()

    fig, ax = plt.subplots(figsize=(14, 12), facecolor="white")

    pos = nx.spring_layout(G, k=2.5 / np.sqrt(max(len(G.nodes), 1)), seed=42)

    # Node sizes by degree
    degrees = dict(G.degree())
    node_sizes = [max(300, degrees.get(n, 1) * 150) for n in G.nodes]

    # Color hub genes
    hub_set = set(centrality_df.head(top_n)["Gene"].tolist()) if not centrality_df.empty else set()
    node_colors = ["#E74C3C" if n in hub_set else "#3498DB" for n in G.nodes]

    # Draw edges
    edge_weights = []
    if score_col:
        edge_map = {}
        for _, row in ppi_df.iterrows():
            edge_map[(row[src_col], row[tgt_col])] = row[score_col]
            edge_map[(row[tgt_col], row[src_col])] = row[score_col]
        edge_weights = [edge_map.get(e, 400) / 1000.0 for e in G.edges]
    else:
        edge_weights = [0.5] * len(G.edges)

    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.3, width=edge_weights, edge_color="#95A5A6")
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes, node_color=node_colors, alpha=0.85)
    nx.draw_networkx_labels(G, pos, ax=ax, font_size=7, font_color="white", font_weight="bold")

    # Legend
    hub_patch = mpatches.Patch(color="#E74C3C", label=f"核心靶点 (Top {top_n})")
    other_patch = mpatches.Patch(color="#3498DB", label="其他靶点")
    ax.legend(handles=[hub_patch, other_patch], loc="upper left", fontsize=11)

    ax.set_title(title, fontsize=16, fontweight="bold", pad=15)
    ax.axis("off")
    fig.tight_layout()
    return fig


# --- Bubble/Bar plots for GO/KEGG ---

def plot_kegg_bubbles(kegg_df: pd.DataFrame, top_n: int = 20, title: str = "KEGG 通路富集分析") -> plt.Figure:
    """Draw bubble chart for KEGG enrichment results."""
    if kegg_df.empty:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, "无KEGG数据", ha="center", va="center", fontsize=14)
        return fig

    df = kegg_df.head(top_n).copy()
    df["neg_log_p"] = -np.log10(df["P_value"].clip(lower=1e-300))
    df["gene_count"] = df["Genes"].apply(lambda x: len(str(x).split(";")) if pd.notna(x) else 0)

    # Shorten term names
    df["Term_short"] = df["Term"].apply(lambda x: x.split("__")[-1][:60] if "__" in str(x) else str(x)[:60])
    df = df.sort_values("neg_log_p")

    fig, ax = plt.subplots(figsize=(12, max(6, len(df) * 0.45)), facecolor="white")

    scatter = ax.scatter(
        df["neg_log_p"],
        range(len(df)),
        s=df["gene_count"] * 25 + 50,
        c=df["Combined_score"],
        cmap="RdYlBu_r",
        alpha=0.8,
        edgecolors="gray",
        linewidths=0.5,
    )

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["Term_short"], fontsize=9)
    ax.set_xlabel("-log₁₀(P value)", fontsize=12)
    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    ax.axvline(x=-np.log10(0.05), color="red", linestyle="--", alpha=0.7, label="P=0.05")
    ax.legend(fontsize=9)

    plt.colorbar(scatter, ax=ax, label="Combined Score", shrink=0.6)
    fig.tight_layout()
    return fig


def plot_go_barplot(go_df: pd.DataFrame, category: str = "GO_BP", top_n: int = 15) -> plt.Figure:
    """Draw horizontal bar chart for GO enrichment results."""
    if go_df.empty:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, f"无 {category} 数据", ha="center", va="center", fontsize=14)
        return fig

    df = go_df.head(top_n).copy()
    df["neg_log_p"] = -np.log10(df["P_value"].clip(lower=1e-300))
    df["Term_short"] = df["Term"].apply(lambda x: x.split("(")[0][:55].strip())
    df = df.sort_values("neg_log_p")

    color_map = {
        "GO_BP": "#3498DB",
        "GO_CC": "#27AE60",
        "GO_MF": "#E67E22",
    }
    color = color_map.get(category, "#9B59B6")

    fig, ax = plt.subplots(figsize=(11, max(5, len(df) * 0.4)), facecolor="white")
    bars = ax.barh(range(len(df)), df["neg_log_p"], color=color, alpha=0.8, edgecolor="white")

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["Term_short"], fontsize=9)
    ax.set_xlabel("-log₁₀(P value)", fontsize=12)
    ax.set_title(f"{category} 富集分析 (Top {top_n})", fontsize=14, fontweight="bold", pad=12)
    ax.axvline(x=-np.log10(0.05), color="red", linestyle="--", alpha=0.7, label="P=0.05")
    ax.legend(fontsize=9)

    fig.tight_layout()
    return fig


# --- Component-Target-Pathway network (Plotly, interactive) ---

def plot_network_plotly(
    compounds: list,
    compound_targets: dict,
    intersection_targets: list,
    top_pathways: list,
    title: str = "成分-靶点-通路网络图",
) -> go.Figure:
    """
    Build an interactive Plotly network: Herb → Compounds → Targets → Pathways
    """
    G = nx.Graph()
    node_types = {}

    # Add nodes
    herb_node = "中药复方"
    G.add_node(herb_node)
    node_types[herb_node] = "herb"

    for comp in compounds:
        G.add_node(comp)
        node_types[comp] = "compound"
        G.add_edge(herb_node, comp)

        targets = compound_targets.get(comp, [])
        for t in targets:
            if t in intersection_targets:
                G.add_node(t)
                node_types[t] = "target"
                G.add_edge(comp, t)

    for pathway in top_pathways[:15]:
        G.add_node(pathway)
        node_types[pathway] = "pathway"
        # Connect targets to pathways (simplified)
        for t in intersection_targets[:5]:
            if t in G.nodes:
                G.add_edge(t, pathway)

    pos = nx.spring_layout(G, k=3.5, seed=42)

    # Build Plotly traces
    color_map = {
        "herb": "#E74C3C",
        "compound": "#3498DB",
        "target": "#27AE60",
        "pathway": "#F39C12",
    }
    size_map = {
        "herb": 30,
        "compound": 18,
        "target": 14,
        "pathway": 16,
    }

    edge_x, edge_y = [], []
    for e in G.edges:
        x0, y0 = pos[e[0]]
        x1, y1 = pos[e[1]]
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y, mode="lines",
        line=dict(width=0.7, color="#BDC3C7"),
        hoverinfo="none",
    )

    node_traces = []
    for ntype in ["herb", "compound", "target", "pathway"]:
        nodes = [n for n, t in node_types.items() if t == ntype]
        if not nodes:
            continue
        nx_list = [pos[n][0] for n in nodes]
        ny_list = [pos[n][1] for n in nodes]
        labels = [n[:30] for n in nodes]

        trace = go.Scatter(
            x=nx_list, y=ny_list, mode="markers+text",
            marker=dict(size=size_map[ntype], color=color_map[ntype], opacity=0.85, line=dict(width=1, color="white")),
            text=labels, textposition="top center", textfont=dict(size=8),
            name={"herb": "中药", "compound": "活性成分", "target": "靶点", "pathway": "通路"}[ntype],
            hoverinfo="text",
        )
        node_traces.append(trace)

    fig = go.Figure(
        data=[edge_trace] + node_traces,
        layout=go.Layout(
            title=dict(text=title, font=dict(size=18)),
            showlegend=True,
            hovermode="closest",
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            paper_bgcolor="white",
            plot_bgcolor="white",
            height=700,
        ),
    )
    return fig


def fig_to_base64(fig: plt.Figure) -> str:
    """Convert matplotlib figure to base64 PNG string for Streamlit display."""
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=120, bbox_inches="tight")
    buf.seek(0)
    plt.close(fig)
    return base64.b64encode(buf.read()).decode("utf-8")
