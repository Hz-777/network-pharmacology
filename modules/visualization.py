"""
Visualization module — publication-quality figures.
Styles inspired by Nature/Cell paper conventions.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib.patheffects as pe
import networkx as nx
import plotly.graph_objects as go
from matplotlib_venn import venn2, venn3
import io, base64

# ── Font setup ────────────────────────────────────────────────────────────────
import matplotlib.font_manager as fm
fm.fontManager.__init__()
_avail = {f.name for f in fm.fontManager.ttflist}

# 黑体 priority list  (heiti / gothic)
_HEITI = [
    "SimHei",                                   # Windows / most Linux
    "STHeiti", "Heiti TC", "Heiti SC",          # macOS
    "Microsoft YaHei",                          # Windows
    "WenQuanYi Zen Hei", "WenQuanYi Micro Hei",# Linux with wqy packages
    "Noto Sans CJK SC", "Noto Sans SC",         # Noto
    "PingFang SC", "PingFang HK",               # macOS fallback
    "Arial Unicode MS",
]
# 宋体 priority list  (songti / serif)
_SONGTI = [
    "SimSun", "NSimSun",                        # Windows
    "Songti SC", "STSong",                      # macOS
    "FangSong",                                 # Windows
    "Noto Serif CJK SC",                        # Noto
    "AR PL UMing CN",                           # Linux open-source
]

_heiti_avail  = [f for f in _HEITI  if f in _avail]
_songti_avail = [f for f in _SONGTI if f in _avail]
_any_cjk      = (_heiti_avail or _songti_avail or ["DejaVu Sans"])

# Default: 黑体
plt.rcParams.update({
    "font.family":        (_heiti_avail or _any_cjk) + ["DejaVu Sans"],
    "axes.unicode_minus": False,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.grid":          False,
    "figure.dpi":         150,
    "savefig.dpi":        180,
    "savefig.bbox":       "tight",
    "savefig.facecolor":  "white",
})

# ── Palette ───────────────────────────────────────────────────────────────────
C_RED    = "#C0392B"
C_BLUE   = "#2471A3"
C_GREEN  = "#1A7A4A"
C_ORANGE = "#D35400"
C_GRAY   = "#566573"

CMAP_HEAT  = LinearSegmentedColormap.from_list("heat",  ["#FEF9E7","#F39C12","#C0392B"])
CMAP_BLUE  = LinearSegmentedColormap.from_list("blues", ["#EBF5FB","#2E86C1","#1B2631"])
CMAP_GREEN = LinearSegmentedColormap.from_list("greens",["#EAFAF1","#1E8449","#145A32"])


def _spine_off(ax):
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color("#CCCCCC")
    ax.spines["bottom"].set_color("#CCCCCC")


# ─────────────────────────────────────────────────────────────────────────────
# 1. VENN DIAGRAM
# ─────────────────────────────────────────────────────────────────────────────

def plot_venn(sets: dict, title: str = "靶点韦恩图") -> plt.Figure:
    labels = list(sets.keys())
    values = [sets[k] for k in labels]
    palette = ["#2471A3", "#C0392B", "#1A7A4A"]

    fig, ax = plt.subplots(figsize=(8, 7), facecolor="white")
    fig.subplots_adjust(top=0.88)

    if len(sets) == 2:
        v = venn2(values, set_labels=["", ""], set_colors=palette[:2], alpha=0.55, ax=ax)
    elif len(sets) == 3:
        v = venn3(values, set_labels=["", "", ""], set_colors=palette[:3], alpha=0.55, ax=ax)
    else:
        ax.text(0.5, 0.5, "需要 2-3 个集合", ha="center", va="center")
        return fig

    for subset_id in (v.subset_labels or []):
        if subset_id:
            subset_id.set_fontsize(15)
            subset_id.set_fontweight("bold")
            subset_id.set_color("white")
            subset_id.set_path_effects([pe.withStroke(linewidth=2, foreground="#333333")])

    legend_handles = [
        mpatches.Patch(facecolor=palette[i], alpha=0.75, edgecolor="white",
                       linewidth=2, label=labels[i])
        for i in range(len(labels))
    ]
    ax.legend(handles=legend_handles, loc="lower center",
              bbox_to_anchor=(0.5, -0.06), ncol=len(labels),
              fontsize=12, frameon=False)

    inter = set.intersection(*values)
    ax.text(0.5, 0.02,
            f"交集靶点：{len(inter)} 个",
            transform=ax.transAxes, ha="center", va="bottom",
            fontsize=12, color=C_RED, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.4", fc="#FDEDEC", ec=C_RED, lw=1))

    ax.set_title(title, fontsize=17, fontweight="bold", color="#1C2833", pad=18)
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 2. PPI NETWORK  — clearer edges, hub-only labels, no overlap
# ─────────────────────────────────────────────────────────────────────────────

def plot_ppi_network(
    ppi_df: pd.DataFrame,
    centrality_df: pd.DataFrame,
    top_n: int = 20,
    title: str = "PPI 蛋白互作网络",
) -> plt.Figure:
    if ppi_df is None or ppi_df.empty:
        fig, ax = plt.subplots(figsize=(10, 8), facecolor="white")
        ax.text(0.5, 0.5, "无 PPI 数据", ha="center", va="center", fontsize=16, color=C_GRAY)
        ax.axis("off")
        return fig

    cols = ppi_df.columns.tolist()
    src  = next((c for c in cols if "preferredName_A" in c), cols[0])
    tgt  = next((c for c in cols if "preferredName_B" in c), cols[1])
    scr  = next((c for c in cols if "score" in c.lower()), None)

    G = nx.from_pandas_edgelist(ppi_df, source=src, target=tgt)

    # Keep hub subgraph
    if centrality_df is not None and not centrality_df.empty and len(G.nodes) > top_n * 2:
        hubs = centrality_df.head(top_n)["Gene"].tolist()
        keep = set(hubs)
        for g in hubs:
            if g in G:
                keep.update(list(G.neighbors(g))[:4])
        G = G.subgraph(list(keep)).copy()

    if len(G.nodes) == 0:
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.text(0.5, 0.5, "网络为空", ha="center", va="center", fontsize=14)
        ax.axis("off")
        return fig

    # Layout — spread nodes well
    if len(G.nodes) <= 60:
        pos = nx.kamada_kawai_layout(G)
    else:
        pos = nx.spring_layout(G, k=3.0 / np.sqrt(len(G.nodes)), seed=42, iterations=60)

    # Node properties
    hub_set = (set(centrality_df.head(top_n)["Gene"].tolist())
               if centrality_df is not None and not centrality_df.empty else set())
    degree  = dict(G.degree())
    max_deg = max(degree.values()) if degree else 1

    node_list   = list(G.nodes)
    norm_degree = [degree.get(n, 1) / max_deg for n in node_list]
    node_colors = [CMAP_HEAT(v * 0.85 + 0.05) for v in norm_degree]
    node_sizes  = [350 + degree.get(n, 1) * 200 for n in node_list]

    # Edge width & alpha — higher confidence = darker, wider
    if scr:
        score_map = {}
        for _, row in ppi_df.iterrows():
            s = float(row[scr]) / 1000.0
            score_map[(row[src], row[tgt])] = s
            score_map[(row[tgt], row[src])] = s
        edge_w = [max(0.8, score_map.get(e, 0.5) * 3.0)  for e in G.edges]
        edge_a = [max(0.40, score_map.get(e, 0.5) * 0.75) for e in G.edges]
    else:
        edge_w = [1.2] * len(G.edges)
        edge_a = [0.50] * len(G.edges)

    fig, ax = plt.subplots(figsize=(14, 12), facecolor="white")

    # Draw edges — darker, more visible  (#5D6D7E is a clear medium-gray)
    for idx, (u, v) in enumerate(G.edges):
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color="#4A5568", linewidth=edge_w[idx],
                alpha=edge_a[idx], solid_capstyle="round", zorder=1)

    # Draw nodes
    ax.scatter(
        [pos[n][0] for n in node_list],
        [pos[n][1] for n in node_list],
        s=node_sizes,
        color=node_colors,
        edgecolors="white",
        linewidths=2.0,
        zorder=3,
        alpha=0.93,
    )

    # Labels — only for hub genes; offset outside node center to avoid overlap
    labeled_nodes = set(
        sorted(node_list, key=lambda n: degree.get(n, 0), reverse=True)[:min(top_n, 15)]
    )
    for n in labeled_nodes:
        x, y = pos[n]
        is_hub = n in hub_set
        fs = 9.0 if is_hub else 7.5
        fw = "bold" if is_hub else "normal"
        # Offset direction: push label away from graph center
        cx = np.mean([p[0] for p in pos.values()])
        cy = np.mean([p[1] for p in pos.values()])
        dx, dy = x - cx, y - cy
        norm = max(np.sqrt(dx**2 + dy**2), 1e-6)
        off  = (node_sizes[node_list.index(n)] ** 0.5) * 0.007
        lx, ly = x + dx / norm * off, y + dy / norm * off
        ax.text(lx, ly, n,
                ha="center", va="center",
                fontsize=fs, fontweight=fw,
                color="white" if is_hub else "#1C2833",
                zorder=5,
                path_effects=[pe.withStroke(linewidth=2.5, foreground="#1C2833")] if is_hub
                              else [pe.withStroke(linewidth=1.5, foreground="white")])

    # Colorbar
    sm = cm.ScalarMappable(cmap=CMAP_HEAT, norm=Normalize(0, max_deg))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.42, pad=0.01, aspect=20)
    cbar.set_label("节点连接度 (Degree)", fontsize=10, color=C_GRAY)
    cbar.ax.tick_params(labelsize=9, color=C_GRAY)
    cbar.outline.set_edgecolor("#CCCCCC")

    hub_patch   = mpatches.Patch(color=CMAP_HEAT(0.9), label=f"核心靶点 Top {top_n}")
    other_patch = mpatches.Patch(color=CMAP_HEAT(0.2), label="其他靶点")
    ax.legend(handles=[hub_patch, other_patch], loc="upper left",
              fontsize=10, frameon=True, framealpha=0.9,
              edgecolor="#DDDDDD", fancybox=True)

    ax.set_title(title, fontsize=17, fontweight="bold", color="#1C2833", pad=16)
    ax.axis("off")
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 3. KEGG BUBBLE CHART
# ─────────────────────────────────────────────────────────────────────────────

def plot_kegg_bubbles(kegg_df: pd.DataFrame, top_n: int = 20,
                      title: str = "KEGG 通路富集分析") -> plt.Figure:
    if kegg_df is None or kegg_df.empty:
        fig, ax = plt.subplots(figsize=(10, 6), facecolor="white")
        ax.text(0.5, 0.5, "无 KEGG 数据", ha="center", va="center", fontsize=14, color=C_GRAY)
        ax.axis("off")
        return fig

    df = kegg_df.head(top_n).copy()
    df["neg_log_p"]  = -np.log10(df["P_value"].clip(lower=1e-300))
    df["gene_count"] = df["Genes"].apply(lambda x: len(str(x).split(";")) if pd.notna(x) else 1)
    df["Term_short"] = df["Term"].apply(
        lambda x: (x.split("__")[-1] if "__" in str(x) else str(x))[:50]
    )
    df = df.sort_values("neg_log_p", ascending=True).reset_index(drop=True)

    # Enough row height to prevent label overlap
    row_h  = 0.65
    height = max(6, len(df) * row_h + 2)
    fig, ax = plt.subplots(figsize=(13, height), facecolor="white")

    for i in range(len(df)):
        ax.axhspan(i - 0.48, i + 0.48,
                   color="#F2F3F4" if i % 2 == 0 else "white", zorder=0)

    norm   = Normalize(vmin=df["Combined_score"].min(), vmax=df["Combined_score"].max())
    sizes  = (df["gene_count"] / df["gene_count"].max() * 360 + 80).values
    sc = ax.scatter(df["neg_log_p"], range(len(df)),
                    s=sizes, c=df["Combined_score"].values,
                    cmap=CMAP_HEAT, norm=norm,
                    edgecolors="white", linewidths=1.0,
                    zorder=3, alpha=0.92)

    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(row["neg_log_p"], i, str(int(row["gene_count"])),
                ha="center", va="center", fontsize=7.5, color="white",
                fontweight="bold", zorder=4)

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["Term_short"], fontsize=9.5)
    ax.set_ylim(-0.6, len(df) - 0.4)
    ax.set_xlabel("−log₁₀ (P value)", fontsize=12, color=C_GRAY, labelpad=8)
    ax.tick_params(axis="x", labelsize=10, color="#BBBBBB")
    ax.tick_params(axis="y", left=False)

    p05 = -np.log10(0.05)
    ax.axvline(p05, color=C_RED, linestyle="--", linewidth=1.2, alpha=0.7)
    ax.text(p05 + 0.02, len(df) - 0.8, "P = 0.05",
            color=C_RED, fontsize=8.5, va="top")

    _spine_off(ax)
    ax.spines["left"].set_visible(False)

    cbar = fig.colorbar(sc, ax=ax, shrink=0.5, pad=0.02, aspect=22)
    cbar.set_label("Combined Score", fontsize=10, color=C_GRAY)
    cbar.ax.tick_params(labelsize=9)
    cbar.outline.set_edgecolor("#CCCCCC")

    max_gc = df["gene_count"].max()
    for gc, label in [(5, "5 基因"), (15, "15 基因"), (30, "30 基因")]:
        if gc <= max_gc:
            ax.scatter([], [], s=gc / max_gc * 360 + 80, c=[C_GRAY],
                       alpha=0.6, edgecolors="white", label=label)
    ax.legend(title="基因数", title_fontsize=9, fontsize=8.5,
              loc="lower right", frameon=True, framealpha=0.9,
              edgecolor="#DDDDDD", fancybox=True)

    ax.set_title(title, fontsize=16, fontweight="bold", color="#1C2833", pad=14)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 4. GO BAR CHART
# ─────────────────────────────────────────────────────────────────────────────

def plot_go_barplot(go_df: pd.DataFrame, category: str = "GO_BP",
                    top_n: int = 15) -> plt.Figure:
    label_map = {
        "GO_BP": ("生物过程 (Biological Process)", CMAP_BLUE),
        "GO_CC": ("细胞组分 (Cellular Component)", CMAP_GREEN),
        "GO_MF": ("分子功能 (Molecular Function)",
                  LinearSegmentedColormap.from_list("orange", ["#FEF5E7","#E67E22","#784212"])),
    }
    cat_label, cmap = label_map.get(category, (category, CMAP_BLUE))

    if go_df is None or go_df.empty:
        fig, ax = plt.subplots(figsize=(10, 6), facecolor="white")
        ax.text(0.5, 0.5, f"无 {category} 数据", ha="center", va="center",
                fontsize=14, color=C_GRAY)
        ax.axis("off")
        return fig

    df = go_df.head(top_n).copy()
    df["neg_log_p"]  = -np.log10(df["P_value"].clip(lower=1e-300))
    df["gene_count"] = df["Genes"].apply(lambda x: len(str(x).split(";")) if pd.notna(x) else 1)
    df["Term_short"] = df["Term"].apply(
        lambda x: x.split("(GO:")[0].strip()[:55] if "(GO:" in str(x) else str(x)[:55]
    )
    df = df.sort_values("neg_log_p", ascending=True).reset_index(drop=True)

    row_h  = 0.60
    height = max(5.5, len(df) * row_h + 2)
    fig, ax = plt.subplots(figsize=(13, height), facecolor="white")

    norm   = Normalize(vmin=df["neg_log_p"].min(), vmax=df["neg_log_p"].max())
    colors = [cmap(norm(v) * 0.8 + 0.15) for v in df["neg_log_p"]]

    ax.barh(range(len(df)), df["neg_log_p"],
            color=colors, edgecolor="white", linewidth=0.6,
            height=0.62, zorder=2)

    # n= labels offset slightly right of bar end to avoid overlap with bar
    x_max = df["neg_log_p"].max()
    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(row["neg_log_p"] + x_max * 0.015, i,
                f"n={int(row['gene_count'])}",
                va="center", fontsize=8, color=C_GRAY)

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["Term_short"], fontsize=9.5)
    ax.set_ylim(-0.55, len(df) - 0.45)
    ax.set_xlabel("−log₁₀ (P value)", fontsize=12, color=C_GRAY, labelpad=8)
    ax.tick_params(axis="x", labelsize=10, color="#BBBBBB")
    ax.tick_params(axis="y", left=False)
    # x-axis right margin for n= labels
    ax.set_xlim(right=x_max * 1.18)

    p05 = -np.log10(0.05)
    ax.axvline(p05, color=C_RED, linestyle="--", linewidth=1.2, alpha=0.7)
    ax.text(p05 + x_max * 0.01, len(df) - 0.85, "P = 0.05",
            color=C_RED, fontsize=8.5, va="top")

    _spine_off(ax)
    ax.spines["left"].set_visible(False)

    ax.set_title(f"GO 富集分析 — {cat_label}  (Top {top_n})",
                 fontsize=15, fontweight="bold", color="#1C2833", pad=14)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 5. COMPONENT-TARGET-PATHWAY NETWORK (Plotly interactive)
#    Three-tier: Compound → Target → Pathway  (no herb hub node)
# ─────────────────────────────────────────────────────────────────────────────

def plot_network_plotly(
    compounds: list,
    compound_targets: dict,
    intersection_targets: list,
    top_pathways: list,
    title: str = "成分-靶点-通路网络图",
) -> go.Figure:

    G = nx.Graph()
    node_types = {}

    # Tier 1: compounds
    for comp in compounds:
        G.add_node(comp)
        node_types[comp] = "compound"
        # Tier 2: intersection targets linked to this compound
        for t in compound_targets.get(comp, []):
            if t in intersection_targets:
                G.add_node(t)
                node_types[t] = "target"
                G.add_edge(comp, t)

    # Tier 3: pathways linked to targets
    for path in top_pathways[:12]:
        short = path[:45]
        G.add_node(short)
        node_types[short] = "pathway"
        for t in intersection_targets[:10]:
            if t in G.nodes:
                G.add_edge(t, short)

    if len(G.nodes) == 0:
        return go.Figure()

    # Hierarchical-ish layout using multipartite, fallback to kamada_kawai
    try:
        for n, t in node_types.items():
            G.nodes[n]["layer"] = {"compound": 0, "target": 1, "pathway": 2}[t]
        pos = nx.multipartite_layout(G, subset_key="layer", scale=2.5, align="vertical")
    except Exception:
        pos = nx.kamada_kawai_layout(G) if len(G.nodes) <= 80 \
              else nx.spring_layout(G, k=3.5, seed=42)

    palette = {
        "compound": dict(color="#2471A3", size=24, symbol="circle"),
        "target":   dict(color="#1A7A4A", size=18, symbol="circle"),
        "pathway":  dict(color="#D35400", size=22, symbol="diamond"),
    }
    label_cn = {
        "compound": "活性成分",
        "target":   "交集靶点",
        "pathway":  "KEGG通路",
    }

    # Edges
    edge_x, edge_y = [], []
    for u, v in G.edges:
        edge_x += [pos[u][0], pos[v][0], None]
        edge_y += [pos[u][1], pos[v][1], None]

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y, mode="lines",
        line=dict(width=1.0, color="rgba(120,130,145,0.40)"),
        hoverinfo="none", showlegend=False,
    )

    traces = [edge_trace]
    for ntype, style in palette.items():
        nodes = [n for n, t in node_types.items() if t == ntype]
        if not nodes:
            continue
        deg    = [G.degree(n) for n in nodes]
        max_d  = max(deg) if deg else 1
        sizes  = [style["size"] + d / max_d * style["size"] * 0.9 for d in deg]
        hover  = [f"<b>{n}</b><br>类型: {label_cn[ntype]}<br>连接数: {G.degree(n)}"
                  for n in nodes]

        traces.append(go.Scatter(
            x=[pos[n][0] for n in nodes],
            y=[pos[n][1] for n in nodes],
            mode="markers+text",
            marker=dict(
                size=sizes,
                color=style["color"],
                symbol=style["symbol"],
                opacity=0.90,
                line=dict(width=1.5, color="white"),
            ),
            text=[n[:32] for n in nodes],
            textposition="top center",
            textfont=dict(size=9 if ntype == "pathway" else 8, color="#1C2833"),
            hovertemplate="%{customdata}<extra></extra>",
            customdata=hover,
            name=label_cn[ntype],
        ))

    fig = go.Figure(data=traces, layout=go.Layout(
        title=dict(
            text=f"<b>{title}</b>",
            font=dict(size=17, color="#1C2833", family="Arial"),
            x=0.5, xanchor="center",
        ),
        showlegend=True,
        legend=dict(
            orientation="v", x=1.01, y=0.98,
            bgcolor="rgba(255,255,255,0.9)",
            bordercolor="#DDDDDD", borderwidth=1,
            font=dict(size=12),
        ),
        hovermode="closest",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, visible=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, visible=False),
        paper_bgcolor="white",
        plot_bgcolor="#F8F9FA",
        height=750,
        margin=dict(l=20, r=170, t=70, b=20),
    ))
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Utility
# ─────────────────────────────────────────────────────────────────────────────

def apply_cfg(cfg: dict):
    """Apply user font/style config to matplotlib rcParams globally."""
    if not cfg:
        return

    # Font style switch
    font_style = cfg.get("font_style", "黑体")
    if font_style == "宋体" and _songti_avail:
        plt.rcParams["font.family"] = _songti_avail + ["DejaVu Sans"]
    else:
        plt.rcParams["font.family"] = (_heiti_avail or _any_cjk) + ["DejaVu Sans"]

    fs = cfg.get("font_scale", 1.0)
    plt.rcParams.update({
        "font.size":        10 * fs,
        "axes.titlesize":   15 * fs,
        "axes.labelsize":   12 * fs,
        "xtick.labelsize":  10 * fs,
        "ytick.labelsize":  10 * fs,
        "legend.fontsize":  10 * fs,
        "figure.dpi":       cfg.get("dpi", 150),
        "savefig.dpi":      cfg.get("dpi", 180),
        "axes.unicode_minus": False,
    })


def fig_to_base64(fig: plt.Figure, fmt: str = "png", dpi: int = 180) -> str:
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=dpi, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    buf.seek(0)
    plt.close(fig)
    return base64.b64encode(buf.read()).decode("utf-8")


def fig_to_bytes(fig: plt.Figure, fmt: str = "png", dpi: int = 180) -> bytes:
    buf = io.BytesIO()
    fig.savefig(buf, format=fmt, dpi=dpi, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    buf.seek(0)
    plt.close(fig)
    return buf.getvalue()
