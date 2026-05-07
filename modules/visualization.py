"""
Visualization module — publication-quality figures.
Supports 5 journal-style palettes: NPG · NEJM · Lancet · Science · Cell
Each palette provides consistent colors across Venn / PPI / KEGG / GO / Network.
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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

_HEITI  = ["SimHei","STHeiti","Heiti TC","Microsoft YaHei",
           "WenQuanYi Zen Hei","WenQuanYi Micro Hei",
           "Noto Sans CJK SC","Noto Sans SC","PingFang SC","Arial Unicode MS"]
_SONGTI = ["SimSun","NSimSun","Songti SC","STSong",
           "FangSong","Noto Serif CJK SC","AR PL UMing CN"]

_heiti_avail  = [f for f in _HEITI  if f in _avail]
_songti_avail = [f for f in _SONGTI if f in _avail]
_any_cjk      = _heiti_avail or _songti_avail or ["DejaVu Sans"]

plt.rcParams.update({
    "font.family":        _heiti_avail + ["DejaVu Sans"] if _heiti_avail else ["DejaVu Sans"],
    "axes.unicode_minus": False,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.grid":          False,
    "figure.dpi":         150,
    "savefig.dpi":        180,
    "savefig.bbox":       "tight",
    "savefig.facecolor":  "white",
})

# ── Journal Palettes ──────────────────────────────────────────────────────────
def _cmap(name, colors):
    return LinearSegmentedColormap.from_list(name, colors)

PALETTES = {
    # ── Nature Publishing Group ───────────────────────────────────────────────
    "NPG": {
        "label":   "Nature (NPG)",
        "accent":  ["#E64B35","#4DBBD5","#00A087","#3C5488",
                    "#F39B7F","#8491B4","#91D1C2","#DC0000"],
        "venn":    ["#4DBBD5","#E64B35","#00A087"],
        "ppi":     _cmap("npg_ppi",  ["#FEF0D9","#F4A58A","#E64B35","#3C5488"]),
        "kegg":    _cmap("npg_kegg", ["#FEF0D9","#F4A58A","#E64B35"]),
        "go_bp":   _cmap("npg_bp",   ["#EEF4FF","#8491B4","#3C5488"]),
        "go_cc":   _cmap("npg_cc",   ["#E8FAF5","#91D1C2","#00A087"]),
        "go_mf":   _cmap("npg_mf",   ["#FFF3E0","#F39B7F","#E64B35"]),
        "net":     {"compound":"#3C5488","target":"#00A087","pathway":"#E64B35"},
        "edge":    "#4A5568",
        "ref_line":"#E64B35",
    },
    # ── New England Journal of Medicine ───────────────────────────────────────
    "NEJM": {
        "label":   "NEJM",
        "accent":  ["#BC3C29","#0072B5","#E18727","#20854E",
                    "#7876B1","#6F99AD","#FFDC91","#EE4C97"],
        "venn":    ["#0072B5","#BC3C29","#20854E"],
        "ppi":     _cmap("nejm_ppi",  ["#FFF5F0","#E18727","#BC3C29","#0072B5"]),
        "kegg":    _cmap("nejm_kegg", ["#FFF5F0","#E18727","#BC3C29"]),
        "go_bp":   _cmap("nejm_bp",   ["#EBF4FB","#6F99AD","#0072B5"]),
        "go_cc":   _cmap("nejm_cc",   ["#E8F5E9","#20854E","#145A32"]),
        "go_mf":   _cmap("nejm_mf",   ["#FFFCE0","#FFDC91","#E18727"]),
        "net":     {"compound":"#0072B5","target":"#20854E","pathway":"#BC3C29"},
        "edge":    "#4A5568",
        "ref_line":"#BC3C29",
    },
    # ── Lancet ────────────────────────────────────────────────────────────────
    "Lancet": {
        "label":   "Lancet",
        "accent":  ["#00468B","#ED0000","#42B540","#0099B4",
                    "#925E9F","#FDAF91","#AD002A","#ADB6B6"],
        "venn":    ["#0099B4","#ED0000","#42B540"],
        "ppi":     _cmap("lancet_ppi",  ["#FFF0F0","#FDAF91","#ED0000","#00468B"]),
        "kegg":    _cmap("lancet_kegg", ["#FFF0F0","#FDAF91","#ED0000"]),
        "go_bp":   _cmap("lancet_bp",   ["#E3F2FD","#0099B4","#00468B"]),
        "go_cc":   _cmap("lancet_cc",   ["#E8F5E9","#42B540","#1B5E20"]),
        "go_mf":   _cmap("lancet_mf",   ["#FFF3E0","#FDAF91","#AD002A"]),
        "net":     {"compound":"#00468B","target":"#42B540","pathway":"#ED0000"},
        "edge":    "#4A5568",
        "ref_line":"#ED0000",
    },
    # ── Science (AAAS) ────────────────────────────────────────────────────────
    "Science": {
        "label":   "Science (AAAS)",
        "accent":  ["#3B4992","#EE0000","#008B45","#631879",
                    "#008280","#BB0021","#5F559B","#A20056"],
        "venn":    ["#3B4992","#EE0000","#008B45"],
        "ppi":     _cmap("aaas_ppi",  ["#F0F1F9","#9DA7D6","#EE0000","#3B4992"]),
        "kegg":    _cmap("aaas_kegg", ["#FFF5F5","#FF8080","#EE0000"]),
        "go_bp":   _cmap("aaas_bp",   ["#E8EAF6","#5F559B","#3B4992"]),
        "go_cc":   _cmap("aaas_cc",   ["#E8F5E9","#008B45","#1B5E20"]),
        "go_mf":   _cmap("aaas_mf",   ["#FCE4EC","#BB0021","#631879"]),
        "net":     {"compound":"#3B4992","target":"#008B45","pathway":"#EE0000"},
        "edge":    "#4A5568",
        "ref_line":"#EE0000",
    },
    # ── Cell Press ────────────────────────────────────────────────────────────
    "Cell": {
        "label":   "Cell Press",
        "accent":  ["#D6604D","#4393C3","#5FAD56","#984EA3",
                    "#FF7F00","#E6AB02","#A6761D","#666666"],
        "venn":    ["#4393C3","#D6604D","#5FAD56"],
        "ppi":     _cmap("cell_ppi",  ["#FEF5EC","#F4A582","#D6604D","#4393C3"]),
        "kegg":    _cmap("cell_kegg", ["#FEF5EC","#F4A582","#D6604D"]),
        "go_bp":   _cmap("cell_bp",   ["#DEEBF7","#4393C3","#08519C"]),
        "go_cc":   _cmap("cell_cc",   ["#E5F5E0","#5FAD56","#238B45"]),
        "go_mf":   _cmap("cell_mf",   ["#FEF0D9","#FDBB84","#D7301F"]),
        "net":     {"compound":"#4393C3","target":"#5FAD56","pathway":"#D6604D"},
        "edge":    "#4A5568",
        "ref_line":"#D6604D",
    },
}

# Active palette — set by apply_cfg, default NPG
_PAL = PALETTES["NPG"]


def _get_pal() -> dict:
    return _PAL


def _spine_off(ax):
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color("#CCCCCC")
    ax.spines["bottom"].set_color("#CCCCCC")


# ─────────────────────────────────────────────────────────────────────────────
# 1. VENN DIAGRAM
# ─────────────────────────────────────────────────────────────────────────────

def plot_venn(sets: dict, title: str = "靶点韦恩图") -> plt.Figure:
    pal    = _get_pal()
    labels = list(sets.keys())
    values = [sets[k] for k in labels]
    colors = pal["venn"]

    fig, ax = plt.subplots(figsize=(8, 7), facecolor="white")
    fig.subplots_adjust(top=0.88)

    if len(sets) == 2:
        v = venn2(values, set_labels=["", ""], set_colors=colors[:2], alpha=0.55, ax=ax)
    elif len(sets) == 3:
        v = venn3(values, set_labels=["", "", ""], set_colors=colors[:3], alpha=0.55, ax=ax)
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
        mpatches.Patch(facecolor=colors[i], alpha=0.75, edgecolor="white",
                       linewidth=2, label=labels[i])
        for i in range(len(labels))
    ]
    ax.legend(handles=legend_handles, loc="lower center",
              bbox_to_anchor=(0.5, -0.06), ncol=len(labels),
              fontsize=12, frameon=False)

    inter = set.intersection(*values)
    ref   = pal["ref_line"]
    ax.text(0.5, 0.02, f"交集靶点：{len(inter)} 个",
            transform=ax.transAxes, ha="center", va="bottom",
            fontsize=12, color=ref, fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.4", fc=ref + "18", ec=ref, lw=1))

    ax.set_title(title, fontsize=17, fontweight="bold", color="#1C2833", pad=18)
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 2. PPI NETWORK
# ─────────────────────────────────────────────────────────────────────────────

def plot_ppi_network(
    ppi_df: pd.DataFrame,
    centrality_df: pd.DataFrame,
    top_n: int = 20,
    title: str = "PPI 蛋白互作网络",
) -> plt.Figure:
    pal = _get_pal()
    C_GRAY = "#566573"

    if ppi_df is None or ppi_df.empty:
        fig, ax = plt.subplots(figsize=(10, 8), facecolor="white")
        ax.text(0.5, 0.5, "无 PPI 数据", ha="center", va="center",
                fontsize=16, color=C_GRAY)
        ax.axis("off")
        return fig

    cols = ppi_df.columns.tolist()
    src  = next((c for c in cols if "preferredName_A" in c), cols[0])
    tgt  = next((c for c in cols if "preferredName_B" in c), cols[1])
    scr  = next((c for c in cols if "score" in c.lower()), None)

    G = nx.from_pandas_edgelist(ppi_df, source=src, target=tgt)

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

    if len(G.nodes) <= 60:
        pos = nx.kamada_kawai_layout(G)
    else:
        pos = nx.spring_layout(G, k=3.0 / np.sqrt(len(G.nodes)), seed=42, iterations=60)

    hub_set = (set(centrality_df.head(top_n)["Gene"].tolist())
               if centrality_df is not None and not centrality_df.empty else set())
    degree  = dict(G.degree())
    max_deg = max(degree.values()) if degree else 1
    node_list   = list(G.nodes)
    norm_degree = [degree.get(n, 1) / max_deg for n in node_list]

    cmap_ppi    = pal["ppi"]
    node_colors = [cmap_ppi(v * 0.85 + 0.05) for v in norm_degree]
    node_sizes  = [350 + degree.get(n, 1) * 200 for n in node_list]

    if scr:
        score_map = {}
        for _, row in ppi_df.iterrows():
            s = float(row[scr]) / 1000.0
            score_map[(row[src], row[tgt])] = s
            score_map[(row[tgt], row[src])] = s
        edge_w = [max(0.8, score_map.get(e, 0.5) * 3.0)   for e in G.edges]
        edge_a = [max(0.40, score_map.get(e, 0.5) * 0.75) for e in G.edges]
    else:
        edge_w = [1.2] * len(G.edges)
        edge_a = [0.50] * len(G.edges)

    fig, ax = plt.subplots(figsize=(14, 12), facecolor="white")

    edge_color = pal["edge"]
    for idx, (u, v) in enumerate(G.edges):
        ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                color=edge_color, linewidth=edge_w[idx],
                alpha=edge_a[idx], solid_capstyle="round", zorder=1)

    ax.scatter(
        [pos[n][0] for n in node_list],
        [pos[n][1] for n in node_list],
        s=node_sizes, color=node_colors,
        edgecolors="white", linewidths=2.0,
        zorder=3, alpha=0.93,
    )

    # Only label top hub genes
    labeled = set(sorted(node_list, key=lambda n: degree.get(n, 0), reverse=True)
                  [:min(top_n, 15)])
    cx = np.mean([p[0] for p in pos.values()])
    cy = np.mean([p[1] for p in pos.values()])
    for n in labeled:
        x, y   = pos[n]
        is_hub = n in hub_set
        fs     = 9.0 if is_hub else 7.5
        fw     = "bold" if is_hub else "normal"
        dx, dy = x - cx, y - cy
        norm_v = max(np.sqrt(dx**2 + dy**2), 1e-6)
        off    = (node_sizes[node_list.index(n)] ** 0.5) * 0.007
        lx, ly = x + dx / norm_v * off, y + dy / norm_v * off
        ax.text(lx, ly, n, ha="center", va="center",
                fontsize=fs, fontweight=fw,
                color="white" if is_hub else "#1C2833",
                zorder=5,
                path_effects=[pe.withStroke(linewidth=2.5, foreground="#1C2833")] if is_hub
                              else [pe.withStroke(linewidth=1.5, foreground="white")])

    sm = cm.ScalarMappable(cmap=cmap_ppi, norm=Normalize(0, max_deg))
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, shrink=0.42, pad=0.01, aspect=20)
    cbar.set_label("节点连接度 (Degree)", fontsize=10, color=C_GRAY)
    cbar.ax.tick_params(labelsize=9, color=C_GRAY)
    cbar.outline.set_edgecolor("#CCCCCC")

    accent = pal["accent"]
    hub_patch   = mpatches.Patch(color=cmap_ppi(0.9), label=f"核心靶点 Top {top_n}")
    other_patch = mpatches.Patch(color=cmap_ppi(0.2), label="其他靶点")
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
    pal    = _get_pal()
    C_GRAY = "#566573"

    if kegg_df is None or kegg_df.empty:
        fig, ax = plt.subplots(figsize=(10, 6), facecolor="white")
        ax.text(0.5, 0.5, "无 KEGG 数据", ha="center", va="center",
                fontsize=14, color=C_GRAY)
        ax.axis("off")
        return fig

    df = kegg_df.head(top_n).copy()
    df["neg_log_p"]  = -np.log10(df["P_value"].clip(lower=1e-300))
    df["gene_count"] = df["Genes"].apply(
        lambda x: len(str(x).split(";")) if pd.notna(x) else 1)
    df["Term_short"] = df["Term"].apply(
        lambda x: (x.split("__")[-1] if "__" in str(x) else str(x))[:50])
    df = df.sort_values("neg_log_p", ascending=True).reset_index(drop=True)

    height = max(6, len(df) * 0.65 + 2)
    fig, ax = plt.subplots(figsize=(13, height), facecolor="white")

    for i in range(len(df)):
        ax.axhspan(i - 0.48, i + 0.48,
                   color="#F7F8FA" if i % 2 == 0 else "white", zorder=0)

    cmap_kegg = pal["kegg"]
    norm  = Normalize(vmin=df["Combined_score"].min(), vmax=df["Combined_score"].max())
    sizes = (df["gene_count"] / df["gene_count"].max() * 360 + 80).values
    sc = ax.scatter(df["neg_log_p"], range(len(df)),
                    s=sizes, c=df["Combined_score"].values,
                    cmap=cmap_kegg, norm=norm,
                    edgecolors="white", linewidths=1.0,
                    zorder=3, alpha=0.92)

    for i, (_, row) in enumerate(df.iterrows()):
        ax.text(row["neg_log_p"], i, str(int(row["gene_count"])),
                ha="center", va="center", fontsize=7.5,
                color="white", fontweight="bold", zorder=4)

    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["Term_short"], fontsize=9.5)
    ax.set_ylim(-0.6, len(df) - 0.4)
    ax.set_xlabel("−log₁₀ (P value)", fontsize=12, color=C_GRAY, labelpad=8)
    ax.tick_params(axis="x", labelsize=10, color="#BBBBBB")
    ax.tick_params(axis="y", left=False)

    ref   = pal["ref_line"]
    p05   = -np.log10(0.05)
    ax.axvline(p05, color=ref, linestyle="--", linewidth=1.2, alpha=0.7)
    ax.text(p05 + 0.02, len(df) - 0.8, "P = 0.05",
            color=ref, fontsize=8.5, va="top")

    _spine_off(ax)
    ax.spines["left"].set_visible(False)

    cbar = fig.colorbar(sc, ax=ax, shrink=0.5, pad=0.02, aspect=22)
    cbar.set_label("Combined Score", fontsize=10, color=C_GRAY)
    cbar.ax.tick_params(labelsize=9)
    cbar.outline.set_edgecolor("#CCCCCC")

    max_gc     = df["gene_count"].max()
    demo_sizes = [(gc, lbl) for gc, lbl in [(5,"5 基因"),(15,"15 基因"),(30,"30 基因")]
                  if gc <= max_gc]
    for gc, lbl in demo_sizes:
        ax.scatter([], [], s=gc / max_gc * 360 + 80, c=[C_GRAY],
                   alpha=0.6, edgecolors="white", label=lbl)
    if demo_sizes:
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
    pal    = _get_pal()
    C_GRAY = "#566573"

    _cat_label = {
        "GO_BP": "生物过程 (Biological Process)",
        "GO_CC": "细胞组分 (Cellular Component)",
        "GO_MF": "分子功能 (Molecular Function)",
    }
    _cat_cmap = {
        "GO_BP": pal["go_bp"],
        "GO_CC": pal["go_cc"],
        "GO_MF": pal["go_mf"],
    }
    cat_label = _cat_label.get(category, category)
    cmap      = _cat_cmap.get(category, pal["go_bp"])

    if go_df is None or go_df.empty:
        fig, ax = plt.subplots(figsize=(10, 6), facecolor="white")
        ax.text(0.5, 0.5, f"无 {category} 数据", ha="center", va="center",
                fontsize=14, color=C_GRAY)
        ax.axis("off")
        return fig

    df = go_df.head(top_n).copy()
    df["neg_log_p"]  = -np.log10(df["P_value"].clip(lower=1e-300))
    df["gene_count"] = df["Genes"].apply(
        lambda x: len(str(x).split(";")) if pd.notna(x) else 1)
    df["Term_short"] = df["Term"].apply(
        lambda x: x.split("(GO:")[0].strip()[:55] if "(GO:" in str(x) else str(x)[:55])
    df = df.sort_values("neg_log_p", ascending=True).reset_index(drop=True)

    height = max(5.5, len(df) * 0.60 + 2)
    fig, ax = plt.subplots(figsize=(13, height), facecolor="white")

    norm   = Normalize(vmin=df["neg_log_p"].min(), vmax=df["neg_log_p"].max())
    colors = [cmap(norm(v) * 0.8 + 0.15) for v in df["neg_log_p"]]
    ax.barh(range(len(df)), df["neg_log_p"],
            color=colors, edgecolor="white", linewidth=0.6,
            height=0.62, zorder=2)

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
    ax.set_xlim(right=x_max * 1.18)

    ref = pal["ref_line"]
    p05 = -np.log10(0.05)
    ax.axvline(p05, color=ref, linestyle="--", linewidth=1.2, alpha=0.7)
    ax.text(p05 + x_max * 0.01, len(df) - 0.85, "P = 0.05",
            color=ref, fontsize=8.5, va="top")

    _spine_off(ax)
    ax.spines["left"].set_visible(False)

    ax.set_title(f"GO 富集分析 — {cat_label}  (Top {top_n})",
                 fontsize=15, fontweight="bold", color="#1C2833", pad=14)
    fig.tight_layout()
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# 5. COMPONENT-TARGET-PATHWAY NETWORK (Plotly interactive)
# ─────────────────────────────────────────────────────────────────────────────

def plot_network_plotly(
    compounds: list,
    compound_targets: dict,
    intersection_targets: list,
    top_pathways: list,
    title: str = "成分-靶点-通路网络图",
) -> go.Figure:
    pal = _get_pal()
    net = pal["net"]

    G = nx.Graph()
    node_types = {}

    for comp in compounds:
        G.add_node(comp); node_types[comp] = "compound"
        for t in compound_targets.get(comp, []):
            if t in intersection_targets:
                G.add_node(t); node_types[t] = "target"
                G.add_edge(comp, t)

    for path in top_pathways[:12]:
        short = path[:45]
        G.add_node(short); node_types[short] = "pathway"
        for t in intersection_targets[:10]:
            if t in G.nodes:
                G.add_edge(t, short)

    if len(G.nodes) == 0:
        return go.Figure()

    try:
        for n, t in node_types.items():
            G.nodes[n]["layer"] = {"compound": 0, "target": 1, "pathway": 2}[t]
        pos = nx.multipartite_layout(G, subset_key="layer", scale=2.5, align="vertical")
    except Exception:
        pos = nx.kamada_kawai_layout(G) if len(G.nodes) <= 80 \
              else nx.spring_layout(G, k=3.5, seed=42)

    style = {
        "compound": dict(color=net["compound"], size=24, symbol="circle"),
        "target":   dict(color=net["target"],   size=18, symbol="circle"),
        "pathway":  dict(color=net["pathway"],  size=22, symbol="diamond"),
    }
    label_cn = {"compound":"活性成分","target":"交集靶点","pathway":"KEGG通路"}

    edge_x, edge_y = [], []
    for u, v in G.edges:
        edge_x += [pos[u][0], pos[v][0], None]
        edge_y += [pos[u][1], pos[v][1], None]

    traces = [go.Scatter(
        x=edge_x, y=edge_y, mode="lines",
        line=dict(width=1.0, color="rgba(120,130,145,0.40)"),
        hoverinfo="none", showlegend=False,
    )]

    for ntype, st in style.items():
        nodes = [n for n, t in node_types.items() if t == ntype]
        if not nodes:
            continue
        deg   = [G.degree(n) for n in nodes]
        max_d = max(deg) if deg else 1
        sizes = [st["size"] + d / max_d * st["size"] * 0.9 for d in deg]
        hover = [f"<b>{n}</b><br>类型: {label_cn[ntype]}<br>连接数: {G.degree(n)}"
                 for n in nodes]
        traces.append(go.Scatter(
            x=[pos[n][0] for n in nodes],
            y=[pos[n][1] for n in nodes],
            mode="markers+text",
            marker=dict(size=sizes, color=st["color"], symbol=st["symbol"],
                        opacity=0.90, line=dict(width=1.5, color="white")),
            text=[n[:32] for n in nodes],
            textposition="top center",
            textfont=dict(size=9 if ntype == "pathway" else 8, color="#1C2833"),
            hovertemplate="%{customdata}<extra></extra>",
            customdata=hover,
            name=label_cn[ntype],
        ))

    fig = go.Figure(data=traces, layout=go.Layout(
        title=dict(text=f"<b>{title}</b>",
                   font=dict(size=17, color="#1C2833", family="Arial"),
                   x=0.5, xanchor="center"),
        showlegend=True,
        legend=dict(orientation="v", x=1.01, y=0.98,
                    bgcolor="rgba(255,255,255,0.9)",
                    bordercolor="#DDDDDD", borderwidth=1,
                    font=dict(size=12)),
        hovermode="closest",
        xaxis=dict(showgrid=False, zeroline=False, showticklabels=False, visible=False),
        yaxis=dict(showgrid=False, zeroline=False, showticklabels=False, visible=False),
        paper_bgcolor="white", plot_bgcolor="#F8F9FA",
        height=750, margin=dict(l=20, r=170, t=70, b=20),
    ))
    return fig


# ─────────────────────────────────────────────────────────────────────────────
# Utility
# ─────────────────────────────────────────────────────────────────────────────

def apply_cfg(cfg: dict):
    """Apply user font/palette/style config to matplotlib rcParams globally."""
    global _PAL
    if not cfg:
        return

    # Palette
    theme = cfg.get("color_theme", "NPG")
    _PAL  = PALETTES.get(theme, PALETTES["NPG"])

    # Font
    font_style = cfg.get("font_style", "黑体")
    if font_style == "宋体" and _songti_avail:
        plt.rcParams["font.family"] = _songti_avail + ["DejaVu Sans"]
    else:
        plt.rcParams["font.family"] = (_heiti_avail or _any_cjk) + ["DejaVu Sans"]

    fs = cfg.get("font_scale", 1.0)
    plt.rcParams.update({
        "font.size":          10 * fs,
        "axes.titlesize":     15 * fs,
        "axes.labelsize":     12 * fs,
        "xtick.labelsize":    10 * fs,
        "ytick.labelsize":    10 * fs,
        "legend.fontsize":    10 * fs,
        "figure.dpi":         cfg.get("dpi", 150),
        "savefig.dpi":        cfg.get("dpi", 180),
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
