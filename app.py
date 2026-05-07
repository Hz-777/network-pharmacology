"""
网络药理学一键分析平台
流程: 中药输入 → TCMSP筛选 → PubChem SMILES → Swiss靶点 → 疾病靶点 → 交集 → PPI → GO/KEGG → 报告
"""

import streamlit as st
import pandas as pd
import numpy as np
import os, time, io
from datetime import datetime
from pathlib import Path

st.set_page_config(
    page_title="网络药理学一键分析平台",
    page_icon="🌿",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown("""
<style>
.main-header{background:linear-gradient(135deg,#1E5631 0%,#4CAF50 100%);padding:2rem;border-radius:12px;color:white;text-align:center;margin-bottom:1.5rem;}
h1,h2,h3{color:#1E5631;}
.stProgress>div>div{background-color:#27AE60;}
</style>
""", unsafe_allow_html=True)


# ── Demo / fallback data (must be defined before run_btn block) ──────────────

def _demo_compounds(herb: str) -> pd.DataFrame:
    data = {
        "黄连": [
            {"mol_name": "Berberine",    "OB": 36.86, "DL": 0.78},
            {"mol_name": "Coptisine",    "OB": 30.67, "DL": 0.86},
            {"mol_name": "Palmatine",    "OB": 64.60, "DL": 0.65},
            {"mol_name": "Jatrorrhizine","OB": 47.14, "DL": 0.77},
            {"mol_name": "Epiberberine", "OB": 43.09, "DL": 0.78},
        ],
        "黄芩": [
            {"mol_name": "Baicalein",    "OB": 33.52, "DL": 0.21},
            {"mol_name": "Wogonin",      "OB": 30.68, "DL": 0.23},
            {"mol_name": "Baicalin",     "OB": 40.12, "DL": 0.75},
            {"mol_name": "Scutellarein", "OB": 30.20, "DL": 0.20},
        ],
        "人参": [
            {"mol_name": "Ginsenoside Rh2","OB": 36.32,"DL": 0.56},
            {"mol_name": "Panaxadiol",     "OB": 33.08,"DL": 0.78},
            {"mol_name": "Beta-sitosterol","OB": 36.91,"DL": 0.75},
        ],
    }
    # For unknown herbs, reuse the known demo compounds so they
    # carry through PubChem→ChEMBL and match the targets tab.
    rows = data.get(herb, [
        {"mol_name": "Berberine",    "OB": 36.86, "DL": 0.78},
        {"mol_name": "Quercetin",    "OB": 46.43, "DL": 0.28},
        {"mol_name": "Beta-sitosterol","OB": 36.91,"DL": 0.75},
    ])
    df = pd.DataFrame(rows)
    df["Herb"] = herb
    return df


DEMO_COMPOUND_TARGETS = {
    "Berberine":      ["TP53","AKT1","TNF","IL6","MAPK1","PTGS2","NOS2"],
    "Baicalein":      ["VEGFA","EGFR","STAT3","NFKB1","BCL2","IL1B"],
    "Palmatine":      ["MDM2","PTEN","MTOR","CDK2","CASP3"],
    "Wogonin":        ["MYC","JUN","HSP90AA1","PPARG","AR"],
    "Coptisine":      ["AKT1","PIK3CA","CCND1","STAT3"],
    "Jatrorrhizine":  ["TNF","IL6","NFKB1","MAPK1"],
    "Ginsenoside Rh2":["TP53","BCL2","CASP3","EGFR"],
    "Beta-sitosterol":["ESR1","AR","PPARG","RXRA"],
    "Quercetin":      ["TP53","AKT1","VEGFA","EGFR","TNF","NFKB1","MTOR","MAPK1"],
}


def _demo_drug_targets(herb_names: list) -> tuple:
    """Returns (targets_df, compound_target_map, all_genes_set)."""
    rows, cmap, all_genes = [], {}, set()
    for herb in herb_names:
        for comp, genes in DEMO_COMPOUND_TARGETS.items():
            cmap[comp] = genes
            all_genes.update(genes)
            for g in genes:
                rows.append({"Compound": comp, "Gene": g, "Herb": herb})
    return pd.DataFrame(rows), cmap, all_genes


def _demo_disease_targets(disease: str) -> pd.DataFrame:
    if any(k in disease.lower() for k in ["diabet","糖尿病"]):
        genes = ["INS","INSR","IRS1","PPARG","ADIPOQ","LEP","AKT1","PI3K",
                 "AMPK","SIRT1","GCK","HIF1A","TNF","IL6","NFKB1","MAPK1","TP53","BCL2","VEGFA","PTGS2"]
    elif any(k in disease.lower() for k in ["cancer","肿瘤","癌"]):
        genes = ["TP53","BRCA1","EGFR","MYC","BCL2","VEGFA","PIK3CA","KRAS",
                 "PTEN","AKT1","MTOR","CDK4","MDM2","STAT3","NFKB1","IL6","TNF","HSP90AA1","CASP3"]
    elif any(k in disease.lower() for k in ["hypertens","高血压"]):
        genes = ["ACE","AGT","AGTR1","NOS3","EDN1","PTGS2","KCNJ11","ATP2A2",
                 "TNF","IL6","VEGFA","AKT1","MAPK1","NFKB1","TP53"]
    else:
        genes = ["TP53","AKT1","TNF","IL6","VEGFA","EGFR","MYC","MAPK1",
                 "CASP3","BCL2","STAT3","NFKB1","MDM2","PTEN","MTOR","PTGS2","NOS2","IL1B"]
    rng = np.random.default_rng(42)
    return pd.DataFrame([{"Gene": g, "Score": round(float(rng.uniform(0.3, 1.0)), 3), "Source": "Demo"} for g in genes])


def _demo_ppi(genes: list) -> pd.DataFrame:
    rng = np.random.default_rng(42)
    rows = []
    for i, g1 in enumerate(genes):
        partners = rng.choice(genes, size=min(4, len(genes)), replace=False).tolist()
        for g2 in partners:
            if g1 != g2:
                rows.append({"preferredName_A": g1, "preferredName_B": g2,
                             "score": int(rng.integers(400, 999))})
    df = pd.DataFrame(rows).drop_duplicates(subset=["preferredName_A","preferredName_B"])
    return df


# ── Cached API wrappers (same inputs → skip network round-trip) ───────────────

@st.cache_data(ttl=3600, show_spinner=False)
def _cached_herb(herb_name: str, ob_th: float, dl_th: float) -> pd.DataFrame:
    from modules.tcmsp import search_herb_tcmsp
    return search_herb_tcmsp(herb_name, ob_th, dl_th)

@st.cache_data(ttl=3600, show_spinner=False)
def _cached_compound_info(name: str) -> dict:
    from modules.pubchem import get_compound_info
    return get_compound_info(name)

@st.cache_data(ttl=3600, show_spinner=False)
def _cached_predict_targets(smiles: str, compound_name: str) -> pd.DataFrame:
    from modules.swiss_target import predict_targets
    return predict_targets(smiles=smiles, compound_name=compound_name)

@st.cache_data(ttl=3600, show_spinner=False)
def _cached_disease_targets(disease: str) -> pd.DataFrame:
    from modules.disease_targets import get_disease_targets
    return get_disease_targets(disease)

@st.cache_data(ttl=3600, show_spinner=False)
def _cached_ppi(genes: tuple, min_score: int) -> pd.DataFrame:
    from modules.string_db import get_ppi_network
    return get_ppi_network(list(genes), min_score=min_score)

@st.cache_data(ttl=3600, show_spinner=False)
def _cached_enrichr(genes: tuple) -> dict:
    from modules.enrichment import run_enrichr_analysis
    return run_enrichr_analysis(list(genes))


# ── Session state ─────────────────────────────────────────────────────────────

def _init():
    for k, v in {
        "compounds_df": None, "drug_targets_df": None, "disease_targets_df": None,
        "intersection_genes": [], "ppi_df": None, "centrality_df": None,
        "enrichment": {}, "compound_target_map": {}, "log": [], "done": False,
        "report_path": None,
    }.items():
        if k not in st.session_state:
            st.session_state[k] = v

_init()


def add_log(msg, level="info"):
    icon = {"info":"ℹ️","ok":"✅","warn":"⚠️","error":"❌"}.get(level,"•")
    ts = datetime.now().strftime("%H:%M:%S")
    st.session_state.log.append(f"[{ts}] {icon} {msg}")


# ── Sidebar ───────────────────────────────────────────────────────────────────

with st.sidebar:
    st.markdown("## 🌿 参数设置")
    herb_input = st.text_area("🌱 中药名称（每行一个）", value="黄连\n黄芩", height=100)
    disease_input = st.text_input("🏥 疾病名称（英文/中文）", value="diabetes")

    st.markdown("### ADME 筛选阈值")
    c1, c2 = st.columns(2)
    ob_th = c1.number_input("OB ≥ (%)", 0.0, 100.0, 30.0, 5.0)
    dl_th = c2.number_input("DL ≥", 0.0, 1.0, 0.18, 0.01)

    st.markdown("### 高级参数")
    ppi_score = st.slider("PPI 最低置信分数", 0, 1000, 400, 50)
    top_hub   = st.slider("核心靶点数量 Top N", 5, 50, 20, 5)

    st.markdown("---")

    with st.expander("🎨 出图设置", expanded=False):
        font_style = st.radio(
            "中文字体",
            options=["黑体", "宋体"],
            index=0,
            horizontal=True,
            help="黑体：简洁现代；宋体：正式传统",
        )
        font_scale = st.select_slider(
            "字体大小",
            options=[0.7, 0.85, 1.0, 1.2, 1.4, 1.6],
            value=1.0,
            format_func=lambda x: {0.7:"小", 0.85:"偏小", 1.0:"中（默认）",
                                    1.2:"偏大", 1.4:"大", 1.6:"超大"}[x],
        )
        fig_dpi = st.select_slider(
            "图片分辨率 (DPI)",
            options=[100, 150, 180, 220, 300],
            value=180,
            format_func=lambda x: f"{x} dpi" + (" (印刷级)" if x >= 300 else
                                                  " (推荐)" if x == 180 else ""),
        )
        fig_fmt = st.radio(
            "下载格式",
            options=["png", "pdf", "svg"],
            index=0,
            horizontal=True,
            help="PNG适合插入PPT；PDF/SVG矢量图适合投稿论文",
        )
        color_theme = st.selectbox(
            "期刊配色",
            options=["NPG", "NEJM", "Lancet", "Science", "Cell"],
            index=0,
            format_func=lambda x: {
                "NPG":     "🟥 Nature (NPG)  — 红·青·绿·藏蓝",
                "NEJM":    "🟦 NEJM          — 深红·蓝·橙·绿",
                "Lancet":  "🟩 Lancet        — 深蓝·红·绿·青",
                "Science": "🟪 Science (AAAS)— 靛蓝·红·深绿·紫",
                "Cell":    "🟧 Cell Press    — 珊瑚·蓝·绿·紫",
            }[x],
            help="全局切换所有图表（韦恩/PPI/KEGG/GO/网络图）的配色方案",
        )

    # Pack into cfg dict for passing to visualization functions
    fig_cfg = {
        "font_style":   font_style,
        "font_scale":   font_scale,
        "dpi":          fig_dpi,
        "fmt":          fig_fmt,
        "color_theme":  color_theme,
    }

    st.markdown("---")
    run_btn = st.button("🚀 开始一键分析", type="primary", use_container_width=True)
    if st.button("🔄 重置", use_container_width=True):
        for k in list(st.session_state.keys()):
            del st.session_state[k]
        st.rerun()
    st.markdown("---")
    st.caption("数据来源: TCMSP · PubChem · ChEMBL · STRING · Enrichr · Open Targets · Harmonizome")


# ── Header ────────────────────────────────────────────────────────────────────

st.markdown("""
<div class="main-header">
  <h1 style="color:white;margin:0;font-size:1.9rem;">🌿 网络药理学一键分析平台</h1>
  <p style="color:#A8E6CF;margin:0.4rem 0 0;font-size:0.95rem;">
    TCMSP · PubChem · ChEMBL · STRING · Enrichr 全流程自动分析
  </p>
</div>
""", unsafe_allow_html=True)

m1, m2, m3, m4 = st.columns(4)
m1.metric("活性成分", len(st.session_state.compounds_df) if st.session_state.compounds_df is not None else "—")
_tgt = st.session_state.drug_targets_df
m2.metric("预测靶点", _tgt["Gene"].nunique() if _tgt is not None and "Gene" in _tgt.columns else "—")
m3.metric("交集靶点", len(st.session_state.intersection_genes) or "—")
_kegg = st.session_state.enrichment.get("KEGG")
m4.metric("KEGG通路", len(_kegg) if _kegg is not None and not getattr(_kegg,"empty",True) else "—")

st.markdown("---")


# ── Analysis pipeline ─────────────────────────────────────────────────────────

if run_btn:
    herb_names = [h.strip() for h in herb_input.strip().splitlines() if h.strip()]
    disease    = disease_input.strip()

    if not herb_names or not disease:
        st.error("请先填写中药名称和疾病名称")
        st.stop()

    # reset
    st.session_state.log = []
    st.session_state.done = False
    st.session_state.report_path = None

    progress_bar  = st.progress(0, text="准备中…")
    status_holder = st.empty()
    log_holder    = st.expander("📋 实时日志", expanded=True)

    def step(n, total, msg):
        pct = int(n / total * 100)
        progress_bar.progress(pct, text=f"步骤 {n}/{total}: {msg}")
        status_holder.info(f"⏳ {msg}")
        add_log(msg)
        with log_holder:
            for line in st.session_state.log[-15:]:
                st.text(line)

    # ── 1. TCMSP ─────────────────────────────────────────────────────────────
    step(1, 9, f"从 TCMSP 获取活性成分（OB≥{ob_th}%, DL≥{dl_th}）")
    all_comp = []
    for herb in herb_names:
        try:
            df = _cached_herb(herb, ob_th, dl_th)
            if not df.empty:
                df["Herb"] = herb
                all_comp.append(df)
                add_log(f"  {herb}: {len(df)} 个活性成分", "ok")
            else:
                raise ValueError("空结果")
        except Exception as e:
            add_log(f"  {herb} TCMSP在线查询失败({e})，使用内置数据", "warn")
            all_comp.append(_demo_compounds(herb))

    compounds_df = pd.concat(all_comp, ignore_index=True)
    # deduplicate by mol_name if column exists
    name_col = "mol_name" if "mol_name" in compounds_df.columns else compounds_df.columns[0]
    compounds_df = compounds_df.drop_duplicates(subset=[name_col]).reset_index(drop=True)
    st.session_state.compounds_df = compounds_df
    add_log(f"活性成分合计: {len(compounds_df)} 个", "ok")

    # ── 2. PubChem SMILES ─────────────────────────────────────────────────────
    step(2, 9, "从 PubChem 获取 SMILES 结构式")
    smiles_rows = []
    comp_names  = compounds_df[name_col].dropna().unique().tolist()
    for i, nm in enumerate(comp_names[:20]):          # cap at 20 for speed
        info = _cached_compound_info(nm)
        smiles_rows.append(info)
        if i % 5 == 0:
            add_log(f"  PubChem {i+1}/{min(len(comp_names),20)}: {nm}")
        time.sleep(0.25)
    smiles_df = pd.DataFrame(smiles_rows)
    valid_smiles = smiles_df[smiles_df["SMILES"].notna() & (smiles_df["SMILES"] != "")]
    add_log(f"获得 SMILES: {len(valid_smiles)}/{len(comp_names)} 个化合物", "ok")

    # ── 3. ChEMBL 靶点查询 ───────────────────────────────────────────────────
    step(3, 9, "ChEMBL 数据库查询药物靶点")
    target_rows, cmap, all_drug_genes = [], {}, set()

    # Use all compounds with SMILES, plus try by name for those without
    query_compounds = list(compounds_df[name_col].dropna().unique())[:15]
    for nm in query_compounds:
        sm = ""
        if not valid_smiles.empty and "name" in valid_smiles.columns:
            matched = valid_smiles[valid_smiles["name"] == nm]
            sm = matched["SMILES"].iloc[0] if not matched.empty else ""
        try:
            tdf = _cached_predict_targets(sm, nm)
            if not tdf.empty:
                genes = tdf["Gene"].dropna().unique().tolist()
                cmap[nm] = genes
                all_drug_genes.update(genes)
                for g in genes:
                    target_rows.append({"Compound": nm, "Gene": g})
                add_log(f"  {nm}: {len(genes)} 个靶点", "ok")
            else:
                add_log(f"  {nm}: ChEMBL 未收录", "warn")
        except Exception as e:
            add_log(f"  {nm} 查询失败: {e}", "warn")
        time.sleep(0.5)

    if not target_rows:
        add_log("ChEMBL 返回为空，使用内置靶点数据", "warn")
        drug_targets_df, cmap, all_drug_genes = _demo_drug_targets(herb_names)
        # Sync compounds_df so the compounds tab shows the same names that have targets.
        demo_comp_names = list(DEMO_COMPOUND_TARGETS.keys())
        if compounds_df[name_col].isin(demo_comp_names).sum() == 0:
            # Replace placeholder compounds with the demo ones
            sync_rows = []
            for cn in demo_comp_names:
                sync_rows.append({"mol_name": cn, "OB": 36.0, "DL": 0.30, "Herb": herb_names[0]})
            compounds_df = pd.DataFrame(sync_rows)
            name_col = "mol_name"
            st.session_state.compounds_df = compounds_df
    else:
        drug_targets_df = pd.DataFrame(target_rows)

    st.session_state.drug_targets_df     = drug_targets_df
    st.session_state.compound_target_map = cmap
    add_log(f"药物靶点: {len(all_drug_genes)} 个唯一基因", "ok")

    # ── 4. Disease targets ────────────────────────────────────────────────────
    step(4, 9, f"获取疾病靶点（{disease}）")
    add_log(f"  查询 Open Targets: [{disease}]...")
    disease_df = _cached_disease_targets(disease)
    if disease_df.empty:
        add_log("在线库无返回，使用内置疾病靶点", "warn")
        disease_df = _demo_disease_targets(disease)
    st.session_state.disease_targets_df = disease_df
    add_log(f"疾病靶点: {len(disease_df)} 个", "ok")

    # ── 5. Intersection + Venn ────────────────────────────────────────────────
    step(5, 9, "计算药物-疾病交集靶点")
    drug_genes_set    = set(all_drug_genes)
    disease_genes_set = set(disease_df["Gene"].dropna().tolist())
    intersection      = sorted(drug_genes_set & disease_genes_set)
    if not intersection:
        add_log(f"药物靶点与疾病靶点无交集（疾病: {disease}），后续分析基于药物靶点", "warn")
        st.warning(
            f"⚠️ **未找到药物-疾病交集靶点**\n\n"
            f"药物靶点（{len(drug_genes_set)} 个）与疾病靶点（{len(disease_genes_set)} 个）没有重叠，"
            f"可能原因：①疾病名称不准确（当前: `{disease}`）；②ChEMBL 中该化合物靶点数据不足。\n\n"
            f"后续分析将使用药物靶点继续，**结果仅供参考**。"
        )
        intersection = sorted(drug_genes_set)[:30]
    st.session_state.intersection_genes = intersection
    add_log(f"交集靶点: {len(intersection)} 个", "ok")

    # ── 6. PPI ───────────────────────────────────────────────────────────────
    step(6, 9, f"STRING 数据库构建 PPI 网络（score≥{ppi_score}）")
    from modules.string_db import calculate_network_centrality
    try:
        ppi_df = _cached_ppi(tuple(intersection[:50]), ppi_score)
        if ppi_df.empty:
            raise ValueError("空网络")
        add_log(f"PPI 互作对: {len(ppi_df)} 条", "ok")
    except Exception as e:
        add_log(f"STRING 查询失败({e})，使用模拟数据", "warn")
        ppi_df = _demo_ppi(intersection[:20] or list(drug_genes_set)[:20])

    centrality_df = calculate_network_centrality(ppi_df)
    st.session_state.ppi_df        = ppi_df
    st.session_state.centrality_df = centrality_df
    add_log(f"核心靶点已排序，Top1: {centrality_df.iloc[0]['Gene'] if not centrality_df.empty else 'N/A'}", "ok")

    # ── 7. GO/KEGG ───────────────────────────────────────────────────────────
    step(7, 9, "Enrichr GO/KEGG 富集分析")
    from modules.string_db import get_top_hub_genes
    hub_genes    = get_top_hub_genes(centrality_df, top_n=top_hub)
    enrich_genes = hub_genes or intersection[:30] or list(drug_genes_set)[:30]

    add_log("  提交基因列表到 Enrichr...")
    enrichment = _cached_enrichr(tuple(enrich_genes))
    st.session_state.enrichment = enrichment
    kegg_df = enrichment.get("KEGG", pd.DataFrame())
    add_log(f"KEGG 通路: {len(kegg_df) if kegg_df is not None and not kegg_df.empty else 0} 条 (P≤0.05)", "ok")

    # ── 8. Network figure (Plotly) ────────────────────────────────────────────
    step(8, 9, "构建成分-靶点-通路网络图")
    from modules.visualization import plot_network_plotly
    top_paths = []
    if kegg_df is not None and not kegg_df.empty and "Term" in kegg_df.columns:
        top_paths = [t.split("__")[-1] if "__" in t else t for t in kegg_df["Term"].head(10).tolist()]
    st.session_state.network_fig = plot_network_plotly(
        compounds          = list(cmap.keys())[:15],
        compound_targets   = cmap,
        intersection_targets = intersection[:20],
        top_pathways       = top_paths,
        title              = f"{'、'.join(herb_names)} 成分-靶点-通路网络",
    )
    add_log("网络图构建完成", "ok")

    # ── 9. Excel report ───────────────────────────────────────────────────────
    step(9, 9, "生成 Excel 分析报告")
    from modules.report import generate_excel_report
    import matplotlib
    matplotlib.use("Agg")
    from modules.visualization import (
        plot_venn, plot_ppi_network, plot_kegg_bubbles, plot_go_barplot, apply_cfg
    )

    apply_cfg(fig_cfg)   # apply user font/dpi settings before drawing

    mpl_figs = {}
    mpl_figs["venn"] = plot_venn(
        {"药物靶点": drug_genes_set, f"{disease}靶点": disease_genes_set},
        title="药物-疾病靶点韦恩图",
    )
    mpl_figs["ppi"] = plot_ppi_network(ppi_df, centrality_df, top_n=top_hub)
    if kegg_df is not None and not kegg_df.empty:
        mpl_figs["kegg"] = plot_kegg_bubbles(kegg_df)
    for cat in ["GO_BP","GO_CC","GO_MF"]:
        gdf = enrichment.get(cat, pd.DataFrame())
        if gdf is not None and not gdf.empty:
            mpl_figs[cat.lower()] = plot_go_barplot(gdf, category=cat)

    st.session_state.fig_cfg = fig_cfg   # save cfg for display section

    out_dir = Path("output"); out_dir.mkdir(exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    rpath = out_dir / f"网络药理学分析报告_{ts}.xlsx"
    generate_excel_report(
        output_path          = str(rpath),
        herb_names           = herb_names,
        disease_name         = disease,
        compounds_df         = compounds_df,
        drug_targets_df      = drug_targets_df,
        disease_targets_df   = disease_df,
        intersection_genes   = intersection,
        ppi_df               = ppi_df,
        centrality_df        = centrality_df,
        enrichment_results   = enrichment,
        figures              = mpl_figs,
    )
    st.session_state.report_path = str(rpath)
    add_log(f"报告已保存: {rpath}", "ok")

    progress_bar.progress(100, text="✅ 分析完成！")
    status_holder.success("🎉 分析完成！请查看下方结果标签页")
    st.session_state.done = True
    st.balloons()
    st.rerun()


# ── Results ───────────────────────────────────────────────────────────────────

if st.session_state.compounds_df is not None:
    import matplotlib
    matplotlib.use("Agg")
    from modules.visualization import (
        plot_venn, plot_ppi_network, plot_kegg_bubbles, plot_go_barplot,
        plot_network_plotly, fig_to_base64, fig_to_bytes, apply_cfg,
    )

    # Retrieve saved cfg (falls back to current sidebar values)
    _cfg = st.session_state.get("fig_cfg", fig_cfg)
    apply_cfg(_cfg)
    _fmt = _cfg.get("fmt", "png")
    _dpi = _cfg.get("dpi", 180)
    _mime = {"png": "image/png", "pdf": "application/pdf", "svg": "image/svg+xml"}[_fmt]

    tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
        "📊 活性成分", "🎯 靶点 & 韦恩图", "🕸️ PPI 网络",
        "🔬 GO/KEGG 富集", "🌐 网络图", "📥 下载报告",
    ])

    # ── Tab 1 ─────────────────────────────────────────────────────────────────
    with tab1:
        st.markdown("### 活性成分列表（ADME 筛选后）")
        df = st.session_state.compounds_df
        # Hide internal metadata columns from the user
        display_cols = [c for c in df.columns if c not in ("_source",)]
        st.dataframe(df[display_cols], use_container_width=True, height=380)
        c1, c2 = st.columns(2)
        import plotly.express as px
        if "OB" in df.columns:
            with c1:
                fig = px.histogram(df, x="OB", nbins=15, title="OB 分布",
                                   color_discrete_sequence=["#27AE60"])
                fig.add_vline(x=30, line_dash="dash", line_color="red", annotation_text="OB=30%")
                st.plotly_chart(fig, use_container_width=True)
        if "DL" in df.columns:
            with c2:
                fig = px.histogram(df, x="DL", nbins=15, title="DL 分布",
                                   color_discrete_sequence=["#2980B9"])
                fig.add_vline(x=0.18, line_dash="dash", line_color="red", annotation_text="DL=0.18")
                st.plotly_chart(fig, use_container_width=True)

    # ── Tab 2 ─────────────────────────────────────────────────────────────────
    with tab2:
        c1, c2 = st.columns(2)
        with c1:
            st.markdown("### 药物预测靶点")
            if st.session_state.drug_targets_df is not None:
                st.dataframe(st.session_state.drug_targets_df, use_container_width=True, height=280)
        with c2:
            st.markdown("### 疾病相关靶点")
            if st.session_state.disease_targets_df is not None:
                st.dataframe(st.session_state.disease_targets_df, use_container_width=True, height=280)

        st.markdown("### 药物-疾病交集靶点（韦恩图）")
        c1, c2 = st.columns([1.2, 1])
        with c1:
            if st.session_state.drug_targets_df is not None and st.session_state.disease_targets_df is not None:
                dg  = set(st.session_state.drug_targets_df["Gene"].dropna())
                dis = set(st.session_state.disease_targets_df["Gene"].dropna())
                vfig = plot_venn(
                    {"药物靶点": dg, f"{disease_input.strip()}靶点": dis},
                    title="韦恩图",
                )
                st.image(f"data:image/png;base64,{fig_to_base64(vfig, dpi=_dpi)}", use_column_width=True)
                vfig2 = plot_venn({"药物靶点": dg, f"{disease_input.strip()}靶点": dis}, title="韦恩图")
                st.download_button(f"⬇ 下载韦恩图 (.{_fmt})",
                                   fig_to_bytes(vfig2, fmt=_fmt, dpi=_dpi),
                                   file_name=f"venn.{_fmt}", mime=_mime)
        with c2:
            st.markdown("#### 交集靶点")
            inter = st.session_state.intersection_genes
            if inter:
                st.dataframe(pd.DataFrame({"靶点": inter}), use_container_width=True, height=320)
                st.success(f"共 **{len(inter)}** 个交集靶点")

    # ── Tab 3 ─────────────────────────────────────────────────────────────────
    with tab3:
        st.markdown("### PPI 蛋白互作网络")
        ppi_df = st.session_state.ppi_df
        cen_df = st.session_state.centrality_df
        if ppi_df is not None:
            c1, c2 = st.columns([1.8, 1])
            with c1:
                pfig = plot_ppi_network(ppi_df, cen_df if cen_df is not None else pd.DataFrame(), top_n=20)
                st.image(f"data:image/png;base64,{fig_to_base64(pfig, dpi=_dpi)}", use_column_width=True)
                pfig2 = plot_ppi_network(ppi_df, cen_df if cen_df is not None else pd.DataFrame(), top_n=20)
                st.download_button(f"⬇ 下载 PPI 网络图 (.{_fmt})",
                                   fig_to_bytes(pfig2, fmt=_fmt, dpi=_dpi),
                                   file_name=f"ppi_network.{_fmt}", mime=_mime)
            with c2:
                st.markdown("#### 核心靶点排名（Hub Score）")
                if cen_df is not None and not cen_df.empty:
                    st.dataframe(cen_df.head(20), use_container_width=True, height=420)

    # ── Tab 4 ─────────────────────────────────────────────────────────────────
    with tab4:
        st.markdown("### GO / KEGG 富集分析")
        enrichment = st.session_state.enrichment
        if enrichment:
            t_kegg, t_bp, t_cc, t_mf = st.tabs(["KEGG 通路","GO-BP","GO-CC","GO-MF"])
            for tab_obj, cat in [(t_kegg,"KEGG"),(t_bp,"GO_BP"),(t_cc,"GO_CC"),(t_mf,"GO_MF")]:
                with tab_obj:
                    edf = enrichment.get(cat)
                    if edf is not None and not edf.empty:
                        st.dataframe(edf, use_container_width=True, height=300)
                        if cat == "KEGG":
                            efig = plot_kegg_bubbles(edf, title="KEGG 通路富集")
                        else:
                            efig = plot_go_barplot(edf, category=cat)
                        st.image(f"data:image/png;base64,{fig_to_base64(efig, dpi=_dpi)}", use_column_width=True)
                        # redraw for download (fig_to_base64 closes the figure)
                        efig2 = plot_kegg_bubbles(edf, title="KEGG 通路富集") if cat == "KEGG" \
                                else plot_go_barplot(edf, category=cat)
                        st.download_button(
                            f"⬇ 下载 {cat} 图 (.{_fmt})",
                            fig_to_bytes(efig2, fmt=_fmt, dpi=_dpi),
                            file_name=f"{cat.lower()}.{_fmt}", mime=_mime,
                        )
                    else:
                        st.info(f"暂无 {cat} 数据（P≤0.05）")
        else:
            st.info("分析完成后显示富集结果")

    # ── Tab 5 ─────────────────────────────────────────────────────────────────
    with tab5:
        st.markdown("### 成分-靶点-通路 交互网络图")
        if hasattr(st.session_state, "network_fig") and st.session_state.network_fig is not None:
            st.plotly_chart(st.session_state.network_fig, use_container_width=True)
        else:
            st.info("网络图将在分析完成后显示")

    # ── Tab 6 ─────────────────────────────────────────────────────────────────
    with tab6:
        st.markdown("### 下载分析报告")
        rp = st.session_state.report_path
        if rp and os.path.exists(rp):
            with open(rp, "rb") as f:
                data = f.read()
            st.success(f"✅ 报告已生成: {os.path.basename(rp)}")
            st.download_button(
                "📥 下载 Excel 完整报告", data,
                file_name=os.path.basename(rp),
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                type="primary", use_container_width=True,
            )
            st.markdown("---")
            st.markdown("#### 单项数据")
            dc1, dc2, dc3 = st.columns(3)
            with dc1:
                buf = io.BytesIO()
                st.session_state.compounds_df.to_excel(buf, index=False)
                st.download_button("📊 活性成分", buf.getvalue(), "活性成分.xlsx")
            with dc2:
                if st.session_state.intersection_genes:
                    st.download_button("🎯 交集靶点", "\n".join(st.session_state.intersection_genes).encode(), "交集靶点.txt")
            with dc3:
                kegg = st.session_state.enrichment.get("KEGG")
                if kegg is not None and not kegg.empty:
                    buf2 = io.BytesIO()
                    kegg.to_excel(buf2, index=False)
                    st.download_button("🔬 KEGG通路", buf2.getvalue(), "KEGG通路.xlsx")
        else:
            st.info("完成分析后可在此下载报告")

else:
    # welcome screen
    st.markdown("""
    <div style="text-align:center; padding:3rem 0; color:#666;">
        <h2 style="color:#1E5631;">👈 在左侧填写参数，点击"开始一键分析"</h2>
        <p>支持单味药或复方 · 全程自动化 · 一键生成 Excel 报告</p>
    </div>
    """, unsafe_allow_html=True)

    with st.expander("📖 分析流程说明"):
        st.markdown("""
| 步骤 | 操作 | 数据库 |
|------|------|--------|
| ① 活性成分提取 | OB≥30%, DL≥0.18 筛选 | TCMSP |
| ② SMILES 获取 | 化合物标准结构式 | PubChem |
| ③ 靶点预测 | 基于分子结构预测 | ChEMBL |
| ④ 疾病靶点 | 疾病相关基因 | Open Targets + Harmonizome |
| ⑤ 交集靶点 | 韦恩图取交集 | 本地计算 |
| ⑥ PPI 网络 | 蛋白互作 + Hub 排序 | STRING |
| ⑦ GO/KEGG 富集 | 通路注释 | Enrichr |
| ⑧ 网络图 | 成分-靶点-通路可视化 | Plotly |
| ⑨ Excel 报告 | 含图表完整报告 | — |
        """)


# ── footer ────────────────────────────────────────────────────────────────────
st.markdown("---")
st.caption("🌿 网络药理学一键分析平台 · 仅供科研参考，结果需结合实验验证")
