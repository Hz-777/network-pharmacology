"""
网络药理学一键分析平台
流程: 中药输入 → TCMSP筛选 → PubChem SMILES → Swiss靶点 → 疾病靶点 → 交集 → PPI → GO/KEGG → 报告
"""

import streamlit as st
import pandas as pd
import numpy as np
import os, time, io, yaml
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime, timezone, timedelta
from pathlib import Path

_BJT = timezone(timedelta(hours=8))

def _now():
    return datetime.now(_BJT)

st.set_page_config(
    page_title="网络药理学一键分析平台",
    page_icon="🌿",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ── 登录认证 ──────────────────────────────────────────────────────────────────

_USERS_FILE = Path(__file__).parent / "users.yaml"

def _secrets_to_dict(obj) -> dict:
    """递归把 st.secrets 的 AttrDict 转成普通 dict。"""
    if hasattr(obj, "items"):
        return {k: _secrets_to_dict(v) for k, v in obj.items()}
    return obj

def _load_auth_config() -> dict:
    # 优先读 users.yaml（本地 + 云端均适用，提交到 git 后云端自动同步）
    if _USERS_FILE.exists():
        with open(_USERS_FILE, encoding="utf-8") as f:
            return yaml.safe_load(f)
    # 降级：从 Streamlit Secrets 读（仅当 users.yaml 不存在时）
    if "credentials" in st.secrets:
        return _secrets_to_dict(st.secrets)
    raise FileNotFoundError("users.yaml")

try:
    import streamlit_authenticator as stauth

    _auth_cfg = _load_auth_config()
    _authenticator = stauth.Authenticate(
        _auth_cfg["credentials"],
        _auth_cfg["cookie"]["name"],
        _auth_cfg["cookie"]["key"],
        _auth_cfg["cookie"]["expiry_days"],
    )

    _authenticator.login()
    _auth_status = st.session_state.get("authentication_status")

    if _auth_status is False:
        st.error("❌ 用户名或密码错误")
        st.stop()
    elif _auth_status is None:
        st.markdown("""
        <div style="max-width:420px;margin:6rem auto 0;text-align:center;">
          <h2 style="color:#1E5631;">🌿 网络药理学分析平台</h2>
          <p style="color:#666;">请登录以继续使用</p>
        </div>
        """, unsafe_allow_html=True)
        st.stop()

    # 已登录 — 侧边栏显示用户信息和退出按钮
    _current_user = st.session_state.get("name", "")
    _current_role = _auth_cfg["credentials"]["usernames"].get(
        st.session_state.get("username", ""), {}
    ).get("role", "user")

except FileNotFoundError:
    st.error("⚠️ 未找到 users.yaml，请先创建用户配置文件")
    st.code("cp users.yaml.example users.yaml  # 或运行 python3 tools/add_user.py")
    st.stop()
except Exception as e:
    st.error(f"认证模块加载失败: {e}")
    st.stop()

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
        "report_path": None, "smiles_df": None, "docking_scores": None,
    }.items():
        if k not in st.session_state:
            st.session_state[k] = v

_init()


def add_log(msg, level="info"):
    icon = {"info":"ℹ️","ok":"✅","warn":"⚠️","error":"❌"}.get(level,"•")
    ts = _now().strftime("%H:%M:%S")
    st.session_state.log.append(f"[{ts}] {icon} {msg}")


# ── Sidebar ───────────────────────────────────────────────────────────────────

with st.sidebar:
    # 用户信息 & 退出
    _role_label = {"admin": "管理员 👑", "user": "用户"}.get(_current_role, _current_role)
    st.markdown(f"**👤 {_current_user}** · {_role_label}")
    _authenticator.logout("退出登录", location="sidebar")
    st.markdown("---")
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
    enable_docking = st.checkbox(
        "💊 启用分子对接（可选）",
        value=False,
        help="勾选后分析完成可在「分子对接」Tab 中对关键靶蛋白进行虚拟筛选。\n"
             "对接依赖 CB-Dock2 在线服务，耗时较长，默认关闭。",
    )

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

    with st.expander("🗑️ 缓存管理", expanded=False):
        from modules.cache import cache_clear, cache_info
        _ci = cache_info()
        if _ci["enabled"]:
            st.caption(f"📦 已缓存 **{_ci['items']}** 条 · 占用 **{_ci['size_mb']} MB**")
            st.caption("命中缓存 = 跳过 API 调用，重复分析几乎秒出")
        else:
            st.caption("⚠️ 磁盘缓存未启用（diskcache 未安装）")
        if st.button("🧹 清除全部缓存", use_container_width=True):
            cache_clear()
            st.cache_data.clear()
            st.success("✅ 缓存已清除，下次分析将重新调用所有 API")

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

    # ── 1. TCMSP（多草药并行）+ 疾病靶点同步预取 ─────────────────────────────
    step(1, 9, f"并行获取: TCMSP活性成分（{len(herb_names)}味药）+ 疾病靶点预取")

    # 疾病靶点与草药查询完全独立 —— 立刻放进后台线程，和步骤 1-3 同时跑
    _disease_executor = ThreadPoolExecutor(max_workers=1)
    _disease_future   = _disease_executor.submit(_cached_disease_targets, disease)
    add_log(f"  已在后台启动疾病靶点查询: [{disease}]")

    # 多草药并行查询 TCMSP
    def _fetch_herb(herb: str):
        df = _cached_herb(herb, ob_th, dl_th)
        if df.empty:
            raise ValueError("空结果")
        df["Herb"] = herb
        return herb, df

    all_comp = []
    with ThreadPoolExecutor(max_workers=min(len(herb_names), 4)) as executor:
        herb_futures = {executor.submit(_fetch_herb, h): h for h in herb_names}
        for future in as_completed(herb_futures):
            herb = herb_futures[future]
            try:
                _, df = future.result()
                all_comp.append(df)
                add_log(f"  {herb}: {len(df)} 个活性成分", "ok")
            except Exception as e:
                add_log(f"  {herb} TCMSP查询失败({e})，使用内置数据", "warn")
                all_comp.append(_demo_compounds(herb))

    compounds_df = pd.concat(all_comp, ignore_index=True)
    name_col = "mol_name" if "mol_name" in compounds_df.columns else compounds_df.columns[0]
    compounds_df = compounds_df.drop_duplicates(subset=[name_col]).reset_index(drop=True)
    st.session_state.compounds_df = compounds_df
    add_log(f"活性成分合计: {len(compounds_df)} 个", "ok")

    # ── 2. PubChem SMILES ─────────────────────────────────────────────────────
    step(2, 9, "从 PubChem 并行获取 SMILES 结构式")
    comp_names = compounds_df[name_col].dropna().unique().tolist()
    batch = comp_names[:20]
    smiles_rows = [None] * len(batch)

    with ThreadPoolExecutor(max_workers=5) as executor:
        future_map = {executor.submit(_cached_compound_info, nm): i for i, nm in enumerate(batch)}
        done = 0
        for future in as_completed(future_map):
            idx = future_map[future]
            try:
                smiles_rows[idx] = future.result()
            except Exception:
                smiles_rows[idx] = {"name": batch[idx], "CID": None, "SMILES": None, "formula": None, "MW": None}
            done += 1
            if done % 5 == 0 or done == len(batch):
                add_log(f"  PubChem {done}/{len(batch)} 个化合物已处理")

    smiles_rows = [r for r in smiles_rows if r is not None]
    smiles_df = pd.DataFrame(smiles_rows)
    valid_smiles = smiles_df[smiles_df["SMILES"].notna() & (smiles_df["SMILES"] != "")]
    add_log(f"获得 SMILES: {len(valid_smiles)}/{len(batch)} 个化合物", "ok")
    st.session_state.smiles_df = valid_smiles

    # ── 3. 多数据库靶点查询（ChEMBL + HERB 并行）────────────────────────────
    step(3, 9, "ChEMBL + HERB 双库并行查询药物靶点")
    target_rows, cmap, all_drug_genes = [], {}, set()

    query_compounds = list(compounds_df[name_col].dropna().unique())[:15]
    from modules.herb_targets import get_herb_targets

    def _fetch_one_target(nm: str):
        sm = ""
        if not valid_smiles.empty and "name" in valid_smiles.columns:
            matched = valid_smiles[valid_smiles["name"] == nm]
            sm = matched["SMILES"].iloc[0] if not matched.empty else ""
        chembl_df = _cached_predict_targets(sm, nm)
        herb_df   = get_herb_targets(nm)
        return nm, chembl_df, herb_df

    with ThreadPoolExecutor(max_workers=3) as executor:
        futures = {executor.submit(_fetch_one_target, nm): nm for nm in query_compounds}
        fetch_results = {}
        for future in as_completed(futures):
            try:
                nm, chembl_df, herb_df = future.result()
                fetch_results[nm] = (chembl_df, herb_df)
            except Exception:
                fetch_results[futures[future]] = (pd.DataFrame(), pd.DataFrame())

    for nm in query_compounds:
        chembl_df, herb_df = fetch_results.get(nm, (pd.DataFrame(), pd.DataFrame()))

        genes_chembl = chembl_df["Gene"].dropna().unique().tolist() if not chembl_df.empty else []
        genes_herb   = herb_df["Gene"].dropna().unique().tolist()   if not herb_df.empty   else []
        genes = sorted(set(genes_chembl) | set(genes_herb))

        if genes:
            cmap[nm] = genes
            all_drug_genes.update(genes)
            for g in genes_chembl:
                target_rows.append({"Compound": nm, "Gene": g, "Source": "ChEMBL"})
            for g in genes_herb:
                if g not in set(genes_chembl):
                    target_rows.append({"Compound": nm, "Gene": g, "Source": "HERB"})
            src_info = []
            if genes_chembl: src_info.append(f"ChEMBL {len(genes_chembl)}")
            if genes_herb:   src_info.append(f"HERB {len(genes_herb)}")
            add_log(f"  {nm}: {len(genes)} 个靶点（{'，'.join(src_info)}）", "ok")
        else:
            add_log(f"  {nm}: 两库均未收录", "warn")

    if not target_rows:
        add_log("双库均无结果，使用内置靶点数据", "warn")
        drug_targets_df, cmap, all_drug_genes = _demo_drug_targets(herb_names)
        demo_comp_names = list(DEMO_COMPOUND_TARGETS.keys())
        if compounds_df[name_col].isin(demo_comp_names).sum() == 0:
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
    add_log(f"药物靶点合计: {len(all_drug_genes)} 个唯一基因（ChEMBL + HERB）", "ok")

    # ── 4. 疾病靶点（收取步骤 1 启动的后台结果）─────────────────────────────
    step(4, 9, f"收取疾病靶点结果（{disease}）")
    try:
        disease_df = _disease_future.result(timeout=60)   # 此时大概率已完成
        if disease_df.empty:
            raise ValueError("空结果")
        add_log(f"  Open Targets 返回 {len(disease_df)} 个靶点", "ok")
    except Exception as e:
        add_log(f"  在线库无返回({e})，使用内置疾病靶点", "warn")
        disease_df = _demo_disease_targets(disease)
    finally:
        _disease_executor.shutdown(wait=False)
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
    ts = _now().strftime("%Y%m%d_%H%M%S")
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
        plot_network_plotly, ppi_to_pyvis_html, network_to_pyvis_html,
        fig_to_base64, fig_to_bytes, apply_cfg,
    )

    # Retrieve saved cfg (falls back to current sidebar values)
    _cfg = st.session_state.get("fig_cfg", fig_cfg)
    apply_cfg(_cfg)
    _fmt = _cfg.get("fmt", "png")
    _dpi = _cfg.get("dpi", 180)
    _mime = {"png": "image/png", "pdf": "application/pdf", "svg": "image/svg+xml"}[_fmt]

    tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs([
        "📊 活性成分", "🎯 靶点 & 韦恩图", "🕸️ PPI 网络",
        "🔬 GO/KEGG 富集", "🌐 网络图", "📥 下载报告", "💊 分子对接",
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
            st.markdown("---")
            st.markdown("#### 🖱️ 可拖动 PPI 互作探索图")
            st.caption("节点可自由拖动；稳定后物理引擎自动关闭，拖动不再弹回")
            ppi_html = ppi_to_pyvis_html(
                ppi_df,
                cen_df if cen_df is not None else pd.DataFrame(),
                top_n=30,
                height=600,
            )
            import streamlit.components.v1 as components
            components.html(ppi_html, height=620, scrolling=False)

    # ── Tab 4 ─────────────────────────────────────────────────────────────────
    with tab4:
        st.markdown("### GO / KEGG 富集分析")
        enrichment = st.session_state.enrichment
        if enrichment:
            pval_thresh = st.slider(
                "P-value 阈值（实时筛选）", 0.001, 0.05, 0.05,
                step=0.005, format="%.3f",
                help="拖动滑块动态过滤显著性阈值，图表自动重绘",
            )
            t_kegg, t_bp, t_cc, t_mf = st.tabs(["KEGG 通路","GO-BP","GO-CC","GO-MF"])
            for tab_obj, cat in [(t_kegg,"KEGG"),(t_bp,"GO_BP"),(t_cc,"GO_CC"),(t_mf,"GO_MF")]:
                with tab_obj:
                    edf_raw = enrichment.get(cat)
                    if edf_raw is not None and not edf_raw.empty:
                        edf = edf_raw[edf_raw["P_value"] <= pval_thresh].copy()
                        st.caption(f"显示 {len(edf)} / {len(edf_raw)} 条（P≤{pval_thresh:.3f}）")
                        st.dataframe(edf, use_container_width=True, height=300)
                        if not edf.empty:
                            if cat == "KEGG":
                                efig = plot_kegg_bubbles(edf, title="KEGG 通路富集")
                            else:
                                efig = plot_go_barplot(edf, category=cat)
                            st.image(f"data:image/png;base64,{fig_to_base64(efig, dpi=_dpi)}", use_column_width=True)
                            efig2 = plot_kegg_bubbles(edf, title="KEGG 通路富集") if cat == "KEGG" \
                                    else plot_go_barplot(edf, category=cat)
                            st.download_button(
                                f"⬇ 下载 {cat} 图 (.{_fmt})",
                                fig_to_bytes(efig2, fmt=_fmt, dpi=_dpi),
                                file_name=f"{cat.lower()}.{_fmt}", mime=_mime,
                            )
                        else:
                            st.info(f"当前阈值（P≤{pval_thresh:.3f}）下无显著 {cat} 结果")
                    else:
                        st.info(f"暂无 {cat} 数据（P≤0.05）")
        else:
            st.info("分析完成后显示富集结果")

    # ── Tab 5 ─────────────────────────────────────────────────────────────────
    with tab5:
        st.markdown("### 成分-靶点-通路 交互网络图")
        if hasattr(st.session_state, "network_fig") and st.session_state.network_fig is not None:
            st.plotly_chart(st.session_state.network_fig, use_container_width=True)
            st.markdown("---")
            st.markdown("#### 🖱️ 可拖动同心圆探索图")
            st.caption("节点初始位置为三圈同心布局（成分·靶点·通路），可自由拖动重排")
            # Reconstruct data for pyvis from session state
            _ss         = st.session_state
            _ct_map     = getattr(_ss, "compound_target_map", {}) or {}
            _comp_list  = list(_ct_map.keys())[:15]
            _inter_list = list(_ss.intersection_genes) if _ss.intersection_genes else []
            _pathways: list = []
            if _ss.enrichment and isinstance(_ss.enrichment, dict):
                _kdf = _ss.enrichment.get("KEGG")
                if _kdf is not None and not _kdf.empty:
                    _pathways = [
                        (t.split("__")[-1] if "__" in t else t)
                        for t in _kdf.head(14)["Term"].tolist()
                    ]
            net_html = network_to_pyvis_html(
                _comp_list, _ct_map, _inter_list, _pathways, height=680,
            )
            import streamlit.components.v1 as components
            components.html(net_html, height=700, scrolling=False)
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

    # ── Tab 7: 分子对接 ───────────────────────────────────────────────────────
    with tab7:
        st.markdown("### 💊 分子对接（Molecular Docking）")

        if not enable_docking:
            st.info(
                "**分子对接功能未启用。**\n\n"
                "如需使用，请在左侧边栏勾选 **「💊 启用分子对接（可选）」**，"
                "完成网药分析后再来此 Tab 进行虚拟筛选。\n\n"
                "网药分析（步骤 1~9）可独立完成，无需启用此功能。"
            )
        else:
            st.markdown(
                "基于 **CB-Dock2** 在线服务，对活性化合物与关键靶蛋白进行自动对接，"
                "输出结合能热图。结合能越负（绝对值越大）表示结合亲和力越强。"
            )

            from modules.docking import (
                fetch_pdb_batch, run_docking_matrix,
                plot_docking_heatmap, check_cbdock2_status,
            )

            # ── 服务状态检查 ──────────────────────────────────────────────────
            with st.expander("🔗 CB-Dock2 服务状态", expanded=False):
                if st.button("检测 CB-Dock2 是否可用"):
                    with st.spinner("检测中..."):
                        status = check_cbdock2_status()
                    if status["available"]:
                        st.success(f"✅ CB-Dock2 可用（{status['label']}）: {status['url']}")
                    else:
                        st.error("❌ CB-Dock2 暂不可用（主站和备用IP均无响应）。对接功能需等待服务恢复。")

            # ── 参数选择 ──────────────────────────────────────────────────────
            _smiles_df    = st.session_state.get("smiles_df")
            _compounds_df = st.session_state.get("compounds_df")
            _inter_genes  = st.session_state.get("intersection_genes", [])
            _cent_df      = st.session_state.get("centrality_df")

            if _compounds_df is None:
                st.info("请先点击「开始一键分析」完成分析流程后再进行分子对接。")
            else:
                # 化合物名称来自 compounds_df，SMILES 在对接时按需实时查询
                _name_col = "mol_name" if "mol_name" in _compounds_df.columns else _compounds_df.columns[0]

                # 构建 avail_comps：优先已有 SMILES，其余按名称列出（对接时再查）
                _smiles_map = {}
                if _smiles_df is not None and not _smiles_df.empty and "SMILES" in _smiles_df.columns:
                    for _, r in _smiles_df.iterrows():
                        if r.get("SMILES"):
                            _smiles_map[r.get("name", "")] = r["SMILES"]
                if "SMILES" in _compounds_df.columns:
                    for _, r in _compounds_df.iterrows():
                        nm = r.get(_name_col, "")
                        if nm and r.get("SMILES") and nm not in _smiles_map:
                            _smiles_map[nm] = r["SMILES"]

                all_comp_names = _compounds_df[_name_col].dropna().unique().tolist()
                avail_comps = pd.DataFrame({
                    "name": all_comp_names,
                    "SMILES": [_smiles_map.get(n, "") for n in all_comp_names],
                    "has_smiles": [bool(_smiles_map.get(n)) for n in all_comp_names],
                })

                col_d1, col_d2 = st.columns(2)

                with col_d1:
                    st.markdown("**① 选择化合物**")
                    # 优先推荐已有 SMILES 的化合物
                    has_smiles = avail_comps[avail_comps["has_smiles"]]["name"].tolist()
                    no_smiles  = avail_comps[~avail_comps["has_smiles"]]["name"].tolist()
                    comp_options = has_smiles + no_smiles
                    default_comps = has_smiles[:3] if has_smiles else comp_options[:3]
                    if no_smiles and not has_smiles:
                        st.caption("⚠️ 以下化合物 SMILES 将在对接时从 PubChem 实时获取")
                    elif no_smiles:
                        st.caption(f"✅ {len(has_smiles)} 个已有 SMILES，{len(no_smiles)} 个将实时查询")
                    selected_comps = st.multiselect(
                        "选择参与对接的化合物（建议 ≤5 个）",
                        options=comp_options,
                        default=default_comps,
                        key="dock_comp_sel",
                    )

                with col_d2:
                    st.markdown("**② 选择靶蛋白**")
                    if _cent_df is not None and not _cent_df.empty:
                        suggested_genes = _cent_df.head(10)["Gene"].tolist()
                    elif _inter_genes:
                        suggested_genes = list(_inter_genes)[:10]
                    else:
                        suggested_genes = []

                    target_mode = st.radio(
                        "靶蛋白来源",
                        ["从 RCSB 自动获取（按 Hub 基因）", "手动上传 PDB 文件"],
                        key="dock_target_mode",
                    )

                    if target_mode == "从 RCSB 自动获取（按 Hub 基因）":
                        selected_genes = st.multiselect(
                            "选择靶基因（将自动从 RCSB PDB 下载结构）",
                            options=suggested_genes + [g for g in _inter_genes if g not in suggested_genes],
                            default=suggested_genes[:3] if len(suggested_genes) >= 3 else suggested_genes,
                            key="dock_gene_sel",
                        )
                        uploaded_pdbs = {}
                    else:
                        uploaded_files = st.file_uploader(
                            "上传 PDB 文件（可多选）",
                            type=["pdb"],
                            accept_multiple_files=True,
                            key="dock_pdb_upload",
                        )
                        uploaded_pdbs = {}
                        selected_genes = []
                        for uf in (uploaded_files or []):
                            gene_label = uf.name.replace(".pdb", "")
                            uploaded_pdbs[gene_label] = uf.read().decode("utf-8", errors="replace")
                            selected_genes.append(gene_label)

                st.markdown("---")

                if st.button("🚀 开始分子对接", type="primary",
                             disabled=not selected_comps or not selected_genes,
                             key="dock_run_btn"):
                    dock_log = st.empty()
                    dock_status = []

                    # Thread-safe logger: only appends to list (no Streamlit calls in threads)
                    def _bg_cb(msg):
                        dock_status.append(msg)

                    # UI refresh helper: only called from main thread
                    def _flush_log():
                        dock_log.text("\n".join(dock_status[-10:]))

                    comp_inputs = []
                    from modules.pubchem import get_compound_info as _gci
                    for cn in selected_comps:
                        row = avail_comps[avail_comps["name"] == cn]
                        smiles = row["SMILES"].iloc[0] if not row.empty else ""
                        if not smiles:
                            _bg_cb(f"PubChem 查询 SMILES: {cn} ...")
                            _flush_log()
                            try:
                                info = _gci(cn)
                                smiles = info.get("SMILES", "") if info else ""
                            except Exception:
                                smiles = ""
                        if smiles:
                            comp_inputs.append({"name": cn, "SMILES": smiles})
                            _bg_cb(f"✅ {cn}: SMILES 已获取")
                        else:
                            _bg_cb(f"⚠️ {cn}: 未找到 SMILES，跳过")
                    _flush_log()

                    _bg_cb("正在从 RCSB 获取靶蛋白结构...")
                    _flush_log()
                    if target_mode == "从 RCSB 自动获取（按 Hub 基因）":
                        with st.spinner(f"下载 {len(selected_genes)} 个蛋白结构..."):
                            # Pass thread-safe callback (no Streamlit calls inside)
                            pdb_results = fetch_pdb_batch(
                                selected_genes, max_workers=3, progress_callback=_bg_cb
                            )
                        _flush_log()
                    else:
                        pdb_results = {g: ("uploaded", c) for g, c in uploaded_pdbs.items()}

                    target_inputs = []
                    for gene in selected_genes:
                        result = pdb_results.get(gene)
                        if result:
                            pdb_id, pdb_content = result
                            target_inputs.append({"gene": gene, "pdb_id": pdb_id, "pdb_content": pdb_content})
                            _bg_cb(f"✅ {gene}: PDB {pdb_id}")
                        else:
                            _bg_cb(f"⚠️ {gene}: 未找到 PDB 结构，跳过")
                    _flush_log()

                    if not target_inputs:
                        st.error("未能获取任何靶蛋白结构，请检查网络连接或手动上传 PDB 文件。")
                    elif not comp_inputs:
                        st.error("未找到所选化合物的 SMILES，请检查网络连接后重试。")
                    else:
                        _bg_cb(f"\n开始对接：{len(comp_inputs)} 个化合物 × {len(target_inputs)} 个蛋白...")
                        _flush_log()
                        with st.spinner("对接计算中，请耐心等待..."):
                            try:
                                scores_df = run_docking_matrix(
                                    comp_inputs, target_inputs,
                                    progress_callback=_bg_cb,  # thread-safe
                                )
                                st.session_state.docking_scores = scores_df
                            except Exception as e:
                                err_str = str(e)
                                if "timed out" in err_str.lower() or "timeout" in err_str.lower() or "ConnectionError" in type(e).__name__:
                                    st.error(
                                        "CB-Dock2 服务连接超时。可能原因：\n"
                                        "1. 服务器当前负载较高，请稍后重试\n"
                                        "2. 网络防火墙阻断了对外部服务的访问\n"
                                        "3. CB-Dock2 主站/备用站均暂时不可用\n\n"
                                        f"技术详情：{err_str}"
                                    )
                                else:
                                    st.error(f"对接失败: {err_str}")
                                scores_df = pd.DataFrame()
                        _flush_log()
                        dock_log.empty()

                # ── 显示结果 ──────────────────────────────────────────────────
                _scores = st.session_state.get("docking_scores")
                if _scores is not None and not _scores.empty:
                    st.markdown("#### 对接结果热图")
                    st.caption("结合能（kcal/mol）：数值越负代表结合越强")
                    heatmap_png = plot_docking_heatmap(_scores)
                    st.image(heatmap_png, use_container_width=True)

                    st.markdown("#### 原始数据")
                    st.dataframe(_scores.style.format("{:.2f}"), use_container_width=True)

                    buf_d = io.BytesIO()
                    _scores.to_excel(buf_d, index=True)
                    st.download_button(
                        "📥 下载对接评分表（Excel）",
                        buf_d.getvalue(),
                        file_name="分子对接评分.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                    )
                elif _scores is not None and _scores.empty:
                    st.warning(
                        "对接结果为空。可能原因：\n"
                        "- CB-Dock2 服务暂时不可用\n"
                        "- 网络超时（对接通常需要 1~3 分钟）\n"
                        "- PDB 结构未能成功获取\n\n"
                        "建议：稍后重试，或手动上传 PDB 文件。"
                    )

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
