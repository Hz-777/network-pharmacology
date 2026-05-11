"""
Disease target retrieval.
Sources (queried in parallel, results merged):
  1. Open Targets Platform  — GraphQL API, EFO/MONDO IDs, scored 0–1
  2. Harmonizome/DisGeNET   — REST API, literature evidence, free/no key
  3. UniProt                — REST API, Swiss-Prot disease annotations, free/no key
  4. NCBI Gene              — E-utilities, disease/phenotype curated associations, free/no key
Supports both English and Chinese disease names.
"""

import re
import requests
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from modules.cache import cache_get, cache_set, make_key

OT_API          = "https://api.platform.opentargets.org/api/v4/graphql"
HARMONIZOME_API = "https://maayanlab.cloud/Harmonizome/api/1.0"
OT_HEADERS      = {"Content-Type": "application/json"}
HEADERS         = {"User-Agent": "Mozilla/5.0 (network-pharmacology-app/1.0)"}

# Chinese → English mapping for common diseases in TCM research
_CN_TO_EN = {
    # 心血管
    "高血压": "hypertension",
    "动脉高压": "hypertension",
    "心肌梗死": "myocardial infarction",
    "心肌梗塞": "myocardial infarction",
    "冠心病": "coronary artery disease",
    "冠状动脉疾病": "coronary artery disease",
    "心力衰竭": "heart failure",
    "心衰": "heart failure",
    "心房颤动": "atrial fibrillation",
    "房颤": "atrial fibrillation",
    "动脉粥样硬化": "atherosclerosis",
    "心绞痛": "angina pectoris",
    "脑卒中": "stroke",
    "中风": "stroke",
    # 代谢
    "糖尿病": "diabetes mellitus",
    "2型糖尿病": "type 2 diabetes mellitus",
    "1型糖尿病": "type 1 diabetes mellitus",
    "肥胖": "obesity",
    "肥胖症": "obesity",
    "高脂血症": "hyperlipidemia",
    "血脂异常": "dyslipidemia",
    "痛风": "gout",
    "非酒精性脂肪肝": "non-alcoholic fatty liver disease",
    "脂肪肝": "non-alcoholic fatty liver disease",
    # 肿瘤
    "肿瘤": "cancer",
    "癌症": "cancer",
    "肝癌": "liver cancer",
    "肝细胞癌": "hepatocellular carcinoma",
    "肺癌": "lung cancer",
    "非小细胞肺癌": "non-small cell lung carcinoma",
    "乳腺癌": "breast cancer",
    "胃癌": "gastric cancer",
    "结肠癌": "colon cancer",
    "结直肠癌": "colorectal cancer",
    "宫颈癌": "cervical cancer",
    "前列腺癌": "prostate cancer",
    "卵巢癌": "ovarian cancer",
    "胰腺癌": "pancreatic cancer",
    "食管癌": "esophageal cancer",
    "鼻咽癌": "nasopharyngeal carcinoma",
    "甲状腺癌": "thyroid cancer",
    "白血病": "leukemia",
    "淋巴瘤": "lymphoma",
    # 神经/精神
    "阿尔茨海默病": "Alzheimer disease",
    "阿尔茨海默": "Alzheimer disease",
    "老年痴呆": "Alzheimer disease",
    "帕金森病": "Parkinson disease",
    "帕金森": "Parkinson disease",
    "抑郁症": "major depressive disorder",
    "抑郁": "major depressive disorder",
    "焦虑症": "generalized anxiety disorder",
    "焦虑": "generalized anxiety disorder",
    "精神分裂症": "schizophrenia",
    "癫痫": "epilepsy",
    "脑梗死": "cerebral infarction",
    "失眠": "insomnia",
    "偏头痛": "migraine",
    # 炎症/免疫
    "类风湿关节炎": "rheumatoid arthritis",
    "风湿性关节炎": "rheumatoid arthritis",
    "骨关节炎": "osteoarthritis",
    "系统性红斑狼疮": "systemic lupus erythematosus",
    "狼疮": "systemic lupus erythematosus",
    "克罗恩病": "Crohn disease",
    "溃疡性结肠炎": "ulcerative colitis",
    "炎症性肠病": "inflammatory bowel disease",
    "哮喘": "asthma",
    "慢性阻塞性肺病": "chronic obstructive pulmonary disease",
    "慢阻肺": "chronic obstructive pulmonary disease",
    "银屑病": "psoriasis",
    "牛皮癣": "psoriasis",
    # 肝/肾/消化
    "肝炎": "hepatitis",
    "肝硬化": "liver cirrhosis",
    "肾病综合征": "nephrotic syndrome",
    "慢性肾病": "chronic kidney disease",
    "肾炎": "nephritis",
    "胃炎": "gastritis",
    "消化性溃疡": "peptic ulcer",
    "胃溃疡": "gastric ulcer",
    # 感染
    "新冠": "COVID-19",
    "新型冠状病毒": "COVID-19",
    "新冠肺炎": "COVID-19",
    "肺炎": "pneumonia",
    "脓毒症": "sepsis",
    # 骨骼
    "骨质疏松": "osteoporosis",
    "骨质疏松症": "osteoporosis",
    # 内分泌
    "甲状腺功能亢进": "hyperthyroidism",
    "甲亢": "hyperthyroidism",
    "甲状腺功能减退": "hypothyroidism",
    "多囊卵巢综合征": "polycystic ovary syndrome",
    "多囊卵巢": "polycystic ovary syndrome",
}


def _to_english(disease: str) -> str:
    """Convert Chinese disease name to English. Returns input unchanged if already English."""
    # Check if input contains Chinese characters
    if not re.search(r'[一-鿿]', disease):
        return disease
    # Exact match first
    if disease in _CN_TO_EN:
        return _CN_TO_EN[disease]
    # Partial match (longest matching key wins)
    best_key, best_en = "", ""
    for cn, en in _CN_TO_EN.items():
        if cn in disease and len(cn) > len(best_key):
            best_key, best_en = cn, en
    return best_en if best_en else disease


def _search_disease_id(disease: str) -> object:
    """Search disease name → return best EFO/MONDO ID."""
    query = """
    query($term: String!) {
      search(queryString: $term, entityNames: ["disease"], page: {index:0, size:5}) {
        hits { id name entity }
      }
    }"""
    try:
        r = requests.post(
            OT_API,
            json={"query": query, "variables": {"term": disease}},
            headers=OT_HEADERS,
            timeout=15,
        )
        hits = r.json().get("data", {}).get("search", {}).get("hits", [])
        if hits:
            return hits[0]["id"]
    except Exception:
        pass
    return None


def _get_associated_targets(efo_id: str, size: int = 100) -> pd.DataFrame:
    """Fetch associated gene targets for a disease EFO/MONDO ID."""
    query = """
    query($efoId: String!, $size: Int!) {
      disease(efoId: $efoId) {
        name
        associatedTargets(page: {index: 0, size: $size}) {
          rows {
            target { approvedSymbol approvedName }
            score
          }
        }
      }
    }"""
    try:
        r = requests.post(
            OT_API,
            json={"query": query, "variables": {"efoId": efo_id, "size": size}},
            headers=OT_HEADERS,
            timeout=20,
        )
        data = r.json().get("data", {})
        disease_data = data.get("disease") or {}
        rows = (disease_data.get("associatedTargets") or {}).get("rows", [])
        if rows:
            return pd.DataFrame([
                {
                    "Gene":   row["target"]["approvedSymbol"],
                    "Score":  round(row["score"], 4),
                    "Source": "OpenTargets",
                }
                for row in rows
                if row.get("target", {}).get("approvedSymbol")
            ])
    except Exception:
        pass
    return pd.DataFrame()


def _query_harmonizome(disease: str) -> pd.DataFrame:
    """
    Fetch disease-gene associations from Harmonizome (DisGeNET dataset).
    Returns DataFrame with columns [Gene, Score, Source].
    Score is fixed at 0.3 (binary literature evidence, no quantitative ranking).
    """
    dataset = "DisGeNET+Gene-Disease+Associations"
    try:
        encoded = requests.utils.quote(disease, safe="")
        r = requests.get(
            f"{HARMONIZOME_API}/gene_set/{encoded}/{dataset}",
            headers=HEADERS,
            timeout=20,
        )
        if r.status_code != 200:
            return pd.DataFrame()
        assocs = r.json().get("associations", [])
        if not assocs:
            return pd.DataFrame()
        rows = [
            {"Gene": a["gene"]["symbol"], "Score": 0.3, "Source": "Harmonizome"}
            for a in assocs
            if (a.get("gene") or {}).get("symbol")
        ]
        return pd.DataFrame(rows) if rows else pd.DataFrame()
    except Exception:
        return pd.DataFrame()


def _query_uniprot(disease: str) -> pd.DataFrame:
    """
    Fetch human proteins annotated with a disease from UniProt REST API.
    Uses cc_disease field (Comments → Disease/Phenotype Descriptions).
    Returns DataFrame with columns [Gene, Score, Source].
    Score: 0.5 for Swiss-Prot reviewed, 0.3 for TrEMBL unreviewed.
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    # Query reviewed (Swiss-Prot) first for quality; fall back to include TrEMBL
    params = {
        "query": f"cc_disease:{disease} AND organism_id:9606 AND reviewed:true",
        "fields": "gene_names,reviewed",
        "format": "json",
        "size": 200,
    }
    try:
        resp = requests.get(url, params=params, headers=HEADERS, timeout=20)
        results = resp.json().get("results", []) if resp.ok else []

        # If very few reviewed hits, also fetch unreviewed
        if len(results) < 20:
            params2 = {**params, "query": f"cc_disease:{disease} AND organism_id:9606",
                       "size": 300}
            resp2 = requests.get(url, params=params2, headers=HEADERS, timeout=20)
            results = resp2.json().get("results", []) if resp2.ok else results

        rows = []
        seen = set()
        for entry in results:
            reviewed = entry.get("entryType", "") == "UniProtKB reviewed (Swiss-Prot)"
            score = 0.5 if reviewed else 0.3
            for gene_group in entry.get("genes", []):
                gene_name = (gene_group.get("geneName") or {}).get("value", "")
                if gene_name and gene_name not in seen:
                    seen.add(gene_name)
                    rows.append({"Gene": gene_name, "Score": score, "Source": "UniProt"})
                    break
        return pd.DataFrame(rows) if rows else pd.DataFrame()
    except Exception:
        return pd.DataFrame()


NCBI_EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _query_ncbi_gene(disease: str) -> pd.DataFrame:
    """
    Fetch human genes associated with a disease via NCBI Gene E-utilities.
    Uses the disease/phenotype MeSH annotation index — same source as OMIM/ClinVar.
    Returns DataFrame with columns [Gene, Score, Source].
    Score fixed at 0.4 (curated NCBI association, no quantitative score available).
    """
    try:
        # Search Gene db: disease/phenotype field + human filter
        search_resp = requests.get(
            f"{NCBI_EUTILS}/esearch.fcgi",
            params={
                "db": "gene",
                "term": f"{disease}[disease/phenotype] AND Homo sapiens[Organism]",
                "retmax": 300,
                "retmode": "json",
            },
            headers=HEADERS,
            timeout=15,
        )
        if not search_resp.ok:
            return pd.DataFrame()
        gene_ids = search_resp.json().get("esearchresult", {}).get("idlist", [])
        if not gene_ids:
            return pd.DataFrame()

        # Fetch gene symbols in batches of 100
        rows = []
        for i in range(0, len(gene_ids), 100):
            batch = gene_ids[i : i + 100]
            summary_resp = requests.get(
                f"{NCBI_EUTILS}/esummary.fcgi",
                params={"db": "gene", "id": ",".join(batch), "retmode": "json"},
                headers=HEADERS,
                timeout=20,
            )
            if not summary_resp.ok:
                continue
            result = summary_resp.json().get("result", {})
            for gid in batch:
                info = result.get(gid, {})
                symbol = info.get("name", "")
                if symbol and symbol != "":
                    rows.append({"Gene": symbol, "Score": 0.4, "Source": "NCBI Gene"})
        return pd.DataFrame(rows) if rows else pd.DataFrame()
    except Exception:
        return pd.DataFrame()


def get_disease_targets(disease: str, min_score: float = 0.0,
                        progress_callback=None) -> pd.DataFrame:
    """
    Retrieve disease-associated gene targets from 4 sources in parallel:
      Open Targets, Harmonizome/DisGeNET, UniProt, NCBI Gene.
    Accepts English or Chinese disease names.
    Deduplication priority: Open Targets > UniProt > NCBI Gene > Harmonizome.
    """
    key = make_key("disease_targets", disease, min_score)
    cached = cache_get(key)
    if cached is not None:
        return cached

    en_disease = _to_english(disease)
    if en_disease != disease and progress_callback:
        progress_callback(f"中文疾病名转换: {disease} → {en_disease}")

    def _fetch_ot():
        efo_id = _search_disease_id(en_disease)
        if efo_id:
            return _get_associated_targets(efo_id, size=200)
        return pd.DataFrame()

    tasks = {
        "Open Targets": _fetch_ot,
        "UniProt":      lambda: _query_uniprot(en_disease),
        "NCBI Gene":    lambda: _query_ncbi_gene(en_disease),
        "Harmonizome":  lambda: _query_harmonizome(en_disease),
    }

    if progress_callback:
        progress_callback("并行查询 4 个数据库：Open Targets / UniProt / NCBI Gene / Harmonizome ...")

    source_dfs: dict = {}
    with ThreadPoolExecutor(max_workers=4) as executor:
        future_map = {executor.submit(fn): name for name, fn in tasks.items()}
        for future in as_completed(future_map):
            name = future_map[future]
            try:
                df = future.result()
                source_dfs[name] = df
                if not df.empty and progress_callback:
                    progress_callback(f"{name}: {len(df)} 个靶点")
            except Exception as exc:
                source_dfs[name] = pd.DataFrame()
                if progress_callback:
                    progress_callback(f"{name}: 查询失败 ({exc})")

    ot_df  = source_dfs.get("Open Targets", pd.DataFrame())
    up_df  = source_dfs.get("UniProt",       pd.DataFrame())
    ncbi_df = source_dfs.get("NCBI Gene",    pd.DataFrame())
    hz_df  = source_dfs.get("Harmonizome",   pd.DataFrame())

    if all(df.empty for df in [ot_df, up_df, ncbi_df, hz_df]):
        return pd.DataFrame()

    # Merge with priority: OT > UniProt > NCBI Gene > Harmonizome
    merged = ot_df.copy()
    seen_genes = set(merged["Gene"]) if not merged.empty else set()

    # Cap per-source contribution to avoid one source dominating
    SOURCE_CAP = {"Harmonizome": 400, "NCBI Gene": 200, "UniProt": 300}
    for extra_df in [up_df, ncbi_df, hz_df]:
        if extra_df.empty:
            continue
        new_rows = extra_df[~extra_df["Gene"].isin(seen_genes)]
        source = new_rows["Source"].iloc[0] if not new_rows.empty else ""
        cap = SOURCE_CAP.get(source, 500)
        new_rows = new_rows.head(cap)
        seen_genes.update(new_rows["Gene"].tolist())
        merged = pd.concat([merged, new_rows], ignore_index=True)

    if min_score > 0 and "Score" in merged.columns:
        merged = merged[merged["Score"] >= min_score]

    result = merged.reset_index(drop=True)
    cache_set(key, result, category="disease")
    return result
