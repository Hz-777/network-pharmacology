"""
Disease target retrieval via Open Targets Platform (GraphQL API).
Free, no authentication required.
Supports both English and Chinese disease names.
"""

import re
import requests
import pandas as pd

OT_API = "https://api.platform.opentargets.org/api/v4/graphql"
HEADERS = {"Content-Type": "application/json"}

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
            headers=HEADERS,
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
            headers=HEADERS,
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


def get_disease_targets(disease: str, min_score: float = 0.0,
                        progress_callback=None) -> pd.DataFrame:
    """
    Retrieve disease-associated gene targets from Open Targets Platform.
    Accepts English or Chinese disease names.
    """
    # Translate Chinese to English if needed
    en_disease = _to_english(disease)
    if en_disease != disease and progress_callback:
        progress_callback(f"中文疾病名转换: {disease} → {en_disease}")

    if progress_callback:
        progress_callback(f"Open Targets: 搜索 [{en_disease}]...")

    efo_id = _search_disease_id(en_disease)
    if not efo_id:
        if progress_callback:
            progress_callback("未找到匹配疾病，使用内置数据")
        return pd.DataFrame()

    if progress_callback:
        progress_callback(f"Open Targets: 获取靶点（{efo_id}）...")

    df = _get_associated_targets(efo_id, size=150)
    if df.empty:
        return df

    if min_score > 0 and "Score" in df.columns:
        df = df[df["Score"] >= min_score]

    return df.reset_index(drop=True)
