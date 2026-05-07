"""
Disease target retrieval via Open Targets Platform (GraphQL API).
Free, no authentication required.
"""

import requests
import pandas as pd

OT_API = "https://api.platform.opentargets.org/api/v4/graphql"
HEADERS = {"Content-Type": "application/json"}


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

    Args:
        disease:  Disease name (English or Chinese keywords)
        min_score: Minimum association score [0–1]
        progress_callback: Optional callable(str) for status messages
    """
    if progress_callback:
        progress_callback(f"Open Targets: 搜索疾病 [{disease}]...")

    efo_id = _search_disease_id(disease)
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
