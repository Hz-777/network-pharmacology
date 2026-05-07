"""GO/KEGG enrichment analysis via Enrichr API."""

import requests
import pandas as pd
import time
from typing import Optional


ENRICHR_API = "https://maayanlab.cloud/Enrichr"

ENRICHR_LIBRARIES = {
    "GO_BP": "GO_Biological_Process_2023",
    "GO_CC": "GO_Cellular_Component_2023",
    "GO_MF": "GO_Molecular_Function_2023",
    "KEGG": "KEGG_2021_Human",
}


def enrichr_submit_genes(genes: list) -> Optional[str]:
    """Submit gene list to Enrichr, return user_list_id."""
    genes_str = "\n".join(genes)
    payload = {
        "list": (None, genes_str),
        "description": (None, "NetworkPharmacology"),
    }
    try:
        resp = requests.post(f"{ENRICHR_API}/addList", files=payload, timeout=30)
        resp.raise_for_status()
        return resp.json().get("userListId")
    except Exception as e:
        print(f"Enrichr submit error: {e}")
        return None


def enrichr_get_enrichment(user_list_id: str, library: str) -> pd.DataFrame:
    """Get enrichment results from Enrichr for a given library."""
    try:
        resp = requests.get(
            f"{ENRICHR_API}/enrich",
            params={"userListId": user_list_id, "backgroundType": library},
            timeout=30,
        )
        resp.raise_for_status()
        results = resp.json().get(library, [])
        rows = []
        for item in results:
            if len(item) >= 7:
                rows.append({
                    "Rank":           item[0],
                    "Term":           item[1],
                    "P_value":        item[2],
                    "Z_score":        item[3],
                    "Combined_score": item[4],
                    "Genes":          ";".join(item[5]) if isinstance(item[5], list) else item[5],
                    "Adj_P_value":    item[6],
                    "Library":        library,
                })
        df = pd.DataFrame(rows)
        if not df.empty:
            df = df[df["P_value"] <= 0.05].sort_values("Combined_score", ascending=False)
        return df
    except Exception as e:
        print(f"Enrichr enrichment error ({library}): {e}")
        return pd.DataFrame()


def run_enrichr_analysis(genes: list, progress_callback=None) -> dict:
    """
    Run full GO+KEGG enrichment via Enrichr.
    Returns dict with keys: GO_BP, GO_CC, GO_MF, KEGG.
    """
    results = {}
    if progress_callback:
        progress_callback("提交基因列表到 Enrichr...")
    user_list_id = enrichr_submit_genes(genes)
    if not user_list_id:
        return results
    time.sleep(1)
    for key, library in ENRICHR_LIBRARIES.items():
        if progress_callback:
            progress_callback(f"分析 {key} 通路...")
        results[key] = enrichr_get_enrichment(user_list_id, library)
        time.sleep(0.5)
    return results
