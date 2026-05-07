"""GO/KEGG enrichment analysis via Enrichr and DAVID APIs."""

import requests
import pandas as pd
import time
import json
from typing import Optional


ENRICHR_API = "https://maayanlab.cloud/Enrichr"
DAVID_API = "https://david.ncifcrf.gov/api.jsp"


# --- Enrichr (primary, no auth needed) ---

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
        data = resp.json()
        return data.get("userListId")
    except Exception as e:
        print(f"Enrichr submit error: {e}")
        return None


def enrichr_get_enrichment(user_list_id: str, library: str) -> pd.DataFrame:
    """Get enrichment results from Enrichr for a given library."""
    url = f"{ENRICHR_API}/enrich"
    params = {
        "userListId": user_list_id,
        "backgroundType": library,
    }
    try:
        resp = requests.get(url, params=params, timeout=30)
        resp.raise_for_status()
        data = resp.json()
        results = data.get(library, [])

        rows = []
        for item in results:
            # Enrichr format: [rank, term, p_value, z_score, combined_score, genes, adj_p_value, ...]
            if len(item) >= 7:
                rows.append({
                    "Rank": item[0],
                    "Term": item[1],
                    "P_value": item[2],
                    "Z_score": item[3],
                    "Combined_score": item[4],
                    "Genes": ";".join(item[5]) if isinstance(item[5], list) else item[5],
                    "Adj_P_value": item[6],
                    "Library": library,
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
    Run full GO+KEGG enrichment analysis via Enrichr.
    Returns dict with keys: GO_BP, GO_CC, GO_MF, KEGG
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

        df = enrichr_get_enrichment(user_list_id, library)
        results[key] = df
        time.sleep(0.5)

    return results


# --- DAVID (fallback) ---

def david_add_list(genes: list, id_type: str = "OFFICIAL_GENE_SYMBOL") -> Optional[str]:
    """Submit gene list to DAVID API."""
    genes_str = ",".join(genes)
    url = f"{DAVID_API}"
    params = {
        "tool": "addList",
        "ids": genes_str,
        "idType": id_type,
        "listName": "NP_analysis",
        "listType": "Gene",
    }
    try:
        resp = requests.get(url, params=params, timeout=30)
        if resp.status_code == 200:
            return "david_list"
    except Exception as e:
        print(f"DAVID add list error: {e}")
    return None


def get_enrichment_summary(results: dict) -> pd.DataFrame:
    """Compile enrichment results into a summary table."""
    all_rows = []
    for category, df in results.items():
        if df is not None and not df.empty:
            df = df.copy()
            df["Category"] = category
            all_rows.append(df)

    if all_rows:
        return pd.concat(all_rows, ignore_index=True)
    return pd.DataFrame()


def get_top_kegg_pathways(kegg_df: pd.DataFrame, top_n: int = 20) -> pd.DataFrame:
    """Get top KEGG pathways sorted by combined score."""
    if kegg_df.empty:
        return pd.DataFrame()
    return kegg_df.head(top_n)


def get_top_go_terms(go_df: pd.DataFrame, top_n: int = 10) -> pd.DataFrame:
    """Get top GO terms sorted by combined score."""
    if go_df.empty:
        return pd.DataFrame()
    return go_df.head(top_n)
