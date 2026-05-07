"""
Herb compound retrieval.
Priority: HERB API (herb.ac.cn) → TCMSP scraping → built-in fallback
"""

import requests
import pandas as pd
import json
import time
from pathlib import Path

HEADERS = {
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36",
    "Accept": "application/json, text/html, */*",
}

BUILTIN_PATH = Path(__file__).parent.parent / "data" / "tcmsp_builtin.json"


# ── 1. HERB database API (herb.ac.cn) ────────────────────────────────────────

def _query_herb_api(herb_name: str) -> pd.DataFrame:
    """Query HERB database API for herb compounds."""
    session = requests.Session()
    session.headers.update(HEADERS)

    # Search herb
    search_url = "https://herb.ac.cn/api/herb/search"
    try:
        resp = session.get(search_url, params={"keyword": herb_name, "page": 1, "limit": 5}, timeout=15)
        data = resp.json()
        items = data.get("data", {}).get("items", data.get("items", []))
        if not items:
            return pd.DataFrame()

        herb_id = items[0].get("herb_id") or items[0].get("id")
        if not herb_id:
            return pd.DataFrame()
    except Exception:
        return pd.DataFrame()

    # Get compounds for herb
    comp_url = f"https://herb.ac.cn/api/herb/{herb_id}/molecules"
    try:
        resp = session.get(comp_url, params={"page": 1, "limit": 200}, timeout=15)
        data = resp.json()
        mols = data.get("data", {}).get("items", data.get("items", []))
        if not mols:
            return pd.DataFrame()

        rows = []
        for m in mols:
            rows.append({
                "mol_name": m.get("mol_name") or m.get("name") or m.get("molecule_name", ""),
                "OB":       float(m.get("ob") or m.get("OB") or 0),
                "DL":       float(m.get("dl") or m.get("DL") or 0),
                "MW":       float(m.get("mw") or m.get("MW") or 0),
                "SMILES":   m.get("smiles") or m.get("SMILES") or "",
            })
        return pd.DataFrame(rows)
    except Exception:
        return pd.DataFrame()


# ── 2. TCMSP scraping (fallback) ─────────────────────────────────────────────

def _query_tcmsp(herb_name: str) -> pd.DataFrame:
    """Scrape TCMSP for herb compounds."""
    session = requests.Session()
    session.headers.update(HEADERS)

    for base in ["https://www.tcmsp-e.com", "https://old.tcmsp-e.com"]:
        try:
            resp = session.get(
                f"{base}/tcmsp.php",
                params={"qr": herb_name, "qsr": "herb_cn_name", "token": ""},
                timeout=12,
            )
            if resp.status_code != 200:
                continue

            # Try JSON API endpoint
            api_resp = session.post(
                f"{base}/api/getmolecules.php",
                data={"herb_cn_name": herb_name, "token": ""},
                timeout=12,
            )
            data = api_resp.json()
            if isinstance(data, list) and data:
                df = pd.DataFrame(data)
                return _normalize_cols(df)
        except Exception:
            continue

    return pd.DataFrame()


def _normalize_cols(df: pd.DataFrame) -> pd.DataFrame:
    col_map = {}
    for c in df.columns:
        lo = c.lower()
        if "mol_name" in lo or "molecule name" in lo:
            col_map[c] = "mol_name"
        elif lo in ("ob", "oral bioavailability"):
            col_map[c] = "OB"
        elif lo in ("dl", "drug-likeness"):
            col_map[c] = "DL"
        elif "smiles" in lo:
            col_map[c] = "SMILES"
        elif lo in ("mw", "mol_weight", "molecular weight"):
            col_map[c] = "MW"
    return df.rename(columns=col_map)


# ── 3. Built-in data ──────────────────────────────────────────────────────────

def _query_builtin(herb_name: str) -> pd.DataFrame:
    """Load pre-bundled TCMSP data for common herbs."""
    try:
        with open(BUILTIN_PATH, encoding="utf-8") as f:
            db = json.load(f)
        rows = db.get(herb_name)
        if rows:
            return pd.DataFrame(rows)
    except Exception:
        pass
    return pd.DataFrame()


# ── Public API ────────────────────────────────────────────────────────────────

def search_herb_tcmsp(herb_name: str, ob_threshold: float = 30.0,
                      dl_threshold: float = 0.18) -> pd.DataFrame:
    """
    Retrieve herb compounds with ADME filtering.
    Tries HERB API → TCMSP scraping → built-in data, in order.
    Returns empty DataFrame only if all three fail.
    """
    df = pd.DataFrame()

    # 1. HERB API
    try:
        df = _query_herb_api(herb_name)
        if not df.empty:
            df["_source"] = "HERB"
    except Exception:
        pass

    # 2. TCMSP scraping
    if df.empty:
        try:
            df = _query_tcmsp(herb_name)
            if not df.empty:
                df["_source"] = "TCMSP"
        except Exception:
            pass

    # 3. Built-in
    if df.empty:
        df = _query_builtin(herb_name)
        if not df.empty:
            df["_source"] = "内置数据"

    if df.empty:
        return df

    # ADME filter
    for col in ["OB", "DL"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)

    if "OB" in df.columns:
        df = df[df["OB"] >= ob_threshold]
    if "DL" in df.columns:
        df = df[df["DL"] >= dl_threshold]

    return df.reset_index(drop=True)
