"""TCMSP database scraper for herb compound retrieval."""

import requests
import pandas as pd
import time
import json
from bs4 import BeautifulSoup


TCMSP_BASE = "https://old.tcmsp-e.com"
HERB_BASE = "https://herb.ac.cn"

HEADERS = {
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36",
    "Accept": "application/json, text/javascript, */*",
}


def search_herb_tcmsp(herb_name: str, ob_threshold: float = 30.0, dl_threshold: float = 0.18) -> pd.DataFrame:
    """
    Query TCMSP for herb compounds. Returns filtered DataFrame.
    Uses the TCMSP JSON API endpoint discovered from network inspection.
    """
    session = requests.Session()
    session.headers.update(HEADERS)

    # Step 1: search for herb to get its ID
    search_url = f"{TCMSP_BASE}/tcmsp.php"
    params = {
        "qr": herb_name,
        "qsr": "herb_cn_name",
        "token": "",
    }

    try:
        resp = session.get(search_url, params=params, timeout=30)
        resp.raise_for_status()
    except Exception as e:
        raise ConnectionError(f"TCMSP连接失败: {e}")

    # Try JSON API approach
    api_url = f"{TCMSP_BASE}/api/getmolecules.php"
    payload = {
        "herb_cn_name": herb_name,
        "token": "",
    }

    try:
        resp = session.post(api_url, data=payload, timeout=30)
        data = resp.json()
        if isinstance(data, list) and len(data) > 0:
            df = pd.DataFrame(data)
            return _filter_adme(df, ob_threshold, dl_threshold)
    except Exception:
        pass

    # Fallback: scrape HTML
    return _scrape_tcmsp_html(session, herb_name, ob_threshold, dl_threshold)


def _scrape_tcmsp_html(session, herb_name: str, ob_threshold: float, dl_threshold: float) -> pd.DataFrame:
    """Fallback HTML scraping for TCMSP."""
    url = f"{TCMSP_BASE}/tcmsp.php"
    params = {"qr": herb_name, "qsr": "herb_cn_name"}

    resp = session.get(url, params=params, timeout=30)
    soup = BeautifulSoup(resp.content, "lxml")

    # Find molecule table
    tables = soup.find_all("table")
    for table in tables:
        headers = [th.get_text(strip=True) for th in table.find_all("th")]
        if "Molecule Name" in headers or "mol_name" in headers:
            rows = []
            for tr in table.find_all("tr")[1:]:
                tds = [td.get_text(strip=True) for td in tr.find_all("td")]
                if tds:
                    rows.append(tds)
            if rows:
                df = pd.DataFrame(rows, columns=headers[:len(rows[0])])
                return _filter_adme(df, ob_threshold, dl_threshold)

    return pd.DataFrame()


def _filter_adme(df: pd.DataFrame, ob_threshold: float, dl_threshold: float) -> pd.DataFrame:
    """Filter compounds by OB and DL thresholds."""
    # Normalize column names
    col_map = {}
    for col in df.columns:
        lower = col.lower()
        if "ob" in lower or "oral bioavailability" in lower:
            col_map[col] = "OB"
        elif "dl" in lower or "drug-likeness" in lower:
            col_map[col] = "DL"
        elif "mol_name" in lower or "molecule name" in lower:
            col_map[col] = "mol_name"
        elif "smiles" in lower:
            col_map[col] = "SMILES"
        elif "mw" in lower or "mol_weight" in lower:
            col_map[col] = "MW"

    df = df.rename(columns=col_map)

    for num_col in ["OB", "DL", "MW"]:
        if num_col in df.columns:
            df[num_col] = pd.to_numeric(df[num_col], errors="coerce")

    filtered = df.copy()
    if "OB" in df.columns:
        filtered = filtered[filtered["OB"] >= ob_threshold]
    if "DL" in df.columns:
        filtered = filtered[filtered["DL"] >= dl_threshold]

    return filtered.reset_index(drop=True)


def get_herb_targets_from_tcmsp(herb_name: str) -> pd.DataFrame:
    """Get herb-target associations directly from TCMSP targets page."""
    session = requests.Session()
    session.headers.update(HEADERS)

    api_url = f"{TCMSP_BASE}/api/gettargets.php"
    try:
        resp = session.post(api_url, data={"herb_cn_name": herb_name}, timeout=30)
        data = resp.json()
        if data:
            return pd.DataFrame(data)
    except Exception:
        pass
    return pd.DataFrame()
