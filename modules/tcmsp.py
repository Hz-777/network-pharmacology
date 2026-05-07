"""
Herb compound retrieval.
- 本地运行: 优先 TCMSP → HERB API → 内置数据
- 云端运行: 优先 HERB API → TCMSP → 内置数据  (TCMSP 在云端 IP 受限)
"""

import os
import requests
import pandas as pd
import json
import time
from pathlib import Path
from bs4 import BeautifulSoup

HEADERS = {
    "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
                  "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36",
    "Accept": "application/json, text/html, */*",
    "Accept-Language": "zh-CN,zh;q=0.9,en;q=0.8",
}

BUILTIN_PATH = Path(__file__).parent.parent / "data" / "tcmsp_builtin.json"

# 检测是否在 Streamlit Cloud 运行
_ON_CLOUD = os.environ.get("STREAMLIT_SHARING_MODE") == "streamlit" \
            or os.path.exists("/mount/src")


# ── TCMSP ─────────────────────────────────────────────────────────────────────

def _query_tcmsp(herb_name: str) -> pd.DataFrame:
    """Query TCMSP database (works best on local network)."""
    session = requests.Session()
    session.headers.update(HEADERS)

    for base_url in [
        "https://www.tcmsp-e.com",
        "https://old.tcmsp-e.com",
    ]:
        try:
            # First visit the main page to get session cookies
            session.get(f"{base_url}/tcmsp.php", timeout=10)

            # Try JSON API
            resp = session.post(
                f"{base_url}/api/getmolecules.php",
                data={"herb_cn_name": herb_name, "token": ""},
                timeout=15,
            )
            if resp.status_code == 200:
                try:
                    data = resp.json()
                    if isinstance(data, list) and data:
                        return _normalize(pd.DataFrame(data))
                except Exception:
                    pass

            # Fallback: HTML scraping
            resp = session.get(
                f"{base_url}/tcmsp.php",
                params={"qr": herb_name, "qsr": "herb_cn_name"},
                timeout=15,
            )
            df = _parse_tcmsp_html(resp.text)
            if not df.empty:
                return df

        except Exception:
            continue

    return pd.DataFrame()


def _parse_tcmsp_html(html: str) -> pd.DataFrame:
    """Parse TCMSP HTML table into DataFrame."""
    soup = BeautifulSoup(html, "lxml")
    for table in soup.find_all("table"):
        headers = [th.get_text(strip=True) for th in table.find_all("th")]
        if any(k in " ".join(headers) for k in ["Molecule", "mol_name", "OB", "DL"]):
            rows = []
            for tr in table.find_all("tr")[1:]:
                tds = [td.get_text(strip=True) for td in tr.find_all("td")]
                if tds:
                    rows.append(dict(zip(headers, tds)))
            if rows:
                return _normalize(pd.DataFrame(rows))
    return pd.DataFrame()


# ── HERB API ──────────────────────────────────────────────────────────────────

def _query_herb_api(herb_name: str) -> pd.DataFrame:
    """Query HERB database (herb.ac.cn) — accessible from cloud servers."""
    session = requests.Session()
    session.headers.update(HEADERS)

    try:
        # Search for herb
        resp = session.get(
            "https://herb.ac.cn/api/herb/search",
            params={"keyword": herb_name, "page": 1, "limit": 5},
            timeout=15,
        )
        data = resp.json()
        items = (data.get("data") or {}).get("items") or data.get("items") or []
        if not items:
            return pd.DataFrame()

        herb_id = items[0].get("herb_id") or items[0].get("id")
        if not herb_id:
            return pd.DataFrame()

        # Get molecules
        resp = session.get(
            f"https://herb.ac.cn/api/herb/{herb_id}/molecules",
            params={"page": 1, "limit": 500},
            timeout=15,
        )
        data = resp.json()
        mols = (data.get("data") or {}).get("items") or data.get("items") or []
        if not mols:
            return pd.DataFrame()

        rows = []
        for m in mols:
            rows.append({
                "mol_name": m.get("mol_name") or m.get("name") or m.get("molecule_name", ""),
                "OB":    float(m.get("ob") or m.get("OB") or 0),
                "DL":    float(m.get("dl") or m.get("DL") or 0),
                "MW":    float(m.get("mw") or m.get("MW") or 0),
                "SMILES": m.get("smiles") or m.get("SMILES") or "",
            })
        return pd.DataFrame(rows)

    except Exception:
        return pd.DataFrame()


# ── Built-in fallback ─────────────────────────────────────────────────────────

def _query_builtin(herb_name: str) -> pd.DataFrame:
    try:
        with open(BUILTIN_PATH, encoding="utf-8") as f:
            db = json.load(f)
        rows = db.get(herb_name)
        if rows:
            df = pd.DataFrame(rows)
            df["_source"] = "内置数据"
            return df
    except Exception:
        pass
    return pd.DataFrame()


# ── Column normalizer ─────────────────────────────────────────────────────────

def _normalize(df: pd.DataFrame) -> pd.DataFrame:
    col_map = {}
    for c in df.columns:
        lo = c.lower().replace(" ", "_")
        if "mol_name" in lo or "molecule_name" in lo:
            col_map[c] = "mol_name"
        elif lo in ("ob",):
            col_map[c] = "OB"
        elif lo in ("dl",):
            col_map[c] = "DL"
        elif "smiles" in lo:
            col_map[c] = "SMILES"
        elif lo in ("mw", "mol_weight", "molecular_weight"):
            col_map[c] = "MW"
    df = df.rename(columns=col_map)
    for col in ["OB", "DL", "MW"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0)
    return df


# ── Public entry point ────────────────────────────────────────────────────────

def search_herb_tcmsp(herb_name: str, ob_threshold: float = 30.0,
                      dl_threshold: float = 0.18) -> pd.DataFrame:
    """
    Retrieve & filter herb compounds by ADME criteria.

    Source priority:
      Local:  TCMSP → HERB API → Built-in
      Cloud:  HERB API → TCMSP → Built-in
    """
    sources = (
        [("TCMSP",    _query_tcmsp),
         ("HERB API", _query_herb_api)]
        if not _ON_CLOUD else
        [("HERB API", _query_herb_api),
         ("TCMSP",    _query_tcmsp)]
    )

    df = pd.DataFrame()
    for name, fn in sources:
        try:
            df = fn(herb_name)
            if not df.empty:
                df["_source"] = name
                break
        except Exception:
            continue

    # Final fallback
    if df.empty:
        df = _query_builtin(herb_name)

    if df.empty:
        return df

    # ADME filter
    if "OB" in df.columns:
        df = df[df["OB"] >= ob_threshold]
    if "DL" in df.columns:
        df = df[df["DL"] >= dl_threshold]

    return df.reset_index(drop=True)
