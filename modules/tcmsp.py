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
    """
    Query TCMSP via its Kendo-Grid JSON endpoint.

    Flow:
      1. GET tcmsp.php  → extract session token
      2. GET tcmspsearch.php?qs=herb_all_name&q=<name>  → list of herbs + en_name
      3. GET tcmspsearch.php?qr=<en_name>&qsr=herb_en_name  → compound data (JSON in page)
    """
    import re
    import json as _json

    session = requests.Session()
    session.headers.update(HEADERS)

    try:
        # Step 1 – get session token
        r = session.get("https://www.tcmsp-e.com/tcmsp.php", timeout=12)
        tokens = re.findall(r"value=['\"]([a-f0-9]{32})['\"]", r.text)
        token = tokens[0] if tokens else ""

        # Step 2 – search for herb by name to get English name
        r2 = session.get(
            "https://www.tcmsp-e.com/tcmspsearch.php",
            params={"qs": "herb_all_name", "q": herb_name, "token": token},
            timeout=15,
        )
        m = re.search(r"data:\s*(\[(?:[^\[\]]*|\[[^\[\]]*\])*\])", r2.text, re.DOTALL)
        if not m or len(m.group(1).strip()) <= 2:
            return pd.DataFrame()

        herbs_list = _json.loads(m.group(1))
        # Prefer exact Chinese name match
        herb_info = next((h for h in herbs_list if h.get("herb_cn_name") == herb_name),
                         herbs_list[0] if herbs_list else None)
        if not herb_info:
            return pd.DataFrame()
        en_name = herb_info.get("herb_en_name", "")
        if not en_name:
            return pd.DataFrame()

        # Step 3 – get compound detail page
        r3 = session.get(
            "https://www.tcmsp-e.com/tcmspsearch.php",
            params={"qr": en_name, "qsr": "herb_en_name", "token": token},
            timeout=30,
        )
        start = r3.text.find('data: [{"')
        if start == -1:
            return pd.DataFrame()
        bracket = r3.text.index("[", start)
        decoder = _json.JSONDecoder()
        raw, _ = decoder.raw_decode(r3.text[bracket:])
        if not isinstance(raw, list) or not raw:
            return pd.DataFrame()

        rows = []
        for c in raw:
            rows.append({
                "mol_name": (c.get("molecule_name") or "").strip(),
                "OB":  float(c.get("ob") or 0),
                "DL":  float(c.get("dl") or 0),
                "MW":  float(c.get("mw") or 0),
                "SMILES": (c.get("smiles") or c.get("SMILES") or "").strip(),
            })
        return pd.DataFrame(rows)

    except Exception:
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
