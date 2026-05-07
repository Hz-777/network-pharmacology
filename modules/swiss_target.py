"""SwissTargetPrediction API wrapper for drug target prediction."""

import requests
import pandas as pd
import time
import json
from typing import Optional

SWISS_API = "https://www.swisstargetprediction.ch"


def predict_targets(smiles: str, organism: str = "Homo sapiens") -> pd.DataFrame:
    """
    Submit SMILES to SwissTargetPrediction and retrieve predicted targets.
    Returns DataFrame with columns: Target, Probability, UniprotID, etc.
    """
    # POST to prediction endpoint
    url = f"{SWISS_API}/predict.php"
    data = {
        "smiles": smiles,
        "organism": organism,
    }
    headers = {
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36",
        "Content-Type": "application/x-www-form-urlencoded",
        "Referer": f"{SWISS_API}/",
    }

    try:
        session = requests.Session()
        resp = session.post(url, data=data, headers=headers, timeout=60)
        if resp.status_code == 200:
            result_url = _extract_result_url(resp.text, smiles)
            if result_url:
                time.sleep(3)
                return _fetch_result_table(session, result_url, headers)
    except Exception:
        pass

    return pd.DataFrame()


def _extract_result_url(html: str, smiles: str) -> Optional[str]:
    """Extract result download URL from SwissTargetPrediction response."""
    from bs4 import BeautifulSoup
    soup = BeautifulSoup(html, "lxml")

    # Look for download link or result table link
    for a in soup.find_all("a", href=True):
        href = a["href"]
        if "target" in href.lower() or "result" in href.lower() or ".csv" in href.lower():
            if href.startswith("http"):
                return href
            return f"{SWISS_API}/{href.lstrip('/')}"

    # Try to find result ID in the page
    for script in soup.find_all("script"):
        if script.string and "resultID" in script.string:
            import re
            match = re.search(r'resultID\s*=\s*["\']?(\w+)["\']?', script.string)
            if match:
                return f"{SWISS_API}/result.php?resultID={match.group(1)}"

    return None


def _fetch_result_table(session: requests.Session, url: str, headers: dict) -> pd.DataFrame:
    """Fetch and parse the SwissTargetPrediction result table."""
    from bs4 import BeautifulSoup
    import re

    # Try CSV download first
    csv_url = url.replace("result.php", "download.php") if "result.php" in url else url
    if ".csv" not in csv_url:
        csv_url = csv_url + ("&" if "?" in csv_url else "?") + "format=csv"

    try:
        resp = session.get(csv_url, headers=headers, timeout=30)
        if "text/csv" in resp.headers.get("content-type", "") or resp.text.startswith("Target"):
            from io import StringIO
            return pd.read_csv(StringIO(resp.text))
    except Exception:
        pass

    # Fallback: parse HTML table
    resp = session.get(url, headers=headers, timeout=30)
    soup = BeautifulSoup(resp.content, "lxml")
    tables = soup.find_all("table")
    for table in tables:
        headers_row = [th.get_text(strip=True) for th in table.find_all("th")]
        if any(col in headers_row for col in ["Target", "Probability", "Gene name"]):
            rows = []
            for tr in table.find_all("tr")[1:]:
                tds = [td.get_text(strip=True) for td in tr.find_all("td")]
                if tds:
                    rows.append(tds)
            if rows:
                return pd.DataFrame(rows, columns=headers_row[:len(rows[0])])

    return pd.DataFrame()


def predict_targets_batch(
    compounds: pd.DataFrame,
    smiles_col: str = "SMILES",
    name_col: str = "mol_name",
    min_probability: float = 0.0,
    progress_callback=None,
) -> pd.DataFrame:
    """
    Predict targets for multiple compounds.
    Returns merged DataFrame with compound-target associations.
    """
    all_targets = []
    valid = compounds[compounds[smiles_col].notna() & (compounds[smiles_col] != "")]
    total = len(valid)

    for i, (_, row) in enumerate(valid.iterrows()):
        smiles = row[smiles_col]
        name = row.get(name_col, f"compound_{i}")

        if progress_callback:
            progress_callback(i + 1, total, name)

        targets_df = predict_targets(smiles)

        if not targets_df.empty:
            # Normalize probability column
            prob_col = next((c for c in targets_df.columns if "prob" in c.lower()), None)
            if prob_col:
                targets_df[prob_col] = pd.to_numeric(targets_df[prob_col], errors="coerce")
                targets_df = targets_df[targets_df[prob_col] > min_probability]

            targets_df["Compound"] = name
            all_targets.append(targets_df)

        time.sleep(2)  # Rate limiting

    if all_targets:
        return pd.concat(all_targets, ignore_index=True)
    return pd.DataFrame()


def get_unique_targets(targets_df: pd.DataFrame) -> list:
    """Extract unique gene symbols from targets DataFrame."""
    gene_cols = [c for c in targets_df.columns if any(k in c.lower() for k in ["gene", "symbol", "target"])]
    if not gene_cols:
        return []

    col = gene_cols[0]
    return sorted(targets_df[col].dropna().unique().tolist())
