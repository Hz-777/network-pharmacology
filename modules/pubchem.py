"""PubChem API wrapper for compound SMILES retrieval."""

import requests
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

from modules.cache import cache_get, cache_set, make_key

PUBCHEM_API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


@retry(
    retry=retry_if_exception_type((requests.Timeout, requests.ConnectionError)),
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=1, max=8),
    reraise=False,
)
def _fetch_pubchem(url: str) -> dict:
    resp = requests.get(url, timeout=15)
    if resp.status_code == 404:
        return {}
    resp.raise_for_status()
    return resp.json()


def get_compound_info(compound_name: str) -> dict:
    """Get compound info including SMILES, CID, formula and MW from PubChem."""
    result = {"name": compound_name, "CID": None, "SMILES": None, "formula": None, "MW": None}

    key = make_key("pubchem", compound_name)
    cached = cache_get(key)
    # Only use cache when SMILES was successfully retrieved; skip stale empty results
    if cached is not None and cached.get("SMILES"):
        return cached

    # PubChem returns CanonicalSMILES as "SMILES" and IsomericSMILES as "IsomericSMILES"
    url = (
        f"{PUBCHEM_API}/compound/name/{requests.utils.quote(compound_name)}"
        "/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES/JSON"
    )
    try:
        data = _fetch_pubchem(url)
        props = data.get("PropertyTable", {}).get("Properties", [])
        if props:
            p = props[0]
            result["CID"] = p.get("CID")
            # Field name varies: "IsomericSMILES", "SMILES", or "CanonicalSMILES"
            result["SMILES"] = (
                p.get("IsomericSMILES")
                or p.get("SMILES")
                or p.get("CanonicalSMILES")
                or p.get("ConnectivitySMILES")
            )
            result["formula"] = p.get("MolecularFormula")
            result["MW"] = p.get("MolecularWeight")
    except Exception:
        pass

    if result.get("SMILES"):  # only cache successful results
        cache_set(key, result, category="chemical")
    return result
