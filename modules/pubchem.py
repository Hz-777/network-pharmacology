"""PubChem API wrapper for compound SMILES retrieval."""

import requests
import time
import pandas as pd
from typing import Optional

PUBCHEM_API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def get_smiles_by_name(compound_name: str) -> Optional[str]:
    """Fetch canonical SMILES from PubChem by compound name."""
    url = f"{PUBCHEM_API}/compound/name/{requests.utils.quote(compound_name)}/property/CanonicalSMILES,IsomericSMILES/JSON"
    try:
        resp = requests.get(url, timeout=15)
        if resp.status_code == 200:
            data = resp.json()
            props = data["PropertyTable"]["Properties"]
            if props:
                return props[0].get("IsomericSMILES") or props[0].get("CanonicalSMILES")
    except Exception:
        pass
    return None


def get_cid_by_name(compound_name: str) -> Optional[int]:
    """Get PubChem CID for a compound name."""
    url = f"{PUBCHEM_API}/compound/name/{requests.utils.quote(compound_name)}/cids/JSON"
    try:
        resp = requests.get(url, timeout=15)
        if resp.status_code == 200:
            data = resp.json()
            cids = data.get("IdentifierList", {}).get("CID", [])
            return cids[0] if cids else None
    except Exception:
        pass
    return None


def get_compound_info(compound_name: str) -> dict:
    """Get full compound info including SMILES and CID."""
    url = (
        f"{PUBCHEM_API}/compound/name/{requests.utils.quote(compound_name)}"
        "/property/MolecularFormula,MolecularWeight,CanonicalSMILES,IsomericSMILES,IUPACName/JSON"
    )
    result = {"name": compound_name, "CID": None, "SMILES": None, "formula": None, "MW": None}
    try:
        resp = requests.get(url, timeout=15)
        if resp.status_code == 200:
            data = resp.json()
            props = data["PropertyTable"]["Properties"]
            if props:
                p = props[0]
                result["CID"] = p.get("CID")
                result["SMILES"] = p.get("IsomericSMILES") or p.get("CanonicalSMILES")
                result["formula"] = p.get("MolecularFormula")
                result["MW"] = p.get("MolecularWeight")
    except Exception:
        pass
    return result


def batch_get_smiles(compound_names: list, progress_callback=None, delay: float = 0.3) -> pd.DataFrame:
    """Retrieve SMILES for multiple compounds with rate limiting."""
    results = []
    total = len(compound_names)
    for i, name in enumerate(compound_names):
        info = get_compound_info(name)
        results.append(info)
        if progress_callback:
            progress_callback(i + 1, total, name)
        time.sleep(delay)
    return pd.DataFrame(results)
