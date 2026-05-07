"""PubChem API wrapper for compound SMILES retrieval."""

import requests

PUBCHEM_API = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def get_compound_info(compound_name: str) -> dict:
    """Get compound info including SMILES, CID, formula and MW from PubChem."""
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
