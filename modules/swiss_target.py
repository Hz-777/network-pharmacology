"""
Drug target prediction via ChEMBL database.
Replaces the defunct SwissTargetPrediction REST API.

Workflow:
  compound name/SMILES → ChEMBL ID → bioactivity records → gene symbols
"""

import requests
import pandas as pd
import time

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
HEADERS = {"User-Agent": "Mozilla/5.0 (network-pharmacology-app/1.0)"}


def _chembl_id_from_name(name: str) -> "str | None":
    """Look up ChEMBL ID by compound preferred name."""
    try:
        r = requests.get(
            f"{CHEMBL_API}/molecule/search",
            params={"q": name, "format": "json", "limit": 1},
            headers=HEADERS,
            timeout=15,
        )
        mols = r.json().get("molecules", [])
        return mols[0]["molecule_chembl_id"] if mols else None
    except Exception:
        return None


def _chembl_id_from_smiles(smiles: str) -> "str | None":
    """Look up ChEMBL ID by SMILES (exact structure match)."""
    try:
        r = requests.get(
            f"{CHEMBL_API}/molecule",
            params={
                "molecule_structures__canonical_smiles": smiles,
                "format": "json",
                "limit": 1,
            },
            headers=HEADERS,
            timeout=15,
        )
        mols = r.json().get("molecules", [])
        return mols[0]["molecule_chembl_id"] if mols else None
    except Exception:
        return None


def _get_gene_symbols(chembl_id: str, max_targets: int = 50) -> list:
    """
    Return human gene symbols for a ChEMBL compound ID.
    Uses bioactivity records (requires pChEMBL ≥ 5, i.e., affinity ≤ 10 µM).
    """
    try:
        # Step 1: get activities with measured affinity
        r = requests.get(
            f"{CHEMBL_API}/activity",
            params={
                "molecule_chembl_id": chembl_id,
                "target_organism": "Homo sapiens",
                "pchembl_value__isnull": False,
                "format": "json",
                "limit": 200,
            },
            headers=HEADERS,
            timeout=20,
        )
        acts = r.json().get("activities", [])
        target_ids = list({a["target_chembl_id"] for a in acts if a.get("target_chembl_id")})
        if not target_ids:
            return []

        # Step 2: get gene symbols for those targets
        r2 = requests.get(
            f"{CHEMBL_API}/target",
            params={
                "target_chembl_id__in": ",".join(target_ids[:max_targets]),
                "organism": "Homo sapiens",
                "target_type": "SINGLE PROTEIN",
                "format": "json",
                "limit": max_targets,
            },
            headers=HEADERS,
            timeout=20,
        )
        genes = []
        for t in r2.json().get("targets", []):
            for comp in t.get("target_components", []):
                for syn in comp.get("target_component_synonyms", []):
                    if syn.get("syn_type") == "GENE_SYMBOL":
                        genes.append(syn["component_synonym"])
                        break
        return sorted(set(genes))

    except Exception:
        return []


def predict_targets(smiles: str, compound_name: str = "") -> pd.DataFrame:
    """
    Main entry point: given a SMILES (and optionally a name),
    return a DataFrame with columns [Gene, Source].
    """
    chembl_id = None

    # Try name first (faster), then SMILES
    if compound_name:
        chembl_id = _chembl_id_from_name(compound_name)
    if not chembl_id and smiles:
        chembl_id = _chembl_id_from_smiles(smiles)

    if not chembl_id:
        return pd.DataFrame()

    genes = _get_gene_symbols(chembl_id)
    if not genes:
        return pd.DataFrame()

    return pd.DataFrame({"Gene": genes, "ChEMBL_ID": chembl_id, "Source": "ChEMBL"})
