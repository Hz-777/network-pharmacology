"""
Drug target prediction via ChEMBL database.
Replaces the defunct SwissTargetPrediction REST API.

Improvements over v1:
- No pChEMBL filter: includes all activity types (IC50, Ki, Kd, Inhibition…)
- No target_type filter: includes PROTEIN COMPLEX in addition to SINGLE PROTEIN
- Paginated fetch: retrieves up to 500 activity records (vs 200)
- Tries multiple ChEMBL IDs (salt forms, related compounds)
"""

import requests
import pandas as pd
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

from modules.cache import cache_get, cache_set, make_key

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
HEADERS = {"User-Agent": "Mozilla/5.0 (network-pharmacology-app/1.0)"}

_RETRY = dict(
    retry=retry_if_exception_type((requests.Timeout, requests.ConnectionError)),
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=1, max=10),
    reraise=False,
)


@retry(**_RETRY)
def _search_molecules(query: str, by: str = "name") -> list:
    """Search ChEMBL molecules by name or SMILES, return list of ChEMBL IDs."""
    if by == "name":
        r = requests.get(
            f"{CHEMBL_API}/molecule/search.json",
            params={"q": query, "limit": 5},
            headers=HEADERS, timeout=15,
        )
        return [m["molecule_chembl_id"] for m in r.json().get("molecules", [])]
    else:  # SMILES exact match
        r = requests.get(
            f"{CHEMBL_API}/molecule.json",
            params={"molecule_structures__canonical_smiles": query, "limit": 3},
            headers=HEADERS, timeout=15,
        )
        return [m["molecule_chembl_id"] for m in r.json().get("molecules", [])]


@retry(**_RETRY)
def _fetch_activities(chembl_id: str, offset: int = 0) -> list:
    """Fetch human bioactivity records for a ChEMBL ID (paginated)."""
    r = requests.get(
        f"{CHEMBL_API}/activity.json",
        params={
            "molecule_chembl_id": chembl_id,
            "target_organism": "Homo sapiens",
            "limit": 500,
            "offset": offset,
        },
        headers=HEADERS, timeout=25,
    )
    return r.json().get("activities", [])


@retry(**_RETRY)
def _fetch_targets(target_ids: list) -> list:
    """Fetch target details and extract gene symbols."""
    if not target_ids:
        return []
    r = requests.get(
        f"{CHEMBL_API}/target.json",
        params={
            "target_chembl_id__in": ",".join(target_ids[:100]),
            "organism": "Homo sapiens",
            "target_type__in": "SINGLE PROTEIN,PROTEIN COMPLEX",
            "limit": 100,
        },
        headers=HEADERS, timeout=20,
    )
    genes = []
    for t in r.json().get("targets", []):
        for comp in t.get("target_components", []):
            for syn in comp.get("target_component_synonyms", []):
                if syn.get("syn_type") == "GENE_SYMBOL":
                    genes.append(syn["component_synonym"])
                    break
    return genes


def predict_targets(smiles: str, compound_name: str = "") -> pd.DataFrame:
    """
    Return a DataFrame[Gene, ChEMBL_ID, Source] for the given compound.
    Queries ChEMBL with relaxed filters to maximise target coverage.
    """
    key = make_key("chembl_v2", smiles, compound_name)
    cached = cache_get(key)
    if cached is not None:
        return cached

    # Collect candidate ChEMBL IDs from name + SMILES
    chembl_ids: list[str] = []
    try:
        if compound_name:
            chembl_ids += _search_molecules(compound_name, by="name")
        if smiles and len(chembl_ids) < 3:
            chembl_ids += _search_molecules(smiles, by="smiles")
    except Exception:
        pass

    chembl_ids = list(dict.fromkeys(chembl_ids))  # deduplicate, preserve order

    if not chembl_ids:
        result = pd.DataFrame()
        cache_set(key, result)
        return result

    # Fetch activities for each ChEMBL ID and collect target IDs
    all_target_ids: set[str] = set()
    best_id = chembl_ids[0]
    for cid in chembl_ids[:3]:
        try:
            acts = _fetch_activities(cid)
            ids = {a["target_chembl_id"] for a in acts if a.get("target_chembl_id")}
            all_target_ids.update(ids)
        except Exception:
            pass

    if not all_target_ids:
        result = pd.DataFrame()
        cache_set(key, result)
        return result

    # Resolve target IDs → gene symbols
    try:
        genes = _fetch_targets(list(all_target_ids))
    except Exception:
        genes = []

    if not genes:
        result = pd.DataFrame()
        cache_set(key, result)
        return result

    result = pd.DataFrame({
        "Gene": sorted(set(genes)),
        "ChEMBL_ID": best_id,
        "Source": "ChEMBL",
    })
    cache_set(key, result)
    return result
