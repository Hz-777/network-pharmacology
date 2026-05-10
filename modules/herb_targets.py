"""
HERB database target query module.
Retrieves compound→target associations from herb.ac.cn.

Workflow:
  compound name → HERB ingredient search → ingredient_id
               → HERB target list → gene symbols
"""

import requests
import pandas as pd
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

from modules.cache import cache_get, cache_set, make_key

HERB_API = "https://herb.ac.cn/api"
HEADERS = {
    "User-Agent": "Mozilla/5.0 (network-pharmacology-app/1.0)",
    "Referer": "https://herb.ac.cn/",
}

_RETRY = dict(
    retry=retry_if_exception_type((requests.Timeout, requests.ConnectionError)),
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=1, max=8),
    reraise=False,
)


@retry(**_RETRY)
def _search_ingredient(name: str) -> "str | None":
    """Search HERB for an ingredient by name, return ingredient_id."""
    r = requests.get(
        f"{HERB_API}/search/",
        params={"database": "Ingredient", "term": name, "limit": 3},
        headers=HEADERS, timeout=12,
    )
    items = r.json().get("data", {}).get("items", r.json().get("items", []))
    if items:
        item = items[0]
        return item.get("ingredient_id") or item.get("id") or item.get("HBIN")
    return None


@retry(**_RETRY)
def _get_ingredient_targets(ingredient_id: str) -> list:
    """Fetch targets for a HERB ingredient_id, return list of gene symbols."""
    # Try multiple possible endpoint patterns
    for endpoint in [
        f"{HERB_API}/ingredient/{ingredient_id}/target/",
        f"{HERB_API}/ingredient/{ingredient_id}/targets/",
        f"{HERB_API}/target/?ingredient_id={ingredient_id}&limit=500",
    ]:
        try:
            r = requests.get(endpoint, headers=HEADERS, timeout=15)
            if r.status_code != 200:
                continue
            data = r.json()
            items = data.get("data", {}).get("items", data.get("items", []))
            if not items:
                continue
            genes = []
            for item in items:
                gene = (item.get("gene_symbol") or item.get("Gene_Symbol")
                        or item.get("target_name") or item.get("Target"))
                if gene and isinstance(gene, str) and gene.strip():
                    genes.append(gene.strip())
            if genes:
                return genes
        except Exception:
            continue
    return []


def get_herb_targets(compound_name: str) -> pd.DataFrame:
    """
    Return DataFrame[Gene, Source] with HERB-sourced targets for compound_name.
    Returns empty DataFrame if HERB doesn't have data for this compound.
    """
    key = make_key("herb_targets", compound_name)
    cached = cache_get(key)
    if cached is not None:
        return cached

    result = pd.DataFrame()
    try:
        ingredient_id = _search_ingredient(compound_name)
        if not ingredient_id:
            cache_set(key, result)
            return result

        genes = _get_ingredient_targets(ingredient_id)
        if genes:
            result = pd.DataFrame({
                "Gene": sorted(set(genes)),
                "Source": "HERB",
            })
    except Exception:
        pass

    cache_set(key, result)
    return result


def get_herb_targets_batch(compound_names: list, max_workers: int = 3) -> pd.DataFrame:
    """
    Fetch HERB targets for multiple compounds in parallel.
    Returns merged DataFrame[Compound, Gene, Source].
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    rows = []

    def _fetch(name):
        df = get_herb_targets(name)
        if not df.empty:
            df = df.copy()
            df["Compound"] = name
        return df

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_fetch, nm): nm for nm in compound_names}
        for future in as_completed(futures):
            try:
                df = future.result()
                if not df.empty:
                    rows.append(df)
            except Exception:
                pass

    if not rows:
        return pd.DataFrame(columns=["Compound", "Gene", "Source"])
    return pd.concat(rows, ignore_index=True)
