"""Disease target retrieval from GeneCards, DisGeNET, and OMIM."""

import requests
import pandas as pd
import time
from bs4 import BeautifulSoup
from typing import Optional


DISGENET_API = "https://www.disgenet.org/api"
GENECARDS_BASE = "https://www.genecards.org"
OMIM_API = "https://api.omim.org/api"


def get_disease_targets_genecards(disease: str, min_score: float = 0) -> pd.DataFrame:
    """
    Scrape GeneCards for disease-associated genes.
    Note: GeneCards may require cookies/session; uses best-effort scraping.
    """
    session = requests.Session()
    session.headers.update({
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 Chrome/120.0.0.0 Safari/537.36",
        "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
        "Accept-Language": "en-US,en;q=0.5",
    })

    search_url = f"{GENECARDS_BASE}/Search/Search.aspx"
    params = {
        "queryString": disease,
        "type": "GENE",
        "startPage": 1,
        "pageSize": 200,
    }

    genes = []
    try:
        resp = session.get(search_url, params=params, timeout=30)
        soup = BeautifulSoup(resp.content, "lxml")

        # Find gene cards / table
        for article in soup.find_all(["article", "div"], class_=lambda c: c and "gene" in c.lower()):
            symbol = article.find(class_=lambda c: c and "symbol" in str(c).lower())
            if symbol:
                text = symbol.get_text(strip=True)
                if text and text.isupper() and len(text) <= 20:
                    genes.append({"Gene": text, "Source": "GeneCards"})

        # Also check table rows
        for table in soup.find_all("table"):
            for tr in table.find_all("tr")[1:]:
                tds = tr.find_all("td")
                if tds:
                    gene_text = tds[0].get_text(strip=True)
                    if gene_text and len(gene_text) <= 20:
                        score = tds[-1].get_text(strip=True) if len(tds) > 1 else "0"
                        genes.append({"Gene": gene_text, "Score": score, "Source": "GeneCards"})
    except Exception as e:
        print(f"GeneCards scraping error: {e}")

    if genes:
        return pd.DataFrame(genes).drop_duplicates(subset=["Gene"])
    return pd.DataFrame()


def get_disease_targets_disgenet(disease: str, min_score: float = 0.0) -> pd.DataFrame:
    """
    Query DisGeNET API for disease-gene associations.
    Uses the public API (no authentication needed for basic queries).
    """
    # Search for disease ID first
    search_url = f"{DISGENET_API}/disease/search/"
    params = {"q": disease, "limit": 10}
    headers = {
        "accept": "application/json",
        "Authorization": "Bearer ",  # Public access
    }

    try:
        resp = requests.get(search_url, params=params, headers=headers, timeout=20)
        if resp.status_code == 200:
            diseases = resp.json().get("payload", {}).get("results", [])
            if not diseases:
                diseases = resp.json() if isinstance(resp.json(), list) else []

            if diseases:
                disease_id = diseases[0].get("diseaseId") or diseases[0].get("diseaseid")
                if disease_id:
                    return _fetch_disgenet_genes(disease_id, min_score)
    except Exception as e:
        print(f"DisGeNET search error: {e}")

    # Try direct gene-disease search
    return _search_disgenet_genes_direct(disease, min_score)


def _fetch_disgenet_genes(disease_id: str, min_score: float) -> pd.DataFrame:
    """Fetch gene associations for a DisGeNET disease ID."""
    url = f"{DISGENET_API}/gda/disease/{disease_id}"
    params = {"min_score": min_score, "limit": 500}
    headers = {"accept": "application/json"}

    try:
        resp = requests.get(url, params=params, headers=headers, timeout=30)
        if resp.status_code == 200:
            data = resp.json()
            results = data.get("payload", {}).get("results", data if isinstance(data, list) else [])
            genes = []
            for item in results:
                gene = item.get("gene_symbol") or item.get("geneName") or item.get("symbol")
                score = item.get("score", 0)
                if gene:
                    genes.append({"Gene": gene, "Score": score, "Source": "DisGeNET"})
            if genes:
                return pd.DataFrame(genes).drop_duplicates(subset=["Gene"])
    except Exception as e:
        print(f"DisGeNET fetch error: {e}")
    return pd.DataFrame()


def _search_disgenet_genes_direct(disease: str, min_score: float) -> pd.DataFrame:
    """Direct gene-disease search on DisGeNET."""
    url = f"{DISGENET_API}/gda/search/"
    params = {
        "q": disease,
        "min_score": min_score,
        "limit": 500,
        "source": "ALL",
    }
    headers = {"accept": "application/json"}

    try:
        resp = requests.get(url, params=params, headers=headers, timeout=30)
        if resp.status_code == 200:
            data = resp.json()
            results = data if isinstance(data, list) else data.get("payload", {}).get("results", [])
            genes = []
            for item in results:
                gene = item.get("gene_symbol") or item.get("geneName")
                score = item.get("score", 0)
                if gene:
                    genes.append({"Gene": gene, "Score": score, "Source": "DisGeNET"})
            if genes:
                return pd.DataFrame(genes).drop_duplicates(subset=["Gene"])
    except Exception as e:
        print(f"DisGeNET direct search error: {e}")
    return pd.DataFrame()


def get_disease_targets_omim(disease: str) -> pd.DataFrame:
    """Query OMIM for disease gene associations (requires API key)."""
    # OMIM requires registration; skip if no key
    return pd.DataFrame()


def get_disease_targets(disease: str, min_score: float = 0.0, progress_callback=None) -> pd.DataFrame:
    """
    Get disease targets from multiple sources, merge and deduplicate.
    Tries DisGeNET first, then GeneCards.
    """
    all_genes = []

    if progress_callback:
        progress_callback("正在查询 DisGeNET 数据库...")

    disgenet_df = get_disease_targets_disgenet(disease, min_score)
    if not disgenet_df.empty:
        all_genes.append(disgenet_df)

    if progress_callback:
        progress_callback("正在查询 GeneCards 数据库...")

    genecards_df = get_disease_targets_genecards(disease, min_score)
    if not genecards_df.empty:
        all_genes.append(genecards_df)

    if all_genes:
        merged = pd.concat(all_genes, ignore_index=True)
        return merged.drop_duplicates(subset=["Gene"]).reset_index(drop=True)

    return pd.DataFrame(columns=["Gene", "Score", "Source"])
