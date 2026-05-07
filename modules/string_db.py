"""STRING database API wrapper for PPI network analysis."""

import requests
import pandas as pd
import io
from typing import Optional

STRING_API = "https://string-db.org/api"
SPECIES_HUMAN = 9606


def get_ppi_network(genes: list, min_score: int = 400) -> pd.DataFrame:
    """
    Fetch protein-protein interaction network from STRING database.
    min_score: interaction score threshold (0-1000), 400 = medium confidence
    Returns DataFrame with columns: preferredName_A, preferredName_B, score
    """
    if not genes:
        return pd.DataFrame()

    # First map gene names to STRING IDs
    genes_str = "%0d".join(genes)
    map_url = f"{STRING_API}/json/get_string_ids"
    params = {
        "identifiers": genes_str,
        "species": SPECIES_HUMAN,
        "limit": 1,
        "echo_query": 1,
        "caller_identity": "network_pharmacology_tool",
    }

    try:
        resp = requests.post(map_url, data=params, timeout=60)
        resp.raise_for_status()
    except Exception as e:
        raise ConnectionError(f"STRING数据库连接失败: {e}")

    mapped = resp.json()
    if not mapped:
        return pd.DataFrame()

    string_ids = [item["stringId"] for item in mapped if "stringId" in item]
    if not string_ids:
        return pd.DataFrame()

    # Get interaction network
    network_url = f"{STRING_API}/tsv/network"
    network_params = {
        "identifiers": "%0d".join(string_ids),
        "species": SPECIES_HUMAN,
        "required_score": min_score,
        "network_type": "functional",
        "caller_identity": "network_pharmacology_tool",
    }

    try:
        resp = requests.post(network_url, data=network_params, timeout=60)
        resp.raise_for_status()
        df = pd.read_csv(io.StringIO(resp.text), sep="\t")
        return df
    except Exception as e:
        raise ConnectionError(f"STRING网络获取失败: {e}")


def get_network_enrichment(genes: list) -> dict:
    """Get functional enrichment from STRING."""
    if not genes:
        return {}

    genes_str = "%0d".join(genes)
    url = f"{STRING_API}/json/enrichment"
    params = {
        "identifiers": genes_str,
        "species": SPECIES_HUMAN,
        "caller_identity": "network_pharmacology_tool",
    }

    try:
        resp = requests.post(url, data=params, timeout=60)
        resp.raise_for_status()
        return resp.json()
    except Exception:
        return {}


def calculate_network_centrality(ppi_df: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate degree centrality for each node in PPI network.
    Returns DataFrame sorted by degree (hub genes first).
    """
    import networkx as nx

    if ppi_df.empty:
        return pd.DataFrame()

    # Find column names
    node_a_col = next((c for c in ppi_df.columns if "preferredName_A" in c or "node1" in c.lower()), None)
    node_b_col = next((c for c in ppi_df.columns if "preferredName_B" in c or "node2" in c.lower()), None)
    score_col = next((c for c in ppi_df.columns if "score" in c.lower()), None)

    if not node_a_col or not node_b_col:
        # Try first two columns
        cols = ppi_df.columns.tolist()
        node_a_col, node_b_col = cols[0], cols[1]

    G = nx.from_pandas_edgelist(
        ppi_df,
        source=node_a_col,
        target=node_b_col,
        edge_attr=score_col if score_col else None,
    )

    degree = dict(G.degree())
    betweenness = nx.betweenness_centrality(G)
    closeness = nx.closeness_centrality(G)

    centrality_df = pd.DataFrame({
        "Gene": list(degree.keys()),
        "Degree": list(degree.values()),
        "Betweenness": [betweenness.get(n, 0) for n in degree.keys()],
        "Closeness": [closeness.get(n, 0) for n in degree.keys()],
    })

    centrality_df["HubScore"] = (
        centrality_df["Degree"] / centrality_df["Degree"].max() * 0.4 +
        centrality_df["Betweenness"] / max(centrality_df["Betweenness"].max(), 1e-10) * 0.4 +
        centrality_df["Closeness"] / max(centrality_df["Closeness"].max(), 1e-10) * 0.2
    )

    return centrality_df.sort_values("HubScore", ascending=False).reset_index(drop=True)


def get_top_hub_genes(centrality_df: pd.DataFrame, top_n: int = 20) -> list:
    """Return top hub gene names."""
    if centrality_df.empty:
        return []
    return centrality_df.head(top_n)["Gene"].tolist()
