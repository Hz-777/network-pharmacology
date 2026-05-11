"""
Molecular docking integration module.

Flow:
  1. Auto-fetch best PDB from RCSB for given gene symbols
  2. Submit compound SMILES + protein PDB to CB-Dock2 web service
  3. Return a score matrix (compound × target) for heatmap visualization

CB-Dock2 references:
  Main:   https://cadd.labshare.cn/cb-dock2/
  Backup: http://183.56.231.194:8001/cb-dock2
"""

import io
import time
import json
import requests
import tempfile
import os
from pathlib import Path
from typing import Optional
import pandas as pd
import numpy as np

from modules.cache import cache_get, cache_set, make_key

RCSB_SEARCH  = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_ENTRY   = "https://data.rcsb.org/rest/v1/core/entry/{}"
RCSB_DOWNLOAD = "https://files.rcsb.org/download/{}.pdb"

CBDOCK2_MAIN   = "https://cadd.labshare.cn/cb-dock2"
CBDOCK2_BACKUP = "http://183.56.231.194:8001/cb-dock2"

HEADERS = {"User-Agent": "Mozilla/5.0 (network-pharmacology-app/1.0)"}
TIMEOUT = 30
UPLOAD_TIMEOUT = (12, 120)   # (connect, read) for large PDB file uploads
SUBMIT_TIMEOUT = (12, 45)
POLL_TIMEOUT   = (12, 30)


# ── RCSB PDB ─────────────────────────────────────────────────────────────────



def fetch_pdb_for_gene(gene_symbol: str) -> Optional[tuple]:
    """
    Fetch PDB file content for a gene symbol.
    Returns (pdb_id, pdb_text) or None if not found.
    Uses disk cache (30d); failed lookups are not cached.
    """
    key = make_key("pdb_for_gene", gene_symbol)
    cached = cache_get(key)
    if cached is not None:  # only truthy (pdb_id, content) tuples are stored
        return cached

    # Get up to 3 candidate PDB IDs and try each until download succeeds
    candidates = _get_pdb_candidates(gene_symbol, top_n=3)
    for pdb_id in candidates:
        try:
            resp = requests.get(RCSB_DOWNLOAD.format(pdb_id), timeout=60)
            if resp.status_code == 200 and resp.text.strip():
                content = resp.text
                cache_set(key, (pdb_id, content), category="chemical")
                return pdb_id, content
        except Exception:
            continue

    return None


def _get_pdb_candidates(gene_symbol: str, top_n: int = 3) -> list:
    """Return up to top_n PDB IDs for a gene from RCSB search."""
    for q in [
        {
            "query": {
                "type": "group",
                "logical_operator": "and",
                "nodes": [
                    {"type": "terminal", "service": "full_text",
                     "parameters": {"value": gene_symbol}},
                    {"type": "terminal", "service": "text", "parameters": {
                        "attribute": "rcsb_entity_source_organism.scientific_name",
                        "operator": "exact_match", "value": "Homo sapiens"}},
                ],
            },
            "request_options": {"paginate": {"start": 0, "rows": top_n},
                                 "sort": [{"sort_by": "score", "direction": "desc"}]},
            "return_type": "entry",
        },
        {
            "query": {"type": "terminal", "service": "full_text",
                      "parameters": {"value": f"{gene_symbol} homo sapiens"}},
            "request_options": {"paginate": {"start": 0, "rows": top_n},
                                 "sort": [{"sort_by": "score", "direction": "desc"}]},
            "return_type": "entry",
        },
    ]:
        try:
            resp = requests.post(RCSB_SEARCH, json=q, timeout=TIMEOUT)
            if resp.status_code == 200:
                hits = resp.json().get("result_set", [])
                if hits:
                    return [h["identifier"] for h in hits[:top_n]]
        except Exception:
            continue
    return []


def fetch_pdb_batch(gene_list: list, max_workers: int = 3,
                    progress_callback=None) -> dict:
    """
    Fetch PDB files for multiple genes in parallel.
    Returns dict: {gene: (pdb_id, pdb_content)} or {gene: None}
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    results = {}

    def _fetch(gene):
        result = fetch_pdb_for_gene(gene)
        if progress_callback:
            progress_callback(f"PDB: {gene} → {'找到' if result else '未找到'}")
        return gene, result

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_fetch, g): g for g in gene_list}
        for future in as_completed(futures):
            gene, result = future.result()
            results[gene] = result

    return results


# ── CB-Dock2 ──────────────────────────────────────────────────────────────────

def _probe_cbdock2(base: str) -> bool:
    """
    Probe whether a CB-Dock2 base URL is truly serving its API
    (not just the landing page). Tries the upload endpoint with a HEAD/GET.
    """
    try:
        r = requests.get(
            f"{base}/php/upload_protein.php",
            timeout=(8, 10),
            allow_redirects=False,
        )
        # 200 (endpoint exists) or 405 (method not allowed = endpoint exists)
        # are both good. 301/302 to maintenance page is bad.
        if r.status_code in (200, 405):
            return True
        if r.status_code in (301, 302):
            loc = r.headers.get("Location", "")
            return "maintenance" not in loc.lower()
        # 404/5xx → endpoint missing or server error
        return False
    except Exception:
        return False


def _get_cbdock2_base() -> Optional[str]:
    """Return the first CB-Dock2 endpoint whose API is actually responding."""
    for base in (CBDOCK2_MAIN, CBDOCK2_BACKUP):
        if _probe_cbdock2(base):
            return base
    return None


def _upload_pdb(base_url: str, pdb_path: str) -> Optional[str]:
    """Upload PDB file to CB-Dock2. Returns protein_id or None."""
    upload_url = f"{base_url}/php/upload_protein.php"
    last_exc = None
    for attempt in range(3):
        try:
            with open(pdb_path, "rb") as pf:
                resp = requests.post(
                    upload_url,
                    files={"protein_file": ("protein.pdb", pf, "chemical/x-pdb")},
                    headers=HEADERS,
                    timeout=UPLOAD_TIMEOUT,
                )
            if resp.status_code != 200:
                return None
            data = resp.json()
            pid = data.get("protein_id") or data.get("id")
            if pid:
                return str(pid)
            return None
        except (requests.Timeout, requests.ConnectionError) as e:
            last_exc = e
            if attempt < 2:
                time.sleep(3 * (attempt + 1))
            continue
        except Exception:
            return None
    raise last_exc  # propagate so caller can try backup URL


def _submit_cbdock2(smiles: str, pdb_content: str,
                    base_url: str) -> Optional[float]:
    """
    Submit one SMILES + PDB to CB-Dock2. Returns binding energy (kcal/mol) or None.
    Raises requests.RequestException on network-level failure (caller may retry with backup).
    """
    with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False, mode="w") as f:
        f.write(pdb_content)
        pdb_path = f.name

    try:
        protein_id = _upload_pdb(base_url, pdb_path)
        if not protein_id:
            return None

        # Submit docking job
        dock_url = f"{base_url}/php/dock.php"
        dock_resp = requests.post(
            dock_url,
            data={"smiles": smiles, "protein_id": protein_id},
            headers=HEADERS,
            timeout=SUBMIT_TIMEOUT,
        )
        if dock_resp.status_code != 200:
            return None
        dock_data = dock_resp.json()
        job_id = dock_data.get("job_id") or dock_data.get("id")
        if not job_id:
            return None

        # Poll for result (up to ~150 s)
        result_url = f"{base_url}/php/get_result.php"
        for _ in range(30):
            time.sleep(5)
            try:
                r = requests.get(
                    result_url, params={"job_id": job_id},
                    headers=HEADERS, timeout=POLL_TIMEOUT,
                )
            except Exception:
                continue
            if r.status_code == 200:
                data = r.json()
                status = data.get("status", "")
                if status == "done":
                    score = data.get("vina_score") or data.get("score")
                    if score is not None:
                        return float(score)
                    results_list = data.get("results") or []
                    if results_list:
                        return float(results_list[0].get("score", 0))
                elif status in ("failed", "error"):
                    return None
        return None

    finally:
        try:
            os.unlink(pdb_path)
        except Exception:
            pass


def dock_compound_target(smiles: str, pdb_id: str, pdb_content: str,
                         compound_name: str = "", gene: str = "") -> Optional[float]:
    """
    Dock one compound against one protein. Returns binding energy or None.
    Tries main URL first, falls back to backup on network errors. Uses disk cache (7d).
    """
    key = make_key("docking", smiles, pdb_id)
    cached = cache_get(key)
    if cached is not None:
        return cached

    for base in (CBDOCK2_MAIN, CBDOCK2_BACKUP):
        try:
            score = _submit_cbdock2(smiles, pdb_content, base)
            if score is not None:
                cache_set(key, score, category="target")
            return score
        except (requests.Timeout, requests.ConnectionError):
            continue
        except Exception:
            return None

    return None


def run_docking_matrix(
    compounds: list,            # list of {"name": str, "SMILES": str}
    targets: list,              # list of {"gene": str, "pdb_id": str, "pdb_content": str}
    progress_callback=None,
) -> pd.DataFrame:
    """
    Run all-vs-all docking. Returns DataFrame rows=compounds, cols=genes.
    Scores are binding energies (kcal/mol, negative = better binding).
    """
    rows = {}
    total = len(compounds) * len(targets)
    done = 0

    for comp in compounds:
        name = comp["name"]
        smiles = comp["SMILES"]
        row = {}
        for tgt in targets:
            gene = tgt["gene"]
            pdb_id = tgt["pdb_id"]
            pdb_content = tgt["pdb_content"]
            score = dock_compound_target(smiles, pdb_id, pdb_content, name, gene)
            row[gene] = score
            done += 1
            if progress_callback:
                progress_callback(
                    f"对接进度 {done}/{total}: {name} × {gene} → "
                    f"{f'{score:.2f} kcal/mol' if score is not None else '失败'}"
                )
        rows[name] = row

    df = pd.DataFrame(rows).T
    df.index.name = "Compound"
    return df


# ── Visualization ─────────────────────────────────────────────────────────────

def plot_docking_heatmap(score_df: pd.DataFrame) -> bytes:
    """
    Plot compound × target docking score heatmap.
    Returns PNG bytes.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import matplotlib.font_manager as fm

    _avail = {f.name for f in fm.fontManager.ttflist}
    _CJK = ["SimHei", "STHeiti", "Heiti TC", "Microsoft YaHei",
            "WenQuanYi Zen Hei", "WenQuanYi Micro Hei",
            "Noto Sans CJK SC", "Noto Sans SC", "PingFang SC", "Arial Unicode MS"]
    _cjk_avail = [f for f in _CJK if f in _avail]
    plt.rcParams.update({
        "font.family": _cjk_avail + ["DejaVu Sans"] if _cjk_avail else ["DejaVu Sans"],
        "axes.unicode_minus": False,
    })

    if score_df.empty:
        fig, ax = plt.subplots(figsize=(6, 3))
        ax.text(0.5, 0.5, "无对接结果", ha="center", va="center", fontsize=14)
        ax.axis("off")
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
        plt.close(fig)
        return buf.getvalue()

    data = score_df.copy().astype(float)
    n_comp, n_gene = data.shape
    fig_w = max(8, n_gene * 1.2)
    fig_h = max(4, n_comp * 0.8)

    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    # Color: lower (more negative) = better binding → deeper blue/purple
    vmin = data.min().min()
    vmax = min(0, data.max().max())
    norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=(vmin + vmax) / 2, vmax=max(vmax, vmin + 0.01))
    cmap = plt.cm.RdYlBu_r

    im = ax.imshow(data.values, aspect="auto", cmap=cmap,
                   norm=norm, interpolation="nearest")

    ax.set_xticks(range(n_gene))
    ax.set_yticks(range(n_comp))
    ax.set_xticklabels(data.columns, rotation=45, ha="right", fontsize=9)
    ax.set_yticklabels(data.index, fontsize=9)
    ax.set_xlabel("靶蛋白 (基因名)", fontsize=11)
    ax.set_ylabel("活性化合物", fontsize=11)
    ax.set_title("分子对接评分热图 (kcal/mol)", fontsize=13, fontweight="bold", pad=12)

    # Annotate cells
    for i in range(n_comp):
        for j in range(n_gene):
            val = data.values[i, j]
            if not np.isnan(val):
                ax.text(j, i, f"{val:.1f}", ha="center", va="center",
                        fontsize=7, color="white" if val < (vmin + vmax) / 2 else "black")
            else:
                ax.text(j, i, "N/A", ha="center", va="center",
                        fontsize=7, color="#888888")

    cbar = fig.colorbar(im, ax=ax, shrink=0.7, pad=0.02)
    cbar.set_label("结合能 (kcal/mol)\n← 值越小结合越强", fontsize=9)

    fig.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    return buf.getvalue()


def check_cbdock2_status() -> dict:
    """Check if CB-Dock2 API is reachable (probes the upload endpoint). Returns status dict."""
    for label, url in [("主站", CBDOCK2_MAIN), ("备用IP", CBDOCK2_BACKUP)]:
        if _probe_cbdock2(url):
            return {"available": True, "url": url, "label": label}
    return {"available": False, "url": None, "label": None}
