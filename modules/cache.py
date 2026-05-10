"""
Persistent disk cache layer using diskcache.
Falls back to a no-op if diskcache is not installed or filesystem is read-only.

TTL 按数据类型分级：
  - chemical  (PubChem/SMILES)  : 30 天  — 化学结构不会变
  - compound  (TCMSP/HERB化合物): 30 天  — 数据库更新慢
  - target    (ChEMBL/HERB靶点) :  7 天
  - disease   (Open Targets)    : 14 天
  - ppi       (STRING)          : 14 天
  - enrichment(Enrichr)         :  7 天
"""

import hashlib
import json
import os
from pathlib import Path

# ── 分级 TTL（秒）────────────────────────────────────────────────────────────
TTL = {
    "chemical":   30 * 86400,
    "compound":   30 * 86400,
    "target":      7 * 86400,
    "disease":    14 * 86400,
    "ppi":        14 * 86400,
    "enrichment":  7 * 86400,
    "default":     7 * 86400,
}

_CACHE_DIR = Path(__file__).parent.parent / "data" / "cache"

# ── 尝试初始化 diskcache ─────────────────────────────────────────────────────
_cache = None
_enabled = False

try:
    import diskcache as dc

    _CACHE_DIR.mkdir(parents=True, exist_ok=True)
    _cache = dc.Cache(str(_CACHE_DIR), size_limit=int(1e9))  # 1 GB 上限
    _enabled = True

except Exception:
    pass  # diskcache 未安装或文件系统只读 → 全部 no-op


# ── 公共 API ─────────────────────────────────────────────────────────────────

def cache_get(key: str):
    """Return cached value, or None on miss / disabled."""
    if not _enabled:
        return None
    try:
        return _cache.get(key, default=None)
    except Exception:
        return None


def cache_set(key: str, value, category: str = "default"):
    """Store value with TTL determined by category."""
    if not _enabled:
        return
    try:
        ttl = TTL.get(category, TTL["default"])
        _cache.set(key, value, expire=ttl)
    except Exception:
        pass


def cache_clear():
    """Delete all cached entries."""
    if not _enabled:
        return
    try:
        _cache.clear()
    except Exception:
        pass


def cache_info() -> dict:
    """Return cache statistics for display in the UI."""
    if not _enabled:
        return {"enabled": False, "items": 0, "size_mb": 0.0}
    try:
        return {
            "enabled": True,
            "items": len(_cache),
            "size_mb": round(_cache.volume() / 1e6, 1),
        }
    except Exception:
        return {"enabled": True, "items": 0, "size_mb": 0.0}


def make_key(*args, **kwargs) -> str:
    """Create a deterministic cache key from arbitrary arguments."""
    payload = json.dumps({"args": args, "kwargs": kwargs}, sort_keys=True, default=str)
    return hashlib.sha256(payload.encode()).hexdigest()
