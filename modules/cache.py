"""
Persistent disk cache layer using diskcache.
Falls back to a no-op if diskcache is not installed.
"""

import hashlib
import json
import os
import pickle
from pathlib import Path

_CACHE_DIR = Path(__file__).parent.parent / "data" / "cache"
_TTL_SECONDS = int(os.environ.get("CACHE_TTL_DAYS", "7")) * 86400

try:
    import diskcache as dc

    _cache = dc.Cache(
        str(_CACHE_DIR),
        size_limit=int(1e9),  # 1 GB
    )

    def cache_get(key: str):
        return _cache.get(key, default=None)

    def cache_set(key: str, value, ttl: int = _TTL_SECONDS):
        _cache.set(key, value, expire=ttl)

    def cache_clear():
        _cache.clear()

    def cache_size() -> int:
        """Return number of cached entries."""
        return len(_cache)

except ImportError:
    # diskcache not installed — silent no-op

    def cache_get(key: str):  # type: ignore[misc]
        return None

    def cache_set(key: str, value, ttl: int = _TTL_SECONDS):  # type: ignore[misc]
        pass

    def cache_clear():  # type: ignore[misc]
        pass

    def cache_size() -> int:  # type: ignore[misc]
        return 0


def make_key(*args, **kwargs) -> str:
    """Create a deterministic cache key from arbitrary arguments."""
    payload = json.dumps({"args": args, "kwargs": kwargs}, sort_keys=True, default=str)
    return hashlib.sha256(payload.encode()).hexdigest()
