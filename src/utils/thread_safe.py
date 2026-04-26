"""
Thread-safe utilities for concurrent API access.

Provides rate limiting and caching primitives that are safe to use
from multiple threads (e.g. ThreadPoolExecutor workers).
"""

import threading
import time
from typing import Any, Callable, Dict, Optional, TypeVar

T = TypeVar('T')


class ThreadSafeRateLimiter:
    """
    Thread-safe rate limiter.

    Ensures that calls to ``wait()`` are spaced at least
    ``min_interval`` seconds apart, regardless of how many threads
    are calling concurrently.
    """

    def __init__(self, min_interval: float):
        self._min_interval = min_interval
        self._last_request_time = 0.0
        self._lock = threading.Lock()

    def wait(self):
        """Block the calling thread until the next request is allowed."""
        with self._lock:
            now = time.time()
            elapsed = now - self._last_request_time
            if elapsed < self._min_interval:
                time.sleep(self._min_interval - elapsed)
            self._last_request_time = time.time()


class ThreadSafeCache:
    """
    Thread-safe dictionary with a get-or-compute pattern.

    Regular ``get`` / ``set`` / ``in`` / ``clear`` operations are
    protected by a single lock.  The ``get_or_compute`` method uses
    per-key locking so that expensive computations for *different* keys
    can run in parallel, while computations for the *same* key are
    deduplicated (only the first caller computes; the rest wait).
    """

    # Class-level registry of all instances for bulk clearing
    _all_instances = []
    _instances_lock = threading.Lock()

    def __init__(self):
        self._data: Dict[str, Any] = {}
        self._lock = threading.Lock()
        self._key_locks: Dict[str, threading.Lock] = {}
        self._key_locks_lock = threading.Lock()

        # Register this instance for bulk clearing
        with ThreadSafeCache._instances_lock:
            ThreadSafeCache._all_instances.append(self)

    @classmethod
    def clear_all_instances(cls):
        """Clear all registered ThreadSafeCache instances."""
        with cls._instances_lock:
            for instance in cls._all_instances:
                instance.clear()
            cls._all_instances.clear()

    # -- basic dict-like interface ----------------------------------------

    def get(self, key: str, default: Any = None) -> Any:
        with self._lock:
            return self._data.get(key, default)

    def __contains__(self, key: str) -> bool:
        with self._lock:
            return key in self._data

    def __setitem__(self, key: str, value: Any) -> None:
        with self._lock:
            self._data[key] = value

    def __getitem__(self, key: str) -> Any:
        with self._lock:
            return self._data[key]

    def clear(self) -> None:
        with self._lock:
            self._data.clear()
        with self._key_locks_lock:
            self._key_locks.clear()

    def keys(self):
        with self._lock:
            return list(self._data.keys())

    def values(self):
        with self._lock:
            return list(self._data.values())

    def items(self):
        with self._lock:
            return list(self._data.items())

    def __len__(self) -> int:
        with self._lock:
            return len(self._data)

    # -- get-or-compute ---------------------------------------------------

    def get_or_compute(self, key: str, compute_fn: Callable[[], T]) -> T:
        """
        Return the cached value for *key*, or compute it exactly once.

        If multiple threads call this with the same key concurrently,
        only one will execute ``compute_fn``; the others will block
        until the result is available.  Different keys do not block
        each other.
        """
        # Fast path — already cached
        with self._lock:
            if key in self._data:
                return self._data[key]

        # Obtain (or create) a per-key lock
        with self._key_locks_lock:
            if key not in self._key_locks:
                self._key_locks[key] = threading.Lock()
            key_lock = self._key_locks[key]

        # Only one thread will compute for this key
        with key_lock:
            # Double-check after acquiring the key lock
            with self._lock:
                if key in self._data:
                    return self._data[key]

            value = compute_fn()

            with self._lock:
                self._data[key] = value

            return value


# ---------------------------------------------------------------------------
# Shared Ensembl rate limiter singleton
# ---------------------------------------------------------------------------
# Both EnsemblClient and VariantDBClient hit rest.ensembl.org (15 req/s).
# They must share a single rate limiter so the combined request rate stays
# within the limit.

_ensembl_rate_limiter: Optional[ThreadSafeRateLimiter] = None
_ensembl_rate_limiter_lock = threading.Lock()


def get_ensembl_rate_limiter() -> ThreadSafeRateLimiter:
    """Return (and lazily create) the shared Ensembl rate limiter."""
    global _ensembl_rate_limiter
    if _ensembl_rate_limiter is None:
        with _ensembl_rate_limiter_lock:
            if _ensembl_rate_limiter is None:
                _ensembl_rate_limiter = ThreadSafeRateLimiter(0.07)
    return _ensembl_rate_limiter
