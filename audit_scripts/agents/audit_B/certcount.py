"""Monkeypatch certify_classical_unlink inside fast_enum to count diagrams silently dropped by failed certification."""
import collections, threading
import surfaceknots.fast_enum as fe
_orig = fe.certify_classical_unlink
COUNTS = collections.Counter()
FAILED = []
_lock = threading.Lock()
def _wrapped(d, *a, **k):
    c = _orig(d, *a, **k)
    with _lock:
        COUNTS["calls"] += 1
        if not c.is_unlink:
            COUNTS["fail:" + str(c.reason)] += 1
    return c
fe.certify_classical_unlink = _wrapped
