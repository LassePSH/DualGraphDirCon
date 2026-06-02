# space_filling_model/tests/test_stub_noise.py
import os
import sys
import hashlib
from collections import Counter

import numpy as np
import pytest

# Make `real_space_filling` importable regardless of the pytest invocation dir.
SPACE_FILLING_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if SPACE_FILLING_DIR not in sys.path:
    sys.path.insert(0, SPACE_FILLING_DIR)

from real_space_filling import run_circular_model  # noqa: E402

# Fixed, fast, deterministic configuration reused across tests.
BASE_KW = dict(area_m2=np.pi * 1000.0 ** 2, model=2, min_length=50.0,
               max_iter=2000, seed=42, verbose=False)


def _geom_hash(segments):
    arr = np.round(np.array([[a[0], a[1], b[0], b[1]] for a, b in segments]), 3)
    return hashlib.sha256(np.ascontiguousarray(arr).tobytes()).hexdigest()[:16]


def _degree_counter(segments):
    c = Counter()
    for a, b in segments:
        c[(round(a[0], 6), round(a[1], 6))] += 1
        c[(round(b[0], 6), round(b[1], 6))] += 1
    return c


def test_baseline_geometry_unchanged():
    """Current default behavior must stay byte-stable through the refactor."""
    res = run_circular_model(**BASE_KW)
    segs = res['segments']
    assert len(segs) == 3728
    assert _geom_hash(segs) == '913fac0a0326c33b'


def test_nearest_segment_snaps_within_tol():
    from real_space_filling import _nearest_segment
    # One horizontal segment from (0,0) to (10,0).
    seg = np.array([[[0.0, 0.0], [10.0, 0.0]]])
    # A point just above the middle, within tolerance.
    i, fx, fy = _nearest_segment(5.0, 1.0, seg, max_dist=2.0)
    assert i == 0
    assert fx == pytest.approx(5.0)
    assert fy == pytest.approx(0.0)


def test_nearest_segment_returns_none_outside_tol():
    from real_space_filling import _nearest_segment
    seg = np.array([[[0.0, 0.0], [10.0, 0.0]]])
    i, fx, fy = _nearest_segment(5.0, 5.0, seg, max_dist=2.0)
    assert i == -1
    assert (fx, fy) == (5.0, 5.0)


def test_nearest_segment_excludes_anchor_edges():
    from real_space_filling import _nearest_segment
    # Segment incident to the anchor (0,0) must be ignored even if nearest.
    seg = np.array([[[0.0, 0.0], [10.0, 0.0]]])
    i, fx, fy = _nearest_segment(5.0, 1.0, seg, max_dist=2.0,
                                 exclude_xy=(0.0, 0.0))
    assert i == -1
