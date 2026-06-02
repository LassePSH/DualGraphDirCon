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


def test_p_stub_zero_is_identical():
    """p_stub=0.0 must reproduce the baseline geometry exactly."""
    res = run_circular_model(p_stub=0.0, **BASE_KW)
    assert len(res['segments']) == 3728
    assert _geom_hash(res['segments']) == '913fac0a0326c33b'


def test_stubs_add_low_degree_deadends():
    """Turning on stubs must raise the count of degree-1 dead-end nodes."""
    base = run_circular_model(p_stub=0.0, **BASE_KW)['segments']
    noisy = run_circular_model(p_stub=0.3, **BASE_KW)['segments']

    base_deg1 = sum(1 for v in _degree_counter(base).values() if v == 1)
    noisy_deg1 = sum(1 for v in _degree_counter(noisy).values() if v == 1)

    assert noisy_deg1 > base_deg1


def test_stub_len_max_bounds_new_deadends():
    """Tightening stub_len_max must not create short dead-ends longer than the
    bound. We compare against the baseline so pre-existing long degree-1
    endpoints (boundary spokes) are excluded."""
    stub_len_max = 150.0

    def short_deadend_lengths(segments):
        counter = _degree_counter(segments)
        out = []
        for a, b in segments:
            ka = (round(a[0], 6), round(a[1], 6))
            kb = (round(b[0], 6), round(b[1], 6))
            if counter[ka] == 1 or counter[kb] == 1:
                out.append(float(np.hypot(b[0] - a[0], b[1] - a[1])))
        return out

    base = short_deadend_lengths(run_circular_model(p_stub=0.0, **BASE_KW)['segments'])
    noisy = short_deadend_lengths(
        run_circular_model(p_stub=0.3, stub_len_max=stub_len_max, **BASE_KW)['segments'])

    base_short = sum(1 for L in base if L <= stub_len_max + 1e-6)
    noisy_short = sum(1 for L in noisy if L <= stub_len_max + 1e-6)
    assert noisy_short > base_short
    for L in noisy:
        if L <= stub_len_max + 1e-6:
            continue  # within bound
        # Any dead-end longer than the bound must already exist at baseline
        # (e.g. boundary spokes), not be a freshly minted stub.
        assert any(abs(L - Lb) < 1e-6 for Lb in base)


def test_reproducible_with_seed():
    a = run_circular_model(p_stub=0.3, **BASE_KW)['segments']
    b = run_circular_model(p_stub=0.3, **BASE_KW)['segments']
    assert _geom_hash(a) == _geom_hash(b)
