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


def _deadend_seg_lengths(segments):
    """Lengths of every segment that has at least one degree-1 (dead-end)
    endpoint."""
    c = _degree_counter(segments)
    out = []
    for a, b in segments:
        ka = (round(a[0], 6), round(a[1], 6))
        kb = (round(b[0], 6), round(b[1], 6))
        if c[ka] == 1 or c[kb] == 1:
            out.append(float(np.hypot(b[0] - a[0], b[1] - a[1])))
    return out


def test_stubs_are_additive():
    """Stubs are an additive layer, so enabling them adds segments, and more
    segments accrue as p_stub rises."""
    n0 = len(run_circular_model(p_stub=0.0, **BASE_KW)['segments'])
    n_lo = len(run_circular_model(p_stub=0.1, **BASE_KW)['segments'])
    n_hi = len(run_circular_model(p_stub=0.3, **BASE_KW)['segments'])
    assert n0 < n_lo < n_hi


def test_stubs_dangle_as_culdesacs_when_snap_tol_small():
    """With a tiny snap_tol a stub free end cannot T into a nearby street, so it
    dangles as a degree-1 cul-de-sac and raises the degree-1 node count."""
    base = run_circular_model(p_stub=0.0, **BASE_KW)['segments']
    noisy = run_circular_model(p_stub=0.3, snap_tol=1.0, **BASE_KW)['segments']

    base_deg1 = sum(1 for v in _degree_counter(base).values() if v == 1)
    noisy_deg1 = sum(1 for v in _degree_counter(noisy).values() if v == 1)
    assert noisy_deg1 > base_deg1


def test_stubs_tee_in_as_degree2_with_default_snap_tol():
    """In a dense domain (default snap_tol) stub free ends snap onto nearby
    streets, creating degree-2 junctions — populating the low-degree shoulder."""
    base = run_circular_model(p_stub=0.0, **BASE_KW)['segments']
    noisy = run_circular_model(p_stub=0.3, **BASE_KW)['segments']

    base_deg2 = sum(1 for v in _degree_counter(base).values() if v == 2)
    noisy_deg2 = sum(1 for v in _degree_counter(noisy).values() if v == 2)
    assert noisy_deg2 > base_deg2


def test_stub_length_respects_stub_len_max():
    """Dangling stubs (forced via a tiny snap_tol) never exceed stub_len_max,
    and some are longer than the baseline's longest dead-end — confirming the
    bound scales stub length rather than being coincidentally satisfied."""
    stub_len_max = 200.0
    base = run_circular_model(p_stub=0.0, **BASE_KW)['segments']
    noisy = run_circular_model(p_stub=0.3, snap_tol=1.0,
                               stub_len_max=stub_len_max, **BASE_KW)['segments']

    base_max = max(_deadend_seg_lengths(base))
    noisy_max = max(_deadend_seg_lengths(noisy))
    assert noisy_max <= stub_len_max + 1e-6
    assert noisy_max > base_max + 1e-6


def test_reproducible_with_seed():
    a = run_circular_model(p_stub=0.3, **BASE_KW)['segments']
    b = run_circular_model(p_stub=0.3, **BASE_KW)['segments']
    assert _geom_hash(a) == _geom_hash(b)
