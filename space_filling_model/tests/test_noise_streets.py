# space_filling_model/tests/test_noise_streets.py
"""Tests for the post-growth low-degree noise decoration (add_noise_streets)."""
import os
import sys

import numpy as np
import pytest

# Make `space_filling` importable regardless of the pytest invocation dir.
SPACE_FILLING_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if SPACE_FILLING_DIR not in sys.path:
    sys.path.insert(0, SPACE_FILLING_DIR)

from space_filling import (  # noqa: E402
    MIN_LENGTH,
    add_noise_streets,
    compute_degrees,
    run_model,
)


def _seg_set(segments):
    return {(round(a[0], 6), round(a[1], 6), round(b[0], 6), round(b[1], 6))
            for a, b in segments}


def _deg1_count(segments):
    return sum(1 for d in compute_degrees(segments).values() if d == 1)


def _point_on_segment(p, a, b, tol=1e-9):
    """True if p lies on segment ab (within tol)."""
    ax, ay = a
    bx, by = b
    px, py = p
    dx, dy = bx - ax, by - ay
    L2 = dx * dx + dy * dy
    if L2 < tol * tol:
        return False
    t = ((px - ax) * dx + (py - ay) * dy) / L2
    if t < -tol or t > 1.0 + tol:
        return False
    qx, qy = ax + t * dx, ay + t * dy
    return float(np.hypot(px - qx, py - qy)) <= tol


def _proper_cross(p, q, a, b, tol=1e-9):
    """True if open segments pq and ab properly cross (shared endpoints OK)."""
    for u in (p, q):
        for v in (a, b):
            if abs(u[0] - v[0]) < tol and abs(u[1] - v[1]) < tol:
                return False

    def cross(o, u, v):
        return (u[0] - o[0]) * (v[1] - o[1]) - (u[1] - o[1]) * (v[0] - o[0])

    d1 = cross(p, q, a)
    d2 = cross(p, q, b)
    d3 = cross(a, b, p)
    d4 = cross(a, b, q)
    return d1 * d2 < 0 and d3 * d4 < 0


@pytest.fixture(scope='module')
def backbone():
    """A converged-enough unit-square network to decorate."""
    return run_model(model=1, seed=42, max_iter=300, verbose=False)


@pytest.fixture(scope='module')
def small_backbone():
    """Smaller network for the O(n^2) planarity scan."""
    return run_model(model=1, seed=42, max_iter=120, verbose=False)


# --------------------------------------------------------------------------- #
# Validation & no-op                                                          #
# --------------------------------------------------------------------------- #

def test_noise_frac_zero_is_a_noop(backbone):
    out = add_noise_streets(backbone, noise_frac=0.0, seed=0)
    assert out == list(backbone)


def test_empty_input_returns_empty():
    assert add_noise_streets([], noise_frac=1.0, seed=0) == []


def test_negative_noise_frac_raises(backbone):
    with pytest.raises(ValueError):
        add_noise_streets(backbone, noise_frac=-0.1, seed=0)


@pytest.mark.parametrize("bad", [-0.1, 1.5])
def test_branch_prob_out_of_range_raises(backbone, bad):
    with pytest.raises(ValueError):
        add_noise_streets(backbone, noise_frac=0.2, seed=0, branch_prob=bad)


@pytest.mark.parametrize("good", [0.0, 1.0])
def test_branch_prob_boundaries_accepted(backbone, good):
    add_noise_streets(backbone, noise_frac=0.1, seed=0, branch_prob=good)


# --------------------------------------------------------------------------- #
# Core behaviour                                                              #
# --------------------------------------------------------------------------- #

def test_determinism(backbone):
    a = add_noise_streets(backbone, noise_frac=0.5, seed=3)
    b = add_noise_streets(backbone, noise_frac=0.5, seed=3)
    assert a == b


def test_added_segment_count_bounds(backbone):
    """Each placed noise street adds 2 rows (host split + noise) or 3 (when
    its free end T-in splits another street). With this backbone/seed no
    placement should be skipped, so the row delta brackets the target."""
    out = add_noise_streets(backbone, noise_frac=0.5, seed=3)
    target = round(0.5 * len(backbone))
    added = len(out) - len(backbone)
    assert 2 * target <= added <= 3 * target


def test_endpoints_stay_inside_unit_square(backbone):
    out = add_noise_streets(backbone, noise_frac=0.8, seed=5)
    for a, b in out:
        for x, y in (a, b):
            assert -1e-9 <= x <= 1.0 + 1e-9
            assert -1e-9 <= y <= 1.0 + 1e-9


def test_planarity(small_backbone):
    """No two output segments properly cross."""
    out = add_noise_streets(small_backbone, noise_frac=0.8, seed=5)
    for i in range(len(out)):
        p, q = out[i]
        for j in range(i + 1, len(out)):
            a, b = out[j]
            assert not _proper_cross(p, q, a, b), (i, j, out[i], out[j])


def test_snap_tol_zero_yields_only_dangling_tips(backbone):
    """With snapping off every noise street dangles: exactly 2 rows per
    placement and one new degree-1 node per placement."""
    out = add_noise_streets(backbone, noise_frac=0.5, seed=3,
                            snap_tol=0.0, branch_prob=0.0)
    added = len(out) - len(backbone)
    assert added % 2 == 0
    placed = added // 2
    assert placed > 0
    assert _deg1_count(out) - _deg1_count(backbone) == placed


def test_noise_heading_respects_min_angle(backbone):
    """Backbone-anchored noise streets keep >= min_angle to their host axis
    (so DGDC at a_threshold=20deg can never merge them into the host)."""
    min_angle = np.pi / 4
    out = add_noise_streets(backbone, noise_frac=0.5, seed=3,
                            snap_tol=0.0, branch_prob=0.0,
                            min_angle=min_angle)
    in_set = _seg_set(backbone)
    checked = 0
    for a, b in out:
        if (round(a[0], 6), round(a[1], 6),
                round(b[0], 6), round(b[1], 6)) in in_set:
            continue
        # Classify: a noise street has exactly one endpoint on some input
        # segment (the anchor); a split half has both endpoints on its parent.
        hosts_a = [s for s in backbone if _point_on_segment(a, s[0], s[1])]
        hosts_b = [s for s in backbone if _point_on_segment(b, s[0], s[1])]
        if hosts_a and hosts_b:
            continue  # a split half of an input segment
        if not hosts_a and not hosts_b:
            continue  # noise anchored on noise (not checked here)
        host = hosts_a[0] if hosts_a else hosts_b[0]
        hphi = np.arctan2(host[1][1] - host[0][1], host[1][0] - host[0][0])
        nphi = np.arctan2(b[1] - a[1], b[0] - a[0])
        d = (nphi - hphi) % np.pi
        assert min(d, np.pi - d) >= min_angle - 1e-9
        checked += 1
    assert checked > 0


def test_branch_prob_one_anchors_on_noise(backbone):
    """With branch_prob=1 only the first noise street may split a backbone
    segment; everything after anchors on (and splits) noise streets."""
    out = add_noise_streets(backbone, noise_frac=0.5, seed=3,
                            snap_tol=0.0, branch_prob=1.0)
    assert len(out) > len(backbone) + 2  # more than one street was placed
    lost_input_rows = _seg_set(backbone) - _seg_set(out)
    assert len(lost_input_rows) == 1


def test_snap_creates_t_junctions(backbone):
    """With a generous snap_tol some free ends T into streets, so fewer new
    degree-1 tips appear than with snapping off."""
    dangling = add_noise_streets(backbone, noise_frac=0.5, seed=3,
                                 snap_tol=0.0, branch_prob=0.0)
    snapped = add_noise_streets(backbone, noise_frac=0.5, seed=3,
                                snap_tol=2.0 * MIN_LENGTH, branch_prob=0.0)
    assert _deg1_count(snapped) < _deg1_count(dangling)
