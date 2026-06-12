# space_filling_model/tests/test_unit_square_stub_noise.py
"""Stub-noise tests for the unit-square model (space_filling.run_model)."""
import os
import sys

import pytest

# Make `space_filling` importable regardless of the pytest invocation dir.
SPACE_FILLING_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if SPACE_FILLING_DIR not in sys.path:
    sys.path.insert(0, SPACE_FILLING_DIR)

from space_filling import run_model, compute_degrees  # noqa: E402


def _seg_set(segments):
    return {(round(a[0], 6), round(a[1], 6), round(b[0], 6), round(b[1], 6))
            for a, b in segments}


def _deg1_count(segments):
    return sum(1 for d in compute_degrees(segments).values() if d == 1)


def test_pstub_zero_is_byte_identical_to_default():
    """The noise is strictly opt-in: p_stub=0 must not change anything."""
    base = run_model(model=1, seed=42, max_iter=400, verbose=False)
    zero = run_model(model=1, seed=42, max_iter=400, verbose=False, p_stub=0.0)
    assert _seg_set(base) == _seg_set(zero)


def test_pstub_is_reproducible_for_a_fixed_seed():
    a = run_model(model=1, seed=7, max_iter=400, verbose=False, p_stub=0.3)
    b = run_model(model=1, seed=7, max_iter=400, verbose=False, p_stub=0.3)
    assert _seg_set(a) == _seg_set(b)


def test_pstub_adds_segments():
    """Stubs are additive: a positive p_stub yields strictly more segments."""
    base = run_model(model=1, seed=42, max_iter=400, verbose=False, p_stub=0.0)
    noisy = run_model(model=1, seed=42, max_iter=400, verbose=False, p_stub=0.3)
    assert len(noisy) > len(base)


def test_small_snap_tol_increases_dead_ends():
    """With snap_tol=0 stubs dangle, raising the count of degree-1 nodes."""
    base = run_model(model=1, seed=42, max_iter=400, verbose=False, p_stub=0.0)
    culs = run_model(model=1, seed=42, max_iter=400, verbose=False,
                     p_stub=0.4, snap_tol=0.0)
    assert _deg1_count(culs) > _deg1_count(base)


def test_stubs_stay_inside_unit_square():
    segs = run_model(model=1, seed=3, max_iter=400, verbose=False, p_stub=0.4)
    for a, b in segs:
        for x, y in (a, b):
            assert -1e-9 <= x <= 1.0 + 1e-9
            assert -1e-9 <= y <= 1.0 + 1e-9


def test_stub_splittable_flag_changes_output():
    """The splittable knob feeds stubs back into growth, changing the network."""
    grow = run_model(model=1, seed=11, max_iter=400, verbose=False,
                     p_stub=0.4, stub_splittable=True)
    dead = run_model(model=1, seed=11, max_iter=400, verbose=False,
                     p_stub=0.4, stub_splittable=False)
    assert _seg_set(grow) != _seg_set(dead)


def test_model2_pstub_runs_and_adds_segments():
    base = run_model(model=2, seed=42, max_iter=400, verbose=False, p_stub=0.0)
    noisy = run_model(model=2, seed=42, max_iter=400, verbose=False, p_stub=0.3)
    assert len(noisy) >= len(base)


def test_stub_at_midpoint_prob_default_is_one_half():
    """Default behavior must equal an explicit stub_at_midpoint_prob=0.5."""
    a = run_model(model=1, seed=11, max_iter=400, verbose=False, p_stub=0.4)
    b = run_model(model=1, seed=11, max_iter=400, verbose=False, p_stub=0.4,
                  stub_at_midpoint_prob=0.5)
    assert _seg_set(a) == _seg_set(b)


def test_stub_at_midpoint_prob_extremes_differ():
    """Anchoring every stub at the midpoint vs at existing degree>=3 nodes
    produces different networks."""
    mid = run_model(model=1, seed=11, max_iter=400, verbose=False,
                    p_stub=0.4, stub_at_midpoint_prob=1.0)
    old = run_model(model=1, seed=11, max_iter=400, verbose=False,
                    p_stub=0.4, stub_at_midpoint_prob=0.0)
    assert _seg_set(mid) != _seg_set(old)


@pytest.mark.parametrize("bad", [-0.1, 1.5, 2.0])
def test_p_stub_out_of_range_raises(bad):
    with pytest.raises(ValueError):
        run_model(model=1, seed=1, max_iter=10, verbose=False, p_stub=bad)


@pytest.mark.parametrize("bad", [-0.1, 1.5, 2.0])
def test_stub_at_midpoint_prob_out_of_range_raises(bad):
    with pytest.raises(ValueError):
        run_model(model=1, seed=1, max_iter=10, verbose=False, p_stub=0.3,
                  stub_at_midpoint_prob=bad)


@pytest.mark.parametrize("good", [0.0, 1.0])
def test_probability_boundaries_are_accepted(good):
    # Both endpoints of [0, 1] are valid and must not raise.
    run_model(model=1, seed=1, max_iter=20, verbose=False, p_stub=good,
              stub_at_midpoint_prob=good)
