# space_filling_model/tests/test_segment_selection.py
import os
import sys

import numpy as np
import pytest

# Make `space_filling` importable regardless of the pytest invocation dir.
SPACE_FILLING_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if SPACE_FILLING_DIR not in sys.path:
    sys.path.insert(0, SPACE_FILLING_DIR)

from space_filling import _select_segment, run_model  # noqa: E402


def test_select_longest_is_deterministic_argmax():
    elig_idx = np.array([3, 7, 9])
    elig_lengths = np.array([0.2, 0.9, 0.5])
    rng = np.random.default_rng(0)
    assert _select_segment(rng, elig_idx, elig_lengths, 'longest') == 7


def test_select_longest_ties_break_to_first_index():
    elig_idx = np.array([4, 8])
    elig_lengths = np.array([0.5, 0.5])
    rng = np.random.default_rng(0)
    assert _select_segment(rng, elig_idx, elig_lengths, 'longest') == 4


def test_select_uniform_returns_eligible_index():
    elig_idx = np.array([2, 5, 6])
    elig_lengths = np.array([0.1, 0.2, 0.3])
    rng = np.random.default_rng(0)
    chosen = _select_segment(rng, elig_idx, elig_lengths, 'uniform')
    assert chosen in {2, 5, 6}


def test_select_power_one_matches_length_for_same_rng():
    elig_idx = np.array([1, 4, 8, 11])
    elig_lengths = np.array([0.1, 0.7, 0.3, 0.4])
    a = _select_segment(np.random.default_rng(123), elig_idx, elig_lengths, 'length')
    b = _select_segment(np.random.default_rng(123), elig_idx, elig_lengths,
                        'power', select_power=1.0)
    assert a == b


def test_select_unknown_raises():
    rng = np.random.default_rng(0)
    with pytest.raises(ValueError):
        _select_segment(rng, np.array([0]), np.array([1.0]), 'bogus')


def _seg_set(segments):
    return {(round(a[0], 6), round(a[1], 6), round(b[0], 6), round(b[1], 6))
            for a, b in segments}


def test_default_select_unchanged_vs_explicit_length():
    a = run_model(model=1, seed=42, max_iter=300, verbose=False)
    b = run_model(model=1, seed=42, max_iter=300, verbose=False, select='length')
    assert _seg_set(a) == _seg_set(b)


def test_longest_is_seed_independent():
    a = run_model(model=1, seed=1, max_iter=300, verbose=False, select='longest')
    b = run_model(model=1, seed=2, max_iter=300, verbose=False, select='longest')
    assert _seg_set(a) == _seg_set(b)


def test_uniform_differs_from_length():
    a = run_model(model=1, seed=42, max_iter=300, verbose=False, select='length')
    b = run_model(model=1, seed=42, max_iter=300, verbose=False, select='uniform')
    assert _seg_set(a) != _seg_set(b)


def test_run_model_rejects_unknown_select():
    with pytest.raises(ValueError):
        run_model(model=1, seed=42, max_iter=10, verbose=False, select='bogus')
