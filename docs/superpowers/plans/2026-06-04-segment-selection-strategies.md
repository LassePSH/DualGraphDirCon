# Segment-Selection Strategies Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the segment-to-split selection in `run_model` pluggable via a `select` string argument (`'length'` default, `'uniform'`, `'power'`, `'longest'`), keeping current behaviour unchanged by default.

**Architecture:** Extract the inline length-weighted inverse-CDF block into a module-level `_select_segment` helper that dispatches on a strategy string. `run_model` gains `select` and `select_power` keyword args (defaults reproduce today's behaviour exactly) and validates `select` once before the loop. Everything flows through `**kwargs` to `run_ensemble`/`get_geodataframe` for free.

**Tech Stack:** Python, NumPy (`np.cumsum`/`np.searchsorted`/`np.argmax`), `np.random.Generator`, pytest.

---

### Task 1: Add the `_select_segment` helper (TDD)

**Files:**
- Modify: `space_filling_model/space_filling.py` (add helper near `_find_best`, after line 304)
- Test: `space_filling_model/tests/test_segment_selection.py` (create)

- [ ] **Step 1: Write the failing test**

Create `space_filling_model/tests/test_segment_selection.py`:

```python
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd space_filling_model && python -m pytest tests/test_segment_selection.py -v`
Expected: FAIL — `ImportError: cannot import name '_select_segment'`

- [ ] **Step 3: Write minimal implementation**

In `space_filling_model/space_filling.py`, add this helper immediately after the
`_find_best` function (after line 304):

```python
def _select_segment(rng, elig_idx, elig_lengths, select='length',
                    select_power=1.0):
    """Return the chosen segment index (an int drawn from elig_idx).

    Strategies
    ----------
    'length'  : probability proportional to length (default; current behaviour).
    'uniform' : every eligible segment equally likely.
    'power'   : probability proportional to length ** select_power.
    'longest' : deterministic argmax (ties → lowest elig_idx).
    """
    if select == 'length':
        cdf = np.cumsum(elig_lengths)
        return int(elig_idx[np.searchsorted(cdf, rng.random() * cdf[-1])])
    if select == 'uniform':
        return int(rng.choice(elig_idx))
    if select == 'power':
        cdf = np.cumsum(elig_lengths ** select_power)
        return int(elig_idx[np.searchsorted(cdf, rng.random() * cdf[-1])])
    if select == 'longest':
        return int(elig_idx[np.argmax(elig_lengths)])
    raise ValueError(f"unknown select={select!r}")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cd space_filling_model && python -m pytest tests/test_segment_selection.py -v`
Expected: PASS (5 passed)

- [ ] **Step 5: Commit**

```bash
git add space_filling_model/space_filling.py space_filling_model/tests/test_segment_selection.py
git commit -m "feat: add _select_segment strategy helper"
```

---

### Task 2: Wire `select`/`select_power` into `run_model` (TDD)

**Files:**
- Modify: `space_filling_model/space_filling.py:154-156` (signature), `:240-243` (call site), and add validation before the loop (before line 219)
- Test: `space_filling_model/tests/test_segment_selection.py` (append)

- [ ] **Step 1: Write the failing test**

Append to `space_filling_model/tests/test_segment_selection.py`:

```python
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
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cd space_filling_model && python -m pytest tests/test_segment_selection.py -v`
Expected: FAIL — `run_model() got an unexpected keyword argument 'select'`

- [ ] **Step 3: Update the signature**

In `space_filling_model/space_filling.py`, change the `run_model` signature
(lines 154-156) from:

```python
def run_model(model=1, min_length=MIN_LENGTH, max_iter=200_000, seed=42,
              min_angle=np.pi / 4, max_angle=None, area_coeff=0.05,
              area_scale=1.0, verbose=True):
```

to:

```python
def run_model(model=1, min_length=MIN_LENGTH, max_iter=200_000, seed=42,
              min_angle=np.pi / 4, max_angle=None, area_coeff=0.05,
              area_scale=1.0, select='length', select_power=1.0, verbose=True):
```

- [ ] **Step 4: Add fail-fast validation before the loop**

In `space_filling_model/space_filling.py`, immediately before
`converged = False` (line 218), insert:

```python
    if select not in ('length', 'uniform', 'power', 'longest'):
        raise ValueError(f"unknown select={select!r}")

```

- [ ] **Step 5: Replace the inline selection block with the helper**

In `space_filling_model/space_filling.py`, replace lines 240-243:

```python
        # Length-weighted sampling via inverse-CDF (faster than np.random.choice).
        elig_lengths = lengths[elig_idx]
        cdf = np.cumsum(elig_lengths)
        seg_i = int(elig_idx[np.searchsorted(cdf, rng.random() * cdf[-1])])
```

with:

```python
        elig_lengths = lengths[elig_idx]
        seg_i = _select_segment(rng, elig_idx, elig_lengths, select,
                                select_power)
```

- [ ] **Step 6: Run tests to verify they pass**

Run: `cd space_filling_model && python -m pytest tests/test_segment_selection.py -v`
Expected: PASS (9 passed)

- [ ] **Step 7: Commit**

```bash
git add space_filling_model/space_filling.py space_filling_model/tests/test_segment_selection.py
git commit -m "feat: expose select/select_power on run_model"
```

---

### Task 3: Document the new strategies

**Files:**
- Modify: `space_filling_model/space_filling.py` (module docstring lines 2-21, and `run_model` docstring lines 157-173)

- [ ] **Step 1: Update the `run_model` docstring**

In `space_filling_model/space_filling.py`, in the Parameters block of
`run_model` (after the `area_scale` entry, before `verbose`), add:

```python
    select     : segment-selection strategy. 'length' (default) picks with
                 probability proportional to length; 'uniform' picks any
                 eligible segment equally; 'power' uses length**select_power;
                 'longest' deterministically splits the longest eligible segment.
    select_power : exponent for select='power' (1.0 ≡ 'length', 0.0 ≡ 'uniform').
```

- [ ] **Step 2: Update the module docstring**

In `space_filling_model/space_filling.py`, after the "Both start from a unit
square..." paragraph (after line 14), add:

```
Segment selection is configurable via ``select`` ('length' default, 'uniform',
'power' with ``select_power``, or 'longest'); see ``run_model``.
```

- [ ] **Step 3: Verify the module still imports and runs**

Run: `cd space_filling_model && python -c "from space_filling import run_model; print(len(run_model(model=1, seed=42, max_iter=200, verbose=False, select='longest')))"`
Expected: prints an integer segment count (no error)

- [ ] **Step 4: Run the full selection test file once more**

Run: `cd space_filling_model && python -m pytest tests/test_segment_selection.py -v`
Expected: PASS (9 passed)

- [ ] **Step 5: Commit**

```bash
git add space_filling_model/space_filling.py
git commit -m "docs: document select/select_power strategies"
```

---

## Notes

- `real_space_filling.py` is intentionally **out of scope** (it duplicates the
  selection block in its own run loop). Per the design decision on 2026-06-04 it
  is left unchanged.
- No CLI or plot-script wiring is needed; callers pass `select=`/`select_power=`
  directly to `run_model`, `run_ensemble`, or `get_geodataframe`.
