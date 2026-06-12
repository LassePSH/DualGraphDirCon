# Segment-selection strategies for the space-filling model

**Date:** 2026-06-04
**Status:** Approved (design)
**Scope:** `space_filling_model/space_filling.py` only

## Problem

`run_model` always picks the segment to split with probability proportional to
its length (inverse-CDF sampling at `space_filling.py:240-243`). This is
hard-wired with no way to try other selection rules. We want to make the
selection strategy pluggable while keeping the current behaviour as the default,
so existing results and callers are unchanged.

## Strategies

Selected via a string `select` argument:

| `select`    | Behaviour                                                            | Stochastic? |
|-------------|---------------------------------------------------------------------|-------------|
| `'length'`  | Probability proportional to length (current default). **Unchanged.**| yes         |
| `'uniform'` | Every eligible segment equally likely.                              | yes         |
| `'power'`   | Probability proportional to `length ** select_power`.              | yes         |
| `'longest'` | Deterministic: always the longest eligible segment (greedy argmax). | no          |

Notes:
- `'power'` generalises the others: `select_power=1.0` ≡ `'length'`,
  `select_power=0.0` ≡ `'uniform'`, `select_power>1` favours long segments more
  aggressively.
- `'length'` vs `'longest'`: `'length'` is weighted-random (any segment can be
  picked, long ones more often); `'longest'` is pure argmax (always the single
  longest). They produce visually different networks.
- `'longest'` ties are broken by first index (lowest `elig_idx`), making it fully
  deterministic and seed-independent.

## API

`run_model` gains two keyword arguments:

```python
def run_model(model=1, min_length=MIN_LENGTH, max_iter=200_000, seed=42,
              min_angle=np.pi/4, max_angle=None, area_coeff=0.05,
              area_scale=1.0, select='length', select_power=1.0, verbose=True):
```

Both default to values that reproduce current behaviour exactly. They flow
through `**kwargs` into `run_ensemble`/`_ensemble_worker` and
`get_geodataframe` automatically (no signature changes needed there).

## Implementation

New module-level helper replacing the inline block at `space_filling.py:240-243`:

```python
def _select_segment(rng, elig_idx, elig_lengths, select='length', select_power=1.0):
    """Return the chosen segment index (an int from elig_idx)."""
    if select == 'length':
        cdf = np.cumsum(elig_lengths)
        return int(elig_idx[np.searchsorted(cdf, rng.random() * cdf[-1])])
    if select == 'uniform':
        return int(rng.choice(elig_idx))
    if select == 'power':
        w = elig_lengths ** select_power
        cdf = np.cumsum(w)
        return int(elig_idx[np.searchsorted(cdf, rng.random() * cdf[-1])])
    if select == 'longest':
        return int(elig_idx[np.argmax(elig_lengths)])
    raise ValueError(f"unknown select={select!r}")
```

The hot-loop call site becomes:

```python
elig_lengths = lengths[elig_idx]
seg_i = _select_segment(rng, elig_idx, elig_lengths, select, select_power)
```

Validation of `select` happens once up front in `run_model` (raise `ValueError`
before the loop) so a bad value fails fast rather than on the first iteration.

## Out of scope

- `real_space_filling.py` (duplicates the same selection block) — left unchanged
  for now, per decision on 2026-06-04.
- Radius-biased selection — not requested.
- Callable/custom strategy objects — string dispatch only.
- CLI / plot-script wiring — not requested; callers can pass `select=` directly.

## Testing

- `select='length'` with a fixed seed produces byte-identical output to the
  pre-change code (regression guard).
- `select='longest'` is deterministic across two different seeds.
- `select='power', select_power=1.0` matches `select='length'` for the same seed.
- `select='power', select_power=0.0` is distributionally uniform (all weights
  equal). Note it consumes the RNG via the inverse-CDF path, so it is not
  byte-identical to `select='uniform'` (which uses `rng.choice`); only the
  distribution matches.
- Unknown `select` raises `ValueError`.
