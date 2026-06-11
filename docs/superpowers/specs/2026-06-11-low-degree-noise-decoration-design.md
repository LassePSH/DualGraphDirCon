# Post-growth low-degree noise decoration for the unit-square space-filling model

**Date:** 2026-06-11
**Component:** `space_filling_model/space_filling.py`
**Status:** Approved (design), pending implementation

## Problem

The space-filling model reproduces the scale-free high-degree (k > 4) regime of
the DGDC dual-graph degree distribution, but under-produces degree 1–3 streets.
Empirically the distribution should be approximately constant up to the
transition at φ = 4. The existing in-loop stub mechanism (`p_stub`) cannot fix
this even when maxed out, for structural reasons:

- It is Bernoulli per growth iteration, so stub streets are capped at ~one per
  iteration (≈50% of street additions), minus rejections.
- Stubs anchor only at split midpoints or degree-≥3 junctions, so they touch
  2+ streets at the anchor and mostly land at dual degree 2–3, rarely 1.
- In `space_filling.py`, splittable stubs are reabsorbed into the backbone by
  later growth, which pumps their degree back up.
- Nothing specifically populates degree 3.

## Goal & success criterion

A post-growth decoration phase with one unbounded intensity knob
(`noise_frac`) that the user sweeps, evaluating the resulting DGDC degree
distribution themselves. Sweeping must be cheap: the backbone is grown once and
re-decorated at many noise levels without re-running growth. Two shape knobs
control the degree-1/2/3 mix (`snap_tol`, `branch_prob`).

`noise_frac = 0.0` (default) leaves `run_model` output byte-identical to
today.

## API

New standalone function in `space_filling.py`:

```python
add_noise_streets(segments, noise_frac=0.5, seed=0,
                  min_length=MIN_LENGTH,
                  noise_len_max=None,    # default 3*min_length
                  snap_tol=None,         # default min_length
                  branch_prob=0.3,
                  min_angle=np.pi/4,
                  max_tries=50, verbose=False) -> list of segments
```

- Input: any segment list `[((x1,y1),(x2,y2)), ...]` (e.g. the return value of
  `run_model`). Output: a new segment list containing the (re-split) input
  streets plus the added noise streets.
- Adds `round(noise_frac * n_input_streets)` noise streets. `noise_frac` is
  unbounded above (0.5, 1.0, 2.0, … all valid); negative values raise
  `ValueError`.
- `snap_tol` shifts the degree-1 ↔ degree-2 mix (dangling cul-de-sac vs T-in);
  `snap_tol=0` disables snapping entirely.
- `branch_prob` pushes mass toward degree 2–3 by anchoring noise on noise
  (small trees). Must be in [0, 1] (else `ValueError`).
- Own RNG: `np.random.default_rng(seed)`.

Convenience wiring on `run_model` (forwarded to `add_noise_streets` after the
growth loop, before the return):

| param               | default | forwarded as  |
|---------------------|---------|---------------|
| `noise_frac`        | `0.0`   | `noise_frac`  |
| `noise_len_max`     | `None`  | `noise_len_max` |
| `noise_snap_tol`    | `None`  | `snap_tol`    |
| `noise_branch_prob` | `0.3`   | `branch_prob` |

Inside `run_model` the decoration RNG seed is derived from the run seed as
`seed ^ 0xBADC0FFE` (distinct from the stub constant `0xDEADBEEF`), so the main
and stub streams are untouched and `noise_frac=0.0` is byte-identical to the
current output. `min_length` and `min_angle` are passed through from the run's
own values. This also makes `run_ensemble(seeds, noise_frac=...)` work
unchanged.

## Mechanism (per noise street)

1. **Parent street.** With probability `branch_prob`, pick a previously added
   noise street (if any exist; otherwise fall back to the general pool);
   otherwise pick from all current streets. Selection is length-weighted within
   the chosen pool.
2. **Anchor.** Uniform position along the parent at `t ~ U(0.1, 0.9)`. The
   parent segment is split at the anchor (mid-block T). The two halves are
   collinear, so DGDC merges them back into one stroke: the noise street
   touches exactly one street at its anchor.
3. **Heading & length.** `theta ~ U(0, 2*pi)`, rejected if the angle to the
   host segment's axis is below `min_angle` (in either direction along the
   axis). With the default `min_angle = pi/4 > a_threshold = 20 deg`, DGDC can
   never merge a noise street into its host. Length
   `~ U(min_length, noise_len_max)`.
4. **Free end.** If the free end lands within `snap_tol` of an existing
   street (via `_nearest_segment`, excluding segments incident to the anchor),
   snap to the foot of the perpendicular and split that street at the contact
   point (T-in, dual degree 2). Otherwise the free end dangles (cul-de-sac,
   dual degree 1).
5. **Validity.** The free end must lie inside the unit square
   (`0 <= x, y <= 1`) and the noise segment must not properly cross any
   existing segment (`_segment_blocks`; the host halves share the anchor
   endpoint and are non-blocking by the existing endpoint-sharing rule). On
   failure, redraw everything (parent, anchor, heading, length) up to
   `max_tries` times; if all tries fail the street is skipped (counted, and
   reported when `verbose=True`).

Resulting dual degrees: dangle = 1, snap = 2, and each child branch anchored
on a noise street adds +1 to that parent — so the degree-1/2/3 population is
directly controllable via `snap_tol` and `branch_prob`.

## Implementation notes

- Reuse existing helpers: `_segment_blocks`, `_nearest_segment`, `nkey`. No
  new dependencies; the growth loop and all `p_stub` machinery are untouched.
- Internal state: the same capacity-doubling numpy segment buffer pattern used
  by `run_model`, plus a boolean `is_noise` mask for branch-parent selection.
  Splitting a parent or snapped street uses the swap-remove + append pattern;
  the `is_noise` flag is inherited by both halves.
- Performance: O(n) blocking checks per placement attempt are fine at
  unit-square sizes (thousands of segments).

## Scope

In scope:

- `space_filling_model/space_filling.py`: `add_noise_streets` plus the four
  forwarded `noise_*` parameters on `run_model`.
- A focused test file `space_filling_model/tests/test_noise_streets.py`.

Out of scope:

- `real_space_filling.py` (per user decision: unit-square model only).
- `to_geodataframe`, exports, `run_ensemble` internals (it forwards kwargs
  already).
- Any change to in-loop stub behaviour.

## Verification

1. `add_noise_streets(segs, noise_frac=0.0)` returns the input unchanged, and
   `run_model(noise_frac=0.0)` output is identical to the pre-change output
   for a fixed seed.
2. The number of added streets equals `round(noise_frac * n_input)` minus the
   (normally zero) skipped count; for a converged backbone and moderate
   `noise_frac` the deficit is 0.
3. Planarity: no two output segments properly cross.
4. All noise endpoints lie inside the unit square.
5. Same `seed` → identical output (determinism).
6. Every noise street's heading differs from its host axis by at least
   `min_angle`.
7. `snap_tol=0` → every noise street has a dangling tip (primal degree-1 free
   end).
8. `branch_prob=1.0` → from the second noise street on, parents are noise
   streets.
9. `noise_frac < 0` and `branch_prob` outside [0, 1] raise `ValueError`.
