# Stochastic stub noise for the unit-square space-filling model

**Date:** 2026-06-05
**Component:** `space_filling_model/space_filling.py`
**Status:** Approved (design), pending implementation

## Problem

`real_space_filling.py` already has a stochastic dead-end "stub" mechanism
(`p_stub`/`stub_len_max`/`snap_tol`) that injects random streets each iteration
to populate the low-degree shoulder of the DGDC dual-graph degree distribution.
The unit-square `space_filling.py` model was explicitly left out of scope by that
original spec and has no noise knob. This spec ports the mechanism into
`space_filling.run_model`.

A key difference from the real-area version: here the random streets should be
**first-class participants in the generative process** by default — future
iterations can split them and connect to them — rather than inert dead-ends.
This is exposed as a `stub_splittable` knob so both behaviors are available.

## Goal & success criterion

A tunable, opt-in knob (`p_stub`) on `run_model` that injects random stub
streets. `p_stub=0.0` (default) reproduces the current output byte-for-byte.
Sweeping `p_stub` upward visibly adds short random streets to the network, and
(with `stub_splittable=True`) those streets seed further growth.

## Mechanism (mirrors `real_space_filling._run`)

Each iteration, after the chosen segment is split at its midpoint and the normal
nearest-visible connection is made, draw `u ~ stub_rng.random()`:

- **`u < p_stub`** → *additionally* emit one stub (the normal connection still
  runs; stubs are additive, not branch-or-replace).

A stub is generated as follows:

1. **Separate RNG.** `stub_rng = np.random.default_rng(seed ^ 0xDEADBEEF)`,
   created only when `p_stub > 0`. The main RNG stream (segment selection and
   normal connections) is untouched, so `p_stub=0` is byte-identical to today.
2. **Anchor.** With probability 0.5 the anchor is the fresh split midpoint;
   otherwise a random existing node with primal degree ≥ 3 (tracked via cheap
   node-degree bookkeeping). Falls back to the midpoint if none exists yet.
3. **Direction & length.** Random heading; length `~ U(min_length,
   stub_len_max)`, `stub_len_max` default `3 * min_length`.
4. **Domain.** The stub free end must lie inside the unit square
   (`0 ≤ x,y ≤ 1`).
5. **Planar / valid.** Reuse `_segment_blocks` (no crossing existing streets);
   in model 2 reject headings whose nearest-incident angle falls outside
   `[min_angle, max_angle]` via `_min_incident_angle`.
6. **Optional T-in.** If the free end lands within `snap_tol` of an existing
   street, snap to the foot of the perpendicular and split that street at the
   contact node (the stub gains a real degree-2 junction). Otherwise the free
   end dangles (degree-1 cul-de-sac).
7. **Feedback (`stub_splittable`).** When `True` (default), stub segments are
   `splittable=True` and their tips are valid connection targets, so the model
   builds on them. When `False`, stubs are `splittable=False` and dangling tips
   are excluded as connection targets — the inert-dead-end behavior of
   `real_space_filling.py`.

## New parameters on `run_model`

| param             | default          | meaning                                                                        |
|-------------------|------------------|--------------------------------------------------------------------------------|
| `p_stub`          | `0.0`            | per-iteration probability of emitting a stub (0 = current behavior)            |
| `stub_len_max`    | `3 * min_length` | upper bound of stub length (unit-square units)                                 |
| `snap_tol`        | `min_length`     | distance under which a stub free end T's into a nearby street                  |
| `stub_splittable` | `True`           | whether stubs participate in further growth (splittable + connectable targets) |
| `stub_at_midpoint_prob` | `0.5`      | probability a stub anchors at the fresh midpoint vs a random degree-≥3 junction |

## Supporting changes to the inner loop

Needed to support stubs; none change `p_stub=0` output:

- Add a `splittable` array (capacity-doubled alongside `seg_arr`, swap-removed on
  split). Eligibility becomes `splittable[:n] & (lengths > min_length)`.
- Add lightweight node-degree bookkeeping (`node_deg` array + `_bump` helper +
  `deg3_idx` list) to support degree-≥3 anchor selection.
- Add a local `_nearest_segment` helper (copied into `space_filling.py` to avoid
  a circular import with `real_space_filling`, which imports from this module).
- When `stub_splittable=False`, exclude dangling stub tips (a `stub_node_set`)
  from the candidate keep-mask.

## Scope

In scope:

- `space_filling_model/space_filling.py` only.
- The four new parameters and the additive-stub logic.
- A focused stub test under `space_filling_model/tests/`.

Out of scope:

- `run_ensemble`, `to_geodataframe`, exports, and the per-iteration KD-tree
  structure stay as-is.
- `real_space_filling.py` is untouched.

## Verification

1. `p_stub=0.0` returns a segment list identical to the pre-change `run_model`
   for a fixed seed.
2. `p_stub>0` adds segments (more than the `p_stub=0` run) and includes short
   degree-1 / degree-2 stub streets.
3. With `stub_splittable=True`, stub segments can be split / connected to in
   later iterations; with `False`, dangling tips are never connection targets.
4. Model 2 stubs respect `[min_angle, max_angle]` and never cross existing
   streets.
5. `p_stub` and `stub_at_midpoint_prob` outside `[0, 1]` raise `ValueError`;
   the boundaries `0.0` and `1.0` are accepted.
