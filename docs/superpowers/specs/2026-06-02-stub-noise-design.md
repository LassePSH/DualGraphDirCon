# Stochastic dead-end stubs for the real-area space-filling model

**Date:** 2026-06-02
**Component:** `space_filling_model/real_space_filling.py`
**Status:** Approved (design), pending implementation plan

## Problem

The DGDC dual-graph degree distribution `p(k)` of real street networks is
roughly flat/noisy up to a crossover around degree φ=4 and then follows a
power law (scale-free tail). The real-area space-filling generator currently
reproduces the power-law tail but **starves the low-degree shoulder** (k ≈ 1–4),
especially degree-2 strokes.

Cause: `_run` does exactly one thing per iteration — pick a segment ∝ length,
split at its midpoint, connect the midpoint to the **nearest visible** valid
target. This greedy rule (a) makes nearly every new node a T-junction and (b)
funnels connections onto a few well-placed early strokes (the tail). Nothing
ever produces a short, peripheral street that touches only 1–2 other strokes
and is then left alone.

## Goal & success criterion

Populate the low-degree shoulder (k ≈ 1–4) of the DGDC dual-graph degree
distribution so that `p(k)` is roughly flat/noisy up to φ=4 and retains its
power-law tail above it.

**Success:** a tunable knob (`p_stub`) that, swept upward from 0, visibly raises
low-`k` mass while leaving the high-`k` tail slope essentially unchanged, as
verified by re-running the `plot_real_space_filling.py` DGDC pipeline.

## Mechanism: branch-or-stub in `_run`

Each iteration, after the chosen segment is split at its midpoint, draw
`u ~ rng.random()`:

- **`u < p_stub`** → emit a **dead-end stub** instead of the normal
  nearest-visible connection.
- **else** → existing behavior (connect midpoint to nearest valid target).

A stub is generated as follows:

1. **Anchor (random).** With probability 0.5 the anchor is the fresh split
   midpoint; otherwise it is a random existing node with primal degree ≥ 3
   (tracked via cheap node-degree bookkeeping). Falls back to the midpoint if no
   degree-≥3 node exists yet. Anchoring at real intersections is what lets stubs
   reach degree 2+ rather than only degree 1.
2. **Direction & length.** Random heading; length `~ U(min_length,
   stub_len_max)` in metres (`stub_len_max` default `3 * min_length`). In
   model 2, reject headings that violate `min_angle` at the anchor (reuse
   `_min_incident_angle`).
3. **Blocking.** Reuse `_segment_blocks` — a stub may not cross an existing
   street.
4. **Optional T-in.** If the free end lands within `snap_tol` of an existing
   segment, snap the endpoint to the foot of the perpendicular and **split that
   segment** at the contact node, so DGDC sees a real junction (the stub stroke
   gains degree 2). Otherwise the free end dangles (degree-1 cul-de-sac).
5. **Non-splittable.** Stubs are added with `splittable=False` so they stay
   short and keep their low degree (they will not later be subdivided into
   hubs).

This gives a natural spread across k=1–4: pure dangling cul-de-sacs (k≈1),
intersection-anchored stubs and T-ins (k≈2+).

## New parameters

Threaded through `_run` → `run_circular_model`.

| param          | default          | meaning                                                                 |
|----------------|------------------|-------------------------------------------------------------------------|
| `p_stub`       | `0.0`            | per-iteration probability of emitting a stub (0 = current behavior)     |
| `stub_len_max` | `3 * min_length` | upper bound of stub length (metres)                                     |
| `snap_tol`     | `min_length`     | distance under which a stub free end T's into a nearby street           |

Defaulting `p_stub=0.0` keeps every existing call byte-for-byte unchanged; the
noise is strictly opt-in. Recommended starting sweep value: `p_stub = 0.15`.

## Scope

In scope:

- `space_filling_model/real_space_filling.py` only.
- The three new parameters and the branch-or-stub logic.

Out of scope (explicitly not doing):

- The unit-square `space_filling.py` model is untouched.
- The other candidate noise types (short connectors, randomized target
  selection, geometry jitter) are not implemented.
- No automatic fitting of `p_stub` to a target curve — the user sweeps it
  manually.

## Verification

Run `plot_real_space_filling.py` (or a small sweep helper) at
`p_stub ∈ {0, 0.1, 0.2, 0.3}`, regenerate the `p(k)` plot, and confirm:

1. The low-`k` shoulder (k ≈ 1–4) rises as `p_stub` increases.
2. The high-`k` tail slope holds roughly constant.
3. The raw-lines plot shows short dangling stub segments.
