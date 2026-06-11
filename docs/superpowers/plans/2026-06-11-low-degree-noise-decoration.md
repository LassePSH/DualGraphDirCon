# Low-Degree Noise Decoration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a post-growth decoration phase `add_noise_streets` to `space_filling_model/space_filling.py` (plus `noise_*` convenience params on `run_model`) so the DGDC dual-graph degree distribution can be tuned flat up to degree 4 by sweeping one unbounded knob, `noise_frac`.

**Architecture:** A standalone function takes a finished segment list and adds `round(noise_frac * n)` short "noise" streets, each anchored mid-block on an existing street (host split at the anchor → the noise street touches exactly one stroke there). Free ends dangle (dual degree 1) or T into a nearby street within `snap_tol` (degree 2); with probability `branch_prob` a noise street anchors on a previously added noise street, pushing parents to degree 2–3. The growth loop and all `p_stub` machinery are untouched; `noise_frac=0.0` is byte-identical to current output.

**Tech Stack:** numpy, existing numba-jitted helpers (`_segment_blocks`), existing `_nearest_segment`, pytest.

**Spec:** `docs/superpowers/specs/2026-06-11-low-degree-noise-decoration-design.md`

---

## File structure

- **Modify** `space_filling_model/space_filling.py`:
  - new public function `add_noise_streets` (insert after `_select_segment`, before the "Parallel ensemble" section);
  - `run_model`: four new keyword params (`noise_frac=0.0`, `noise_len_max=None`, `noise_snap_tol=None`, `noise_branch_prob=0.3`), early validation, and a decoration call at the return site;
  - module + `run_model` docstring updates.
- **Create** `space_filling_model/tests/test_noise_streets.py` (follows the conventions of `tests/test_unit_square_stub_noise.py`: sys.path shim, `_seg_set` helper, small `max_iter` runs).

All commands below run from `/home/lpsha/s154446/fractality/space_filling_model`.

### Geometry notes the implementer must know

- `_segment_blocks(px, py, qx, qy, seg_arr)` treats segments sharing an endpoint (within `EPS=1e-9`) with p or q as non-blocking, and tests *proper* crossing via signed-area products. A point that lies **on** a segment's interior (like our mid-block anchor on its host, or a snapped foot on its target) makes one cross product ≈ ±1e-17, so the host/target row would randomly appear to "block" its own noise street. **Therefore the host row (and a snapped target row) must be hidden during checks.** We hide a row O(1) by temporarily overwriting it with a far-away degenerate sentinel `((-10,-10),(-10,-10))` and restoring it afterwards (a zero-length segment never blocks and is never nearest).
- `_nearest_segment(fx, fy, seg_arr, max_dist)` clips the projection parameter to `[0,1]`, so the returned foot frequently **equals an endpoint** of the matched segment. Splitting there would create a zero-length half; instead we connect to that endpoint without splitting (`endpoint_connect`).
- Snapping moves the free end, so the heading drawn from `theta` is **not** the final direction. The `min_angle` check against the host axis must be re-applied to the final geometry, otherwise a snapped end can land on a node of the host itself, creating a collinear overlap with the host.
- `nkey` rounds to 8 decimals, so any street shorter than ~1e-8 collapses to a single node key. We guard all created geometry with `GUARD = 1e-6`.

---

### Task 1: `add_noise_streets` core function

**Files:**
- Create: `space_filling_model/tests/test_noise_streets.py`
- Modify: `space_filling_model/space_filling.py` (insert function after `_select_segment`, around line 514)

- [ ] **Step 1: Write the failing tests**

Create `space_filling_model/tests/test_noise_streets.py` with exactly:

```python
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
```

- [ ] **Step 2: Run the tests to verify they fail on the import**

Run: `python -m pytest tests/test_noise_streets.py -v`
Expected: collection error — `ImportError: cannot import name 'add_noise_streets' from 'space_filling'`.

- [ ] **Step 3: Implement `add_noise_streets`**

In `space_filling_model/space_filling.py`, insert the following after `_select_segment` (i.e. between the end of `_select_segment` and the `# Parallel ensemble` section banner):

```python
# --------------------------------------------------------------------------- #
# Post-growth low-degree noise decoration                                     #
# --------------------------------------------------------------------------- #

def add_noise_streets(segments, noise_frac=0.5, seed=0,
                      min_length=MIN_LENGTH, noise_len_max=None,
                      snap_tol=None, branch_prob=0.3,
                      min_angle=np.pi / 4, max_tries=50, verbose=False):
    """Decorate a finished network with random low-degree "noise" streets.

    Post-growth counterpart to ``p_stub``: takes the segment list returned by
    `run_model` and adds ``round(noise_frac * len(segments))`` short streets,
    each anchored mid-block on an existing street (the host is split at the
    anchor, so DGDC sees the noise street touch exactly one stroke there).
    Because growth never builds on them, noise streets keep low dual-graph
    degrees: a dangling free end gives degree 1, a snapped (T-in) end gives
    degree 2, and each child anchored on a noise street (``branch_prob``)
    adds +1 to its parent — so the degree-1/2/3 population is directly
    tunable. Sweep ``noise_frac`` on one saved backbone without re-running
    growth.

    Parameters
    ----------
    segments    : list of ((x1,y1),(x2,y2)) street segments (unit square).
    noise_frac  : number of noise streets to add as a fraction of
                  ``len(segments)``; unbounded above, 0.0 is a no-op.
    seed        : RNG seed of the decoration stream.
    min_length  : lower bound of noise street length.
    noise_len_max : upper bound of noise street length (default 3*min_length).
    snap_tol    : free ends within this distance of a street T into it
                  (splitting it, dual degree 2); 0 disables snapping so all
                  noise streets dangle (degree 1). Default min_length.
    branch_prob : probability of anchoring on a previously added noise street
                  (small trees -> parents reach degree 2-3). Falls back to the
                  general pool while no noise street exists yet.
    min_angle   : minimum angle (radians) between a noise street and its
                  host's axis; keep above the DGDC a_threshold so noise is
                  never merged into the host stroke.
    max_tries   : placement attempts per noise street before it is skipped.
    verbose     : print a placement summary.

    Returns a new segment list: the (re-split) input plus the noise streets.
    The result is planar (noise never crosses existing streets) and confined
    to the unit square.

    Raises
    ------
    ValueError
        If ``noise_frac`` is negative or ``branch_prob`` is outside [0, 1].
    """
    if noise_frac < 0.0:
        raise ValueError(f"noise_frac must be >= 0, got {noise_frac!r}")
    if not 0.0 <= branch_prob <= 1.0:
        raise ValueError(f"branch_prob must be in [0, 1], got {branch_prob!r}")
    if noise_len_max is None:
        noise_len_max = 3.0 * min_length
    if snap_tol is None:
        snap_tol = min_length

    target = int(round(noise_frac * len(segments)))
    if target == 0:
        return list(segments)

    rng = np.random.default_rng(seed)

    n = len(segments)
    capacity = max(n + 3 * target, 64)
    seg_arr = np.empty((capacity, 2, 2), dtype=np.float64)
    is_noise = np.zeros(capacity, dtype=bool)
    for i, (a, b) in enumerate(segments):
        seg_arr[i, 0] = a
        seg_arr[i, 1] = b

    def grow(extra):
        nonlocal seg_arr, is_noise
        if n + extra <= seg_arr.shape[0]:
            return
        new_cap = max(seg_arr.shape[0] * 2, n + extra)
        new = np.empty((new_cap, 2, 2), dtype=np.float64)
        new[:n] = seg_arr[:n]
        seg_arr = new
        new_m = np.zeros(new_cap, dtype=bool)
        new_m[:n] = is_noise[:n]
        is_noise = new_m

    def split_row(j, px, py):
        # Overwrite row j with (a, p) and append (p, b): all other row indices
        # stay stable. Both halves inherit the row's noise flag.
        nonlocal n
        bx_, by_ = seg_arr[j, 1]
        grow(1)
        seg_arr[j, 1] = (px, py)
        seg_arr[n, 0] = (px, py)
        seg_arr[n, 1] = (bx_, by_)
        is_noise[n] = is_noise[j]
        n += 1

    # Hiding a row behind this far-away degenerate sentinel excludes it from
    # _segment_blocks / _nearest_segment in O(1): a point lying *on* a row
    # (anchor on host, snapped foot on target) makes the proper-crossing test
    # numerically unstable, so those rows must not take part in the checks.
    _far = ((-10.0, -10.0), (-10.0, -10.0))
    GUARD = 1e-6  # degeneracy guard, safely above nkey's 1e-8 rounding

    placed = 0
    skipped = 0
    for _ in range(target):
        success = False
        for _try in range(max_tries):
            # ---- parent street (length-weighted; optionally noise-only) ----
            use_noise_pool = rng.random() < branch_prob
            if use_noise_pool and is_noise[:n].any():
                pool = np.flatnonzero(is_noise[:n])
            else:
                pool = np.arange(n)
            sa = seg_arr[:n]
            lens = np.hypot(sa[pool, 1, 0] - sa[pool, 0, 0],
                            sa[pool, 1, 1] - sa[pool, 0, 1])
            cdf = np.cumsum(lens)
            h = int(pool[np.searchsorted(cdf, rng.random() * cdf[-1])])
            hax, hay = seg_arr[h, 0]
            hbx, hby = seg_arr[h, 1]
            if np.hypot(hbx - hax, hby - hay) < 10.0 * GUARD:
                continue  # host too short to split cleanly
            host_phi = np.arctan2(hby - hay, hbx - hax)

            # ---- anchor, heading, length ----
            t = 0.1 + 0.8 * rng.random()
            ax0 = hax + t * (hbx - hax)
            ay0 = hay + t * (hby - hay)
            theta = rng.random() * 2.0 * np.pi
            d = (theta - host_phi) % np.pi
            if min(d, np.pi - d) < min_angle:
                continue  # too close to the host axis (DGDC would merge)
            length = min_length + rng.random() * (noise_len_max - min_length)
            fx = ax0 + length * np.cos(theta)
            fy = ay0 + length * np.sin(theta)

            # ---- snap & validity (host, and snapped row, hidden) ----
            saved_h = seg_arr[h].copy()
            seg_arr[h] = _far
            snap_j = -1
            endpoint_connect = False
            saved_j = None
            if snap_tol > 0.0:
                snap_j, fx, fy = _nearest_segment(fx, fy, seg_arr[:n],
                                                  snap_tol)
            if snap_j >= 0:
                ux, uy = seg_arr[snap_j, 0]
                wx, wy = seg_arr[snap_j, 1]
                if np.hypot(fx - ux, fy - uy) < GUARD:
                    fx, fy = float(ux), float(uy)   # join the existing node
                    endpoint_connect = True
                elif np.hypot(fx - wx, fy - wy) < GUARD:
                    fx, fy = float(wx), float(wy)
                    endpoint_connect = True
                else:
                    saved_j = seg_arr[snap_j].copy()
                    seg_arr[snap_j] = _far
            # Snapping moved the free end, so re-check the host-axis angle on
            # the final geometry (otherwise a snapped end can land collinear
            # with — even on — the host).
            df = (np.arctan2(fy - ay0, fx - ax0) - host_phi) % np.pi
            ok = (0.0 <= fx <= 1.0 and 0.0 <= fy <= 1.0
                  and np.hypot(fx - ax0, fy - ay0) >= GUARD
                  and min(df, np.pi - df) >= min_angle
                  and not _segment_blocks(ax0, ay0, fx, fy, seg_arr[:n]))
            seg_arr[h] = saved_h
            if saved_j is not None:
                seg_arr[snap_j] = saved_j
            if not ok:
                continue

            # ---- commit: split host (and snapped street), append noise ----
            split_row(h, ax0, ay0)
            if snap_j >= 0 and not endpoint_connect:
                split_row(snap_j, fx, fy)
            grow(1)
            seg_arr[n, 0] = (ax0, ay0)
            seg_arr[n, 1] = (fx, fy)
            is_noise[n] = True
            n += 1
            success = True
            break
        if success:
            placed += 1
        else:
            skipped += 1

    if verbose:
        print(f"  noise: placed {placed}/{target} noise streets "
              f"({skipped} skipped)")

    return [(tuple(seg_arr[i, 0]), tuple(seg_arr[i, 1])) for i in range(n)]
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python -m pytest tests/test_noise_streets.py -v`
Expected: all tests PASS. If `test_added_segment_count_bounds` fails because placements were skipped, do not loosen the assertion — investigate why placement failed (it indicates a bug in the retry/validity logic at this density).

- [ ] **Step 5: Run the existing suites to verify nothing regressed**

Run: `python -m pytest tests/ -v`
Expected: all tests in `test_segment_selection.py`, `test_stub_noise.py`, `test_unit_square_stub_noise.py`, and `test_noise_streets.py` PASS.

- [ ] **Step 6: Commit**

```bash
git add tests/test_noise_streets.py space_filling.py
git commit -m "feat: add post-growth low-degree noise decoration (add_noise_streets)"
```

---

### Task 2: `run_model` wiring (`noise_*` params)

**Files:**
- Modify: `space_filling_model/space_filling.py` (`run_model` signature ~line 196, validation block ~line 249, return site ~line 479)
- Modify: `space_filling_model/tests/test_noise_streets.py` (append tests)

- [ ] **Step 1: Append the failing wiring tests**

Append to `space_filling_model/tests/test_noise_streets.py`:

```python
# --------------------------------------------------------------------------- #
# run_model wiring                                                            #
# --------------------------------------------------------------------------- #

def test_run_model_noise_frac_zero_is_byte_identical():
    base = run_model(model=1, seed=42, max_iter=300, verbose=False)
    zero = run_model(model=1, seed=42, max_iter=300, verbose=False,
                     noise_frac=0.0)
    assert _seg_set(base) == _seg_set(zero)


def test_run_model_noise_frac_adds_segments_deterministically():
    base = run_model(model=1, seed=7, max_iter=300, verbose=False)
    a = run_model(model=1, seed=7, max_iter=300, verbose=False, noise_frac=0.5)
    b = run_model(model=1, seed=7, max_iter=300, verbose=False, noise_frac=0.5)
    assert _seg_set(a) == _seg_set(b)
    assert len(a) > len(base)


def test_run_model_noise_matches_manual_decoration():
    """run_model(noise_frac=f) == run_model() + add_noise_streets(seed^const)."""
    base = run_model(model=1, seed=7, max_iter=300, verbose=False)
    manual = add_noise_streets(base, noise_frac=0.5, seed=7 ^ 0xBADC0FFE)
    wired = run_model(model=1, seed=7, max_iter=300, verbose=False,
                      noise_frac=0.5)
    assert _seg_set(wired) == _seg_set(manual)


def test_run_model_negative_noise_frac_raises_early():
    with pytest.raises(ValueError):
        run_model(model=1, seed=1, max_iter=10, verbose=False, noise_frac=-0.5)


@pytest.mark.parametrize("bad", [-0.1, 1.5])
def test_run_model_noise_branch_prob_out_of_range_raises(bad):
    with pytest.raises(ValueError):
        run_model(model=1, seed=1, max_iter=10, verbose=False,
                  noise_frac=0.2, noise_branch_prob=bad)
```

- [ ] **Step 2: Run them to verify they fail**

Run: `python -m pytest tests/test_noise_streets.py -k run_model -v`
Expected: FAIL with `TypeError: run_model() got an unexpected keyword argument 'noise_frac'` (the byte-identical test passes trivially only once the kwarg exists; the others must fail now).

- [ ] **Step 3: Wire the parameters into `run_model`**

Three edits in `space_filling_model/space_filling.py`:

(a) Signature — change

```python
def run_model(model=1, min_length=MIN_LENGTH, max_iter=200_000, seed=42,
              min_angle=np.pi / 4, max_angle=None, area_coeff=0.05,
              area_scale=1.0, select='length', select_power=1.0,
              p_stub=0.0, stub_len_max=None, snap_tol=None,
              stub_splittable=True, stub_at_midpoint_prob=0.5, verbose=True):
```

to

```python
def run_model(model=1, min_length=MIN_LENGTH, max_iter=200_000, seed=42,
              min_angle=np.pi / 4, max_angle=None, area_coeff=0.05,
              area_scale=1.0, select='length', select_power=1.0,
              p_stub=0.0, stub_len_max=None, snap_tol=None,
              stub_splittable=True, stub_at_midpoint_prob=0.5,
              noise_frac=0.0, noise_len_max=None, noise_snap_tol=None,
              noise_branch_prob=0.3, verbose=True):
```

(b) Validation — directly after the existing `stub_at_midpoint_prob` ValueError block, add:

```python
    if noise_frac < 0.0:
        raise ValueError(f"noise_frac must be >= 0, got {noise_frac!r}")
    if not 0.0 <= noise_branch_prob <= 1.0:
        raise ValueError(
            f"noise_branch_prob must be in [0, 1], got {noise_branch_prob!r}")
```

(c) Return site — change the last lines of `run_model` from

```python
    return [(tuple(seg_arr[i, 0]), tuple(seg_arr[i, 1])) for i in range(n)]
```

to

```python
    segments = [(tuple(seg_arr[i, 0]), tuple(seg_arr[i, 1])) for i in range(n)]
    if noise_frac > 0.0:
        # Post-growth decoration on its own RNG stream: noise_frac=0.0 is
        # byte-identical to a run without it.
        segments = add_noise_streets(
            segments, noise_frac=noise_frac, seed=seed ^ 0xBADC0FFE,
            min_length=min_length, noise_len_max=noise_len_max,
            snap_tol=noise_snap_tol, branch_prob=noise_branch_prob,
            min_angle=min_angle, verbose=verbose)
    return segments
```

(`add_noise_streets` is defined later in the module; that's fine — the name is resolved at call time.)

- [ ] **Step 4: Run the wiring tests to verify they pass**

Run: `python -m pytest tests/test_noise_streets.py -k run_model -v`
Expected: PASS.

- [ ] **Step 5: Run all suites**

Run: `python -m pytest tests/ -v`
Expected: all PASS (in particular `test_pstub_zero_is_byte_identical_to_default` still passes — the growth loop is untouched).

- [ ] **Step 6: Commit**

```bash
git add tests/test_noise_streets.py space_filling.py
git commit -m "feat: wire noise_frac decoration params into run_model"
```

---

### Task 3: Docstrings

**Files:**
- Modify: `space_filling_model/space_filling.py` (module docstring ~line 28, `run_model` docstring Parameters/Raises sections)

- [ ] **Step 1: Module docstring**

In the module docstring, directly after the `p_stub` paragraph (the one ending `...or inert dead-ends when False.`), insert:

```
Use ``add_noise_streets`` (or ``run_model(noise_frac=...)``) to decorate a
finished network with random low-degree noise streets *after* growth:
mid-block anchored cul-de-sacs (dual degree 1), T-ins (degree 2) and small
branching trees (degree 2-3). The amount is unbounded and the backbone is
unchanged, so one saved run can be re-decorated at many ``noise_frac`` levels
and re-evaluated under DGDC without re-running growth.
```

- [ ] **Step 2: `run_model` docstring**

In the `run_model` docstring Parameters list, after the `stub_at_midpoint_prob` entry and before `verbose`, insert:

```
    noise_frac : post-growth decoration intensity: after the growth loop,
                 add round(noise_frac * n_segments) random low-degree noise
                 streets via `add_noise_streets` (0.0 = off; unbounded above).
                 Runs on RNG seed `seed ^ 0xBADC0FFE`, so 0.0 is
                 byte-identical to a run without it.
    noise_len_max : upper bound of noise street length (default 3*min_length).
    noise_snap_tol : distance under which a noise free end T's into a nearby
                 street (default min_length); 0 yields only dangling
                 degree-1 cul-de-sacs.
    noise_branch_prob : probability a noise street anchors on a previously
                 added noise street (default 0.3), pushing parents to dual
                 degree 2-3.
```

And extend the Raises section sentence to cover the new parameters, e.g.:

```
    ValueError
        If `p_stub` or `stub_at_midpoint_prob` is outside [0, 1], if
        `noise_frac` is negative, if `noise_branch_prob` is outside [0, 1],
        or `select` is not a recognised strategy.
```

- [ ] **Step 3: Run all suites once more**

Run: `python -m pytest tests/ -v`
Expected: all PASS (docstring-only change).

- [ ] **Step 4: Commit**

```bash
git add space_filling.py
git commit -m "docs: document add_noise_streets and run_model noise_* params"
```

---

## Verification checklist (maps to spec)

| Spec verification item | Test |
|---|---|
| 1. no-op / byte-identical | `test_noise_frac_zero_is_a_noop`, `test_run_model_noise_frac_zero_is_byte_identical` |
| 2. added count | `test_added_segment_count_bounds`, `test_snap_tol_zero_yields_only_dangling_tips` |
| 3. planarity | `test_planarity` |
| 4. unit square | `test_endpoints_stay_inside_unit_square` |
| 5. determinism | `test_determinism`, `test_run_model_noise_frac_adds_segments_deterministically` |
| 6. min_angle vs host | `test_noise_heading_respects_min_angle` |
| 7. snap_tol=0 → dangles | `test_snap_tol_zero_yields_only_dangling_tips`, `test_snap_creates_t_junctions` |
| 8. branch_prob=1 | `test_branch_prob_one_anchors_on_noise` |
| 9. ValueErrors | `test_negative_noise_frac_raises`, `test_branch_prob_out_of_range_raises`, `test_run_model_negative_noise_frac_raises_early`, `test_run_model_noise_branch_prob_out_of_range_raises` |
| seed derivation `seed ^ 0xBADC0FFE` | `test_run_model_noise_matches_manual_decoration` |
