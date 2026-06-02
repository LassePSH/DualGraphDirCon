# Stochastic dead-end stubs Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an opt-in `p_stub` knob to the real-area space-filling generator that injects random dead-end stubs, populating the low-degree shoulder (k ≈ 1–4) of the DGDC dual-graph degree distribution without disturbing the power-law tail.

**Architecture:** Modify only `space_filling_model/real_space_filling.py`. Inside the main loop `_run`, after each midpoint split, with probability `p_stub` emit a short stub instead of the usual nearest-visible connection. A stub anchors at the fresh midpoint or (50/50) at a random existing intersection node (primal degree ≥ 3), grows a random short distance, and either dangles (degree-1 cul-de-sac) or T's into a nearby street (degree-2). This requires (a) decoupling "is a street" from "is splittable" so non-splittable stubs are still returned, and (b) cheap monotonic node-degree bookkeeping to sample degree-≥3 anchors.

**Tech Stack:** Python, NumPy, SciPy (`cKDTree`), Shapely, numba (`@njit` hot loops), pytest.

**Reference spec:** `docs/superpowers/specs/2026-06-02-stub-noise-design.md`

---

## File structure

- Modify: `space_filling_model/real_space_filling.py`
  - extend the top-level import from `space_filling` to also bring in `_min_incident_angle` and `_segment_blocks`
  - add module-level helper `_nearest_segment(...)`
  - refactor `_run`: add `is_street` array, node-degree bookkeeping, three new params, and the branch-or-stub block
  - thread the three new params through `run_circular_model`
- Create: `space_filling_model/tests/test_stub_noise.py` (the repo currently has no test suite; this is the first one)

Conventions to follow (already in the file): capacity-doubling NumPy buffers, `nkey()` coordinate keys, numba `@njit` hot loops imported from `space_filling`, metres throughout.

---

## Task 1: Characterization test locking current `_run` behavior

This pins today's geometry so the later refactor is provably behavior-preserving. The golden values were produced on this machine from the current code.

**Files:**
- Create: `space_filling_model/tests/test_stub_noise.py`

- [ ] **Step 1: Write the characterization test**

```python
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
```

- [ ] **Step 2: Run it to confirm it passes against current code**

Run: `cd /home/lpsha/s154446/fractality && python -m pytest space_filling_model/tests/test_stub_noise.py::test_baseline_geometry_unchanged -v`
Expected: PASS (1 passed). If it FAILS, the golden values drifted (e.g. different numba/NumPy build) — regenerate them with the snippet in the spec's Verification section and update the two literals before continuing.

- [ ] **Step 3: Commit**

```bash
cd /home/lpsha/s154446/fractality
git add space_filling_model/tests/test_stub_noise.py
git commit -m "test: characterize current real space-filling geometry"
```

---

## Task 2: Decouple `is_street` from `splittable` + node-degree bookkeeping

No behavior change. After this task the characterization test must still pass. This prepares the two pieces of state the stub logic needs: stubs will be non-splittable yet still returned, and degree-≥3 nodes must be samplable cheaply.

**Files:**
- Modify: `space_filling_model/real_space_filling.py`

- [ ] **Step 1: Extend the `space_filling` import**

Replace (lines ~33-36):

```python
from space_filling import (  # noqa: E402  reuse hot inner loops
    _find_best,
    EPS,
)
```

with:

```python
from space_filling import (  # noqa: E402  reuse hot inner loops
    _find_best,
    _segment_blocks,
    _min_incident_angle,
    EPS,
)
```

- [ ] **Step 2: Allocate the `is_street` buffer alongside `splittable`**

In `_run`, replace the buffer setup (the block beginning `capacity = max((n_init + n_block) * 4, 256)` down through `n = n_init + n_block`):

```python
    capacity = max((n_init + n_block) * 4, 256)
    seg_arr = np.empty((capacity, 2, 2), dtype=np.float64)
    splittable = np.zeros(capacity, dtype=bool)

    for i, (a, b) in enumerate(seg_init):
        seg_arr[i, 0] = a; seg_arr[i, 1] = b
        splittable[i] = True
    for j, (a, b) in enumerate(blocking_segs):
        k = n_init + j
        seg_arr[k, 0] = a; seg_arr[k, 1] = b
        splittable[k] = False
    n = n_init + n_block
```

with (adds the parallel `is_street` array; seeds are streets, blockers are not):

```python
    capacity = max((n_init + n_block) * 4, 256)
    seg_arr = np.empty((capacity, 2, 2), dtype=np.float64)
    splittable = np.zeros(capacity, dtype=bool)
    is_street = np.zeros(capacity, dtype=bool)

    for i, (a, b) in enumerate(seg_init):
        seg_arr[i, 0] = a; seg_arr[i, 1] = b
        splittable[i] = True; is_street[i] = True
    for j, (a, b) in enumerate(blocking_segs):
        k = n_init + j
        seg_arr[k, 0] = a; seg_arr[k, 1] = b
        splittable[k] = False; is_street[k] = False
    n = n_init + n_block
```

- [ ] **Step 3: Add a `node_deg` buffer and grow it inside `add_node`**

Replace the node-table setup + `add_node` (the block from `node_keys = {}` through the end of `add_node`, i.e. its `return n_nodes - 1`):

```python
    node_keys = {}
    node_arr = np.empty((capacity * 2, 2), dtype=np.float64)
    n_nodes = 0

    def add_node(pt):
        nonlocal n_nodes, node_arr
        k = nkey(pt)
        idx = node_keys.get(k)
        if idx is not None:
            return idx
        if n_nodes >= node_arr.shape[0]:
            new = np.empty((node_arr.shape[0] * 2, 2), dtype=np.float64)
            new[:n_nodes] = node_arr[:n_nodes]
            node_arr = new
        node_arr[n_nodes] = pt
        node_keys[k] = n_nodes
        n_nodes += 1
        return n_nodes - 1
```

with (parallel `node_deg`, a monotonic `deg3_idx` list, and a `_bump` helper):

```python
    node_keys = {}
    node_arr = np.empty((capacity * 2, 2), dtype=np.float64)
    node_deg = np.zeros(capacity * 2, dtype=np.int64)
    n_nodes = 0
    deg3_idx = []  # node indices whose primal degree has reached >= 3

    def add_node(pt):
        nonlocal n_nodes, node_arr, node_deg
        k = nkey(pt)
        idx = node_keys.get(k)
        if idx is not None:
            return idx
        if n_nodes >= node_arr.shape[0]:
            new = np.empty((node_arr.shape[0] * 2, 2), dtype=np.float64)
            new[:n_nodes] = node_arr[:n_nodes]
            node_arr = new
            new_deg = np.zeros(node_deg.shape[0] * 2, dtype=np.int64)
            new_deg[:n_nodes] = node_deg[:n_nodes]
            node_deg = new_deg
        node_arr[n_nodes] = pt
        node_keys[k] = n_nodes
        n_nodes += 1
        return n_nodes - 1

    def _bump(idx, amt):
        # Degrees only ever increase in this model, so a node enters deg3_idx
        # exactly once, when it first reaches degree 3.
        old = int(node_deg[idx])
        node_deg[idx] = old + amt
        if old < 3 <= old + amt:
            deg3_idx.append(int(idx))
```

- [ ] **Step 4: Grow `is_street` inside `grow`**

Replace the `grow` helper:

```python
    def grow(extra):
        nonlocal seg_arr, splittable
        if n + extra <= seg_arr.shape[0]:
            return
        new_cap = max(seg_arr.shape[0] * 2, n + extra)
        new = np.empty((new_cap, 2, 2), dtype=np.float64)
        new[:n] = seg_arr[:n]
        seg_arr = new
        new_s = np.zeros(new_cap, dtype=bool)
        new_s[:n] = splittable[:n]
        splittable = new_s
```

with:

```python
    def grow(extra):
        nonlocal seg_arr, splittable, is_street
        if n + extra <= seg_arr.shape[0]:
            return
        new_cap = max(seg_arr.shape[0] * 2, n + extra)
        new = np.empty((new_cap, 2, 2), dtype=np.float64)
        new[:n] = seg_arr[:n]
        seg_arr = new
        new_s = np.zeros(new_cap, dtype=bool)
        new_s[:n] = splittable[:n]
        splittable = new_s
        new_st = np.zeros(new_cap, dtype=bool)
        new_st[:n] = is_street[:n]
        is_street = new_st
```

- [ ] **Step 5: Initialise degrees for seed + blocker endpoints**

Replace the two init loops:

```python
    for a, b in seg_init:
        add_node(a); add_node(b)
    for a, b in blocking_segs:
        add_node(a); add_node(b)
```

with:

```python
    for a, b in seg_init:
        ia = add_node(a); ib = add_node(b); _bump(ia, 1); _bump(ib, 1)
    for a, b in blocking_segs:
        ia = add_node(a); ib = add_node(b); _bump(ia, 1); _bump(ib, 1)
```

- [ ] **Step 6: Maintain degree + `is_street` through the split and connect**

In the main loop, the split currently does a swap-remove then appends two halves. Replace:

```python
        # split: swap-remove the parent, append two halves (both splittable)
        grow(2)
        last_split = splittable[n - 1]
        seg_arr[seg_i] = seg_arr[n - 1]
        splittable[seg_i] = last_split
        n -= 1
        seg_arr[n, 0] = (ax, ay); seg_arr[n, 1] = (midx, midy)
        splittable[n] = True; n += 1
        seg_arr[n, 0] = (midx, midy); seg_arr[n, 1] = (bx, by)
        splittable[n] = True; n += 1

        m_idx = add_node((midx, midy))
        a_idx = node_keys[nkey((ax, ay))]
        b_idx = node_keys[nkey((bx, by))]
```

with (carry `is_street` through the swap-remove; the midpoint gains degree 2 from the two halves):

```python
        # split: swap-remove the parent, append two halves (both splittable)
        grow(2)
        last_split = splittable[n - 1]
        last_street = is_street[n - 1]
        seg_arr[seg_i] = seg_arr[n - 1]
        splittable[seg_i] = last_split
        is_street[seg_i] = last_street
        n -= 1
        seg_arr[n, 0] = (ax, ay); seg_arr[n, 1] = (midx, midy)
        splittable[n] = True; is_street[n] = True; n += 1
        seg_arr[n, 0] = (midx, midy); seg_arr[n, 1] = (bx, by)
        splittable[n] = True; is_street[n] = True; n += 1

        m_idx = add_node((midx, midy))
        a_idx = node_keys[nkey((ax, ay))]
        b_idx = node_keys[nkey((bx, by))]
        _bump(m_idx, 2)  # midpoint now carries the two halves
```

- [ ] **Step 7: Record degree on the normal connection**

Replace the connection-commit block:

```python
        if best is not None:
            grow(1)
            seg_arr[n, 0] = (midx, midy)
            seg_arr[n, 1] = best
            splittable[n] = True
            n += 1
            add_node(best)
```

with:

```python
        if best is not None:
            grow(1)
            seg_arr[n, 0] = (midx, midy)
            seg_arr[n, 1] = best
            splittable[n] = True; is_street[n] = True
            n += 1
            best_idx = add_node(best)
            _bump(m_idx, 1)        # midpoint now a T-junction (degree 3)
            _bump(best_idx, 1)
```

- [ ] **Step 8: Return streets by `is_street` instead of `splittable`**

Replace:

```python
    # return only the splittable (street) segments; boundary stays separate
    return [(tuple(seg_arr[i, 0]), tuple(seg_arr[i, 1]))
            for i in range(n) if splittable[i]]
```

with:

```python
    # return only the street segments; boundary blockers stay separate
    return [(tuple(seg_arr[i, 0]), tuple(seg_arr[i, 1]))
            for i in range(n) if is_street[i]]
```

- [ ] **Step 9: Run the characterization test (behavior must be unchanged)**

Run: `cd /home/lpsha/s154446/fractality && python -m pytest space_filling_model/tests/test_stub_noise.py -v`
Expected: PASS (1 passed). A mismatch means the refactor altered geometry — review Steps 2–8.

- [ ] **Step 10: Commit**

```bash
cd /home/lpsha/s154446/fractality
git add space_filling_model/real_space_filling.py
git commit -m "refactor: decouple is_street from splittable + track node degree"
```

---

## Task 3: `_nearest_segment` helper for T-in snapping

A pure NumPy point-to-segment query: given a free end and the live segment array, return the nearest segment index and the foot of the perpendicular, but only if within `max_dist`. Segments incident to the anchor are excluded so a stub never snaps back onto its own anchor edges.

**Files:**
- Modify: `space_filling_model/real_space_filling.py`
- Modify: `space_filling_model/tests/test_stub_noise.py`

- [ ] **Step 1: Write the failing unit test**

Append to `space_filling_model/tests/test_stub_noise.py`:

```python
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
```

- [ ] **Step 2: Run it to verify it fails**

Run: `cd /home/lpsha/s154446/fractality && python -m pytest space_filling_model/tests/test_stub_noise.py -k nearest_segment -v`
Expected: FAIL with `ImportError: cannot import name '_nearest_segment'`.

- [ ] **Step 3: Implement the helper**

Add this function in `real_space_filling.py` immediately above `_run` (after the `init_spokes` section):

```python
def _nearest_segment(fx, fy, seg_arr, max_dist, exclude_xy=None):
    """Nearest segment to point (fx,fy) within `max_dist`.

    Returns (idx, foot_x, foot_y): the index into `seg_arr` of the closest
    segment and the foot of the perpendicular from the point onto it. If no
    segment lies within `max_dist`, returns (-1, fx, fy). Segments sharing the
    `exclude_xy` endpoint (the stub anchor) are never selected.
    """
    if seg_arr.shape[0] == 0:
        return -1, fx, fy
    ax = seg_arr[:, 0, 0]; ay = seg_arr[:, 0, 1]
    bx = seg_arr[:, 1, 0]; by = seg_arr[:, 1, 1]
    dx = bx - ax; dy = by - ay
    L2 = dx * dx + dy * dy
    L2 = np.where(L2 < EPS, EPS, L2)
    t = ((fx - ax) * dx + (fy - ay) * dy) / L2
    t = np.clip(t, 0.0, 1.0)
    projx = ax + t * dx; projy = ay + t * dy
    dist = np.hypot(projx - fx, projy - fy)
    if exclude_xy is not None:
        ex, ey = exclude_xy
        share = (((np.abs(ax - ex) < EPS) & (np.abs(ay - ey) < EPS)) |
                 ((np.abs(bx - ex) < EPS) & (np.abs(by - ey) < EPS)))
        dist = np.where(share, np.inf, dist)
    i = int(np.argmin(dist))
    if dist[i] <= max_dist:
        return i, float(projx[i]), float(projy[i])
    return -1, fx, fy
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd /home/lpsha/s154446/fractality && python -m pytest space_filling_model/tests/test_stub_noise.py -k nearest_segment -v`
Expected: PASS (3 passed).

- [ ] **Step 5: Commit**

```bash
cd /home/lpsha/s154446/fractality
git add space_filling_model/real_space_filling.py space_filling_model/tests/test_stub_noise.py
git commit -m "feat: add _nearest_segment point-to-segment helper"
```

---

## Task 4: Branch-or-stub logic in `_run`

Add the three parameters and the stub emission path. With probability `p_stub`, replace the normal connect with a stub. Default `p_stub=0.0` keeps existing behavior exactly.

**Files:**
- Modify: `space_filling_model/real_space_filling.py`
- Modify: `space_filling_model/tests/test_stub_noise.py`

- [ ] **Step 1: Write failing behavior tests**

Append to `space_filling_model/tests/test_stub_noise.py`:

```python
def test_p_stub_zero_is_identical():
    """p_stub=0.0 must reproduce the baseline geometry exactly."""
    res = run_circular_model(p_stub=0.0, **BASE_KW)
    assert len(res['segments']) == 3728
    assert _geom_hash(res['segments']) == '913fac0a0326c33b'


def test_stubs_add_low_degree_deadends():
    """Turning on stubs must raise the count of degree-1 dead-end nodes."""
    base = run_circular_model(p_stub=0.0, **BASE_KW)['segments']
    noisy = run_circular_model(p_stub=0.3, **BASE_KW)['segments']

    base_deg1 = sum(1 for v in _degree_counter(base).values() if v == 1)
    noisy_deg1 = sum(1 for v in _degree_counter(noisy).values() if v == 1)

    assert noisy_deg1 > base_deg1


def test_stub_len_max_bounds_new_deadends():
    """Tightening stub_len_max must not create short dead-ends longer than the
    bound. We compare against the baseline so pre-existing long degree-1
    endpoints (boundary spokes) are excluded."""
    stub_len_max = 150.0

    def short_deadend_lengths(segments):
        counter = _degree_counter(segments)
        out = []
        for a, b in segments:
            ka = (round(a[0], 6), round(a[1], 6))
            kb = (round(b[0], 6), round(b[1], 6))
            if counter[ka] == 1 or counter[kb] == 1:
                out.append(float(np.hypot(b[0] - a[0], b[1] - a[1])))
        return out

    base = short_deadend_lengths(run_circular_model(p_stub=0.0, **BASE_KW)['segments'])
    noisy = short_deadend_lengths(
        run_circular_model(p_stub=0.3, stub_len_max=stub_len_max, **BASE_KW)['segments'])

    # New short dead-ends (those at or below the bound) only appear with stubs
    # on, and every one of them must respect the bound.
    base_short = sum(1 for L in base if L <= stub_len_max + 1e-6)
    noisy_short = sum(1 for L in noisy if L <= stub_len_max + 1e-6)
    assert noisy_short > base_short
    for L in noisy:
        if L <= stub_len_max + 1e-6:
            continue  # within bound
        # Any dead-end longer than the bound must already exist at baseline
        # (e.g. boundary spokes), not be a freshly minted stub.
        assert any(abs(L - Lb) < 1e-6 for Lb in base)


def test_reproducible_with_seed():
    a = run_circular_model(p_stub=0.3, **BASE_KW)['segments']
    b = run_circular_model(p_stub=0.3, **BASE_KW)['segments']
    assert _geom_hash(a) == _geom_hash(b)
```

- [ ] **Step 2: Run to verify failure**

Run: `cd /home/lpsha/s154446/fractality && python -m pytest space_filling_model/tests/test_stub_noise.py -k "stub or reproducible" -v`
Expected: FAIL with `TypeError: _run() got an unexpected keyword argument 'p_stub'` (raised via `run_circular_model`).

- [ ] **Step 3: Add the three params to `_run`'s signature**

Replace the `_run` signature:

```python
def _run(seg_init, blocking_segs, area_poly, *,
         model=1, min_length=50.0, max_iter=200_000, seed=42,
         min_angle=np.pi / 4, area_coeff=2500.0, area_scale=None,
         kdtree_rebuild_every=256, progress_every=0, verbose=True):
```

with:

```python
def _run(seg_init, blocking_segs, area_poly, *,
         model=1, min_length=50.0, max_iter=200_000, seed=42,
         min_angle=np.pi / 4, area_coeff=2500.0, area_scale=None,
         p_stub=0.0, stub_len_max=None, snap_tol=None,
         kdtree_rebuild_every=256, progress_every=0, verbose=True):
```

- [ ] **Step 4: Resolve stub defaults near the top of `_run`**

Immediately after `rng = np.random.default_rng(seed)` add:

```python
    if stub_len_max is None:
        stub_len_max = 3.0 * min_length
    if snap_tol is None:
        snap_tol = min_length
```

- [ ] **Step 5: Insert the branch-or-stub block before the candidate search**

In the main loop, the code currently flows from the split/`_bump(m_idx, 2)` directly into `if n_nodes <= 3: continue` and the candidate search. Insert the stub block between them. Replace:

```python
        _bump(m_idx, 2)  # midpoint now carries the two halves

        if n_nodes <= 3:
            continue
```

with:

```python
        _bump(m_idx, 2)  # midpoint now carries the two halves

        # --- branch-or-stub: occasionally emit a dead-end stub instead of the
        # usual nearest-visible connection. Stubs are streets but never split. ---
        if p_stub > 0.0 and rng.random() < p_stub:
            use_mid = (len(deg3_idx) == 0) or (rng.random() < 0.5)
            anchor_idx = m_idx if use_mid else \
                deg3_idx[int(rng.integers(len(deg3_idx)))]
            ax0, ay0 = node_arr[anchor_idx]
            theta = rng.random() * 2.0 * np.pi
            length = min_length + rng.random() * (stub_len_max - min_length)
            fx = ax0 + length * np.cos(theta)
            fy = ay0 + length * np.sin(theta)

            # Optional T-in: snap the free end onto a nearby street (degree-2),
            # excluding edges incident to the anchor.
            snap_i, fx, fy = _nearest_segment(fx, fy, seg_arr[:n], snap_tol,
                                               exclude_xy=(ax0, ay0))

            ok = area_poly.contains(Point(fx, fy))
            if ok and _segment_blocks(ax0, ay0, fx, fy, seg_arr[:n]):
                ok = False
            if ok and model == 2:
                vx, vy = fx - ax0, fy - ay0
                if _min_incident_angle(ax0, ay0, vx, vy, seg_arr[:n]) < min_angle:
                    ok = False

            if ok:
                if snap_i >= 0:
                    # Split the snapped street at the contact point so DGDC sees
                    # a real junction (stub stroke gains degree 2).
                    ux, uy = seg_arr[snap_i, 0]
                    wx, wy = seg_arr[snap_i, 1]
                    s_split = splittable[snap_i]
                    s_street = is_street[snap_i]
                    grow(2)
                    seg_arr[snap_i] = seg_arr[n - 1]
                    splittable[snap_i] = splittable[n - 1]
                    is_street[snap_i] = is_street[n - 1]
                    n -= 1
                    c_idx = add_node((fx, fy))
                    seg_arr[n, 0] = (ux, uy); seg_arr[n, 1] = (fx, fy)
                    splittable[n] = s_split; is_street[n] = s_street; n += 1
                    seg_arr[n, 0] = (fx, fy); seg_arr[n, 1] = (wx, wy)
                    splittable[n] = s_split; is_street[n] = s_street; n += 1
                    _bump(c_idx, 2)
                    far_idx = c_idx
                else:
                    far_idx = add_node((fx, fy))

                grow(1)
                seg_arr[n, 0] = (ax0, ay0); seg_arr[n, 1] = (fx, fy)
                splittable[n] = False; is_street[n] = True; n += 1
                _bump(anchor_idx, 1)
                _bump(far_idx, 1)
            continue  # stub iteration done (whether or not it was placed)

        if n_nodes <= 3:
            continue
```

- [ ] **Step 6: Run the new behavior tests**

Run: `cd /home/lpsha/s154446/fractality && python -m pytest space_filling_model/tests/test_stub_noise.py -k "stub or reproducible" -v`
Expected: PASS. (`test_p_stub_zero_is_identical` confirms the `p_stub>0.0` guard is fully short-circuited when off; the others confirm stubs appear, are bounded, and are deterministic.)

- [ ] **Step 7: Run the whole test file**

Run: `cd /home/lpsha/s154446/fractality && python -m pytest space_filling_model/tests/test_stub_noise.py -v`
Expected: PASS (all tests, including the Task 1 characterization).

- [ ] **Step 8: Commit**

```bash
cd /home/lpsha/s154446/fractality
git add space_filling_model/real_space_filling.py space_filling_model/tests/test_stub_noise.py
git commit -m "feat: emit stochastic dead-end stubs in _run (p_stub)"
```

---

## Task 5: Thread params through `run_circular_model`

Expose the three knobs on the public entry point and forward them to `_run`.

**Files:**
- Modify: `space_filling_model/real_space_filling.py`
- Modify: `space_filling_model/tests/test_stub_noise.py`

- [ ] **Step 1: Write the failing public-API test**

Append to `space_filling_model/tests/test_stub_noise.py`:

```python
def test_run_circular_model_accepts_stub_params():
    res = run_circular_model(area_m2=np.pi * 800.0 ** 2, model=2,
                             min_length=50.0, max_iter=1500, seed=7,
                             p_stub=0.25, stub_len_max=120.0, snap_tol=40.0,
                             verbose=False)
    assert len(res['segments']) > 0
```

- [ ] **Step 2: Run to verify failure**

Run: `cd /home/lpsha/s154446/fractality && python -m pytest space_filling_model/tests/test_stub_noise.py::test_run_circular_model_accepts_stub_params -v`
Expected: FAIL with `TypeError: run_circular_model() got an unexpected keyword argument 'p_stub'`.

- [ ] **Step 3: Add the params to `run_circular_model`'s signature**

Replace:

```python
def run_circular_model(city: Optional[str] = None, *, area_m2: Optional[float] = None,
                       model: int = 1, min_length: float = 50.0,
                       min_angle: float = np.pi / 4, max_iter: int = 200_000,
                       seed: int = 42, n_spokes: int = 8,
                       area_coeff: float = 2500.0,
                       area_scale: Optional[float] = None,
                       kdtree_rebuild_every: int = 256,
                       progress_every: int = 1000,
                       crs=None, verbose: bool = True):
```

with:

```python
def run_circular_model(city: Optional[str] = None, *, area_m2: Optional[float] = None,
                       model: int = 1, min_length: float = 50.0,
                       min_angle: float = np.pi / 4, max_iter: int = 200_000,
                       seed: int = 42, n_spokes: int = 8,
                       area_coeff: float = 2500.0,
                       area_scale: Optional[float] = None,
                       p_stub: float = 0.0,
                       stub_len_max: Optional[float] = None,
                       snap_tol: Optional[float] = None,
                       kdtree_rebuild_every: int = 256,
                       progress_every: int = 1000,
                       crs=None, verbose: bool = True):
```

- [ ] **Step 4: Forward the params in the `_run` call**

Replace:

```python
    segments = _run(seeds, blockers, area,
                    model=model, min_length=min_length, max_iter=max_iter,
                    seed=seed, min_angle=min_angle,
                    area_coeff=area_coeff, area_scale=area_scale,
                    kdtree_rebuild_every=kdtree_rebuild_every,
                    progress_every=progress_every, verbose=verbose)
```

with:

```python
    segments = _run(seeds, blockers, area,
                    model=model, min_length=min_length, max_iter=max_iter,
                    seed=seed, min_angle=min_angle,
                    area_coeff=area_coeff, area_scale=area_scale,
                    p_stub=p_stub, stub_len_max=stub_len_max, snap_tol=snap_tol,
                    kdtree_rebuild_every=kdtree_rebuild_every,
                    progress_every=progress_every, verbose=verbose)
```

- [ ] **Step 5: Update the module docstring knobs note**

In the file's top docstring, after the `min_length=50` m, `min_angle=π/4` sentence, add a sentence:

```
Set `p_stub` (default 0.0) to inject random dead-end stubs that populate the
low-degree shoulder of the dual-graph degree distribution; `stub_len_max`
(default 3*min_length) bounds their length and `snap_tol` (default min_length)
is the distance under which a stub T's into a nearby street.
```

- [ ] **Step 6: Run the public-API test**

Run: `cd /home/lpsha/s154446/fractality && python -m pytest space_filling_model/tests/test_stub_noise.py::test_run_circular_model_accepts_stub_params -v`
Expected: PASS.

- [ ] **Step 7: Run the full test file**

Run: `cd /home/lpsha/s154446/fractality && python -m pytest space_filling_model/tests/test_stub_noise.py -v`
Expected: PASS (all tests).

- [ ] **Step 8: Commit**

```bash
cd /home/lpsha/s154446/fractality
git add space_filling_model/real_space_filling.py space_filling_model/tests/test_stub_noise.py
git commit -m "feat: expose p_stub/stub_len_max/snap_tol on run_circular_model"
```

---

## Task 6: Verify the degree-distribution effect (manual)

Confirms the spec's success criterion. This is a manual visual check, not an automated test.

**Files:** none modified.

- [ ] **Step 1: Sweep `p_stub` and regenerate the dual-graph degree plot**

Edit the call at the top of `space_filling_model/plot_real_space_filling.py` (line ~17) to pass a `p_stub`, e.g.:

```python
res = run_circular_model(area_m2=np.pi*5000**2, area_coeff=1_000, model=2,
                         min_length=50.0, seed=42, p_stub=0.2)
```

Run: `cd /home/lpsha/s154446/fractality/space_filling_model && python plot_real_space_filling.py`
Re-run for `p_stub ∈ {0.0, 0.1, 0.2, 0.3}`, saving `plots/real_space_filling_degree_dist_2.png` under different names.

- [ ] **Step 2: Confirm the success criterion**

Inspect the `p(k)` plots and confirm:
1. The low-`k` shoulder (k ≈ 1–4) rises as `p_stub` increases.
2. The high-`k` tail slope (k > 4) stays roughly constant.
3. The raw-lines plot (`real_space_filling_raw_2.png`) shows short dangling stub segments near the periphery and intersections.

If the tail visibly steepens/flattens with `p_stub`, the stubs are leaking into the hub population — re-check that stub segments are created with `splittable[n] = False` (Task 4, Step 5).

- [ ] **Step 3 (optional): Revert the plot-script edit**

If the `plot_real_space_filling.py` edit was only for the sweep, revert it (or leave `p_stub` as a documented default) and commit if kept.

---

## Self-review notes

- **Spec coverage:** branch-or-stub in main loop (Task 4) ✓; mixed anchor incl. degree-≥3 intersection nodes (Task 2 bookkeeping + Task 4 anchor pick) ✓; optional T-in snap with segment split (Task 4 Step 5) ✓; non-splittable stubs still returned (Task 2 `is_street`) ✓; three params with `p_stub=0.0` backward-compat default (Tasks 4–5) ✓; `min_angle`/blocking reuse (Task 4) ✓; verification sweep (Task 6) ✓.
- **Type/name consistency:** `_nearest_segment(fx, fy, seg_arr, max_dist, exclude_xy=None)` defined in Task 3 and called with those exact args in Task 4. `_bump`, `deg3_idx`, `is_street`, `node_deg` introduced in Task 2 and used consistently thereafter. Imports `_segment_blocks`, `_min_incident_angle` added in Task 2 Step 1 before first use in Task 4.
- **Scope:** single subsystem (`real_space_filling.py`); no decomposition needed.
