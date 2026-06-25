"""E1: Reproduce the KS-elbow at degree 4 from precomputed city degree data.

For every city (and pooled over cities) compute the KS distance D(kmin) of a
discrete power-law fit with fixed lower cutoff kmin, for kmin = 1..15.
The 'transition' phi is where D(kmin) drops / attains its minimum.
Repeats over the full precomputed (t_buffer, a_threshold) grid.
"""
import os
import sys
import numpy as np
from scipy.special import zeta
from scipy.optimize import minimize_scalar

DATA = '/home/lpsha/s154446/fractality/data/city_degrees'
OUT = '/home/lpsha/s154446/fractality/the_answer/results'
XMINS = np.arange(1, 16)


def value_counts(deg):
    deg = np.asarray(deg).astype(np.int64)
    deg = deg[deg >= 1]
    v, c = np.unique(deg, return_counts=True)
    return v, c


def fit_discrete_pl(v, c, xmin):
    """Discrete power-law MLE + KS distance for tail k >= xmin (Clauset 2009)."""
    m = v >= xmin
    v, c = v[m], c[m]
    n = c.sum()
    if n < 50 or len(v) < 3:
        return np.nan, np.nan, n
    slogk = float((c * np.log(v)).sum())

    def negll(alpha):
        return n * np.log(zeta(alpha, xmin)) + alpha * slogk

    res = minimize_scalar(negll, bounds=(1.01, 8.0), method='bounded')
    alpha = res.x
    ecdf = np.cumsum(c) / n
    tcdf = 1.0 - zeta(alpha, v + 1) / zeta(alpha, xmin)
    D = float(np.max(np.abs(ecdf - tcdf)))
    return alpha, D, n


def ks_curve(v, c):
    alphas, Ds, ns = [], [], []
    for xm in XMINS:
        a, d, n = fit_discrete_pl(v, c, xm)
        alphas.append(a); Ds.append(d); ns.append(n)
    return np.array(alphas), np.array(Ds), np.array(ns)


def main():
    os.makedirs(OUT, exist_ok=True)

    # ---- headline combo t10_a20: per-city + pooled --------------------------
    combo = 't10_a20'
    path = os.path.join(DATA, combo)
    cities = sorted(f for f in os.listdir(path) if f.endswith('.npz'))
    print(f'{combo}: {len(cities)} cities')

    per_city_D = np.full((len(cities), len(XMINS)), np.nan)
    per_city_alpha = np.full((len(cities), len(XMINS)), np.nan)
    pooled = {}
    for i, f in enumerate(cities):
        d = np.load(os.path.join(path, f))['degree']
        v, c = value_counts(d)
        a, D, n = ks_curve(v, c)
        per_city_alpha[i], per_city_D[i] = a, D
        for vv, cc in zip(v, c):
            pooled[vv] = pooled.get(vv, 0) + cc

    pv = np.array(sorted(pooled)); pc = np.array([pooled[x] for x in pv])
    pooled_alpha, pooled_D, pooled_n = ks_curve(pv, pc)

    argmins = XMINS[np.nanargmin(per_city_D, axis=1)]
    print('pooled D(kmin):')
    for xm, D, a in zip(XMINS, pooled_D, pooled_alpha):
        print(f'  kmin={xm:2d}  D={D:.4f}  alpha={a:.3f}')
    print('pooled argmin kmin* =', XMINS[np.nanargmin(pooled_D)])
    vals, cnts = np.unique(argmins, return_counts=True)
    print('per-city argmin distribution:', dict(zip(vals.tolist(), cnts.tolist())))

    def elbow(D):
        """kmin with the largest multiplicative drop D(k-1)/D(k)."""
        r = D[:-1] / D[1:]
        return int(XMINS[1:][np.nanargmax(r)])

    print('pooled elbow (largest log-drop) =', elbow(pooled_D))
    per_city_elbow = np.array([elbow(per_city_D[i]) for i in range(len(cities))])
    vals, cnts = np.unique(per_city_elbow, return_counts=True)
    print('per-city elbow distribution:', dict(zip(vals.tolist(), cnts.tolist())))

    # cross-check against the powerlaw package on one city
    import powerlaw
    d = np.load(os.path.join(path, cities[16]))['degree']  # Copenhagen
    d = d[d >= 1]
    v, c = value_counts(d)
    print(f'cross-check on {cities[16]}:')
    for xm in [2, 4, 8]:
        fit = powerlaw.Fit(d, discrete=True, xmin=xm, verbose=False)
        a, D, _ = fit_discrete_pl(v, c, xm)
        print(f'  xmin={xm}: powerlaw alpha={fit.power_law.alpha:.3f} '
              f'D={fit.power_law.D:.4f} | custom alpha={a:.3f} D={D:.4f}')

    # ---- full (t,a) grid, pooled ---------------------------------------------
    t_s = [5, 10, 50, 100]; a_s = [5, 10, 20, 30, 40]
    grid_D = {}
    grid_alpha = {}
    grid_argmin = {}
    for t in t_s:
        for a_th in a_s:
            cb = f't{t}_a{a_th}'
            p = os.path.join(DATA, cb)
            if not os.path.isdir(p):
                continue
            pooled = {}
            for f in os.listdir(p):
                if not f.endswith('.npz'):
                    continue
                d = np.load(os.path.join(p, f))['degree']
                v, c = value_counts(d)
                for vv, cc in zip(v, c):
                    pooled[vv] = pooled.get(vv, 0) + cc
            pv = np.array(sorted(pooled)); pc = np.array([pooled[x] for x in pv])
            al, D, n = ks_curve(pv, pc)
            grid_D[cb] = D; grid_alpha[cb] = al
            r = D[:-1] / D[1:]
            grid_argmin[cb] = int(XMINS[1:][np.nanargmax(r)])
            print(f'{cb}: elbow = {grid_argmin[cb]}, argmin = '
                  f'{int(XMINS[np.nanargmin(D)])}  (D at 1..6: '
                  + ' '.join(f'{x:.3f}' for x in D[:6]) + ')')

    np.savez(os.path.join(OUT, 'e1_ks_elbow.npz'),
             xmins=XMINS,
             cities=np.array([f[:-4] for f in cities]),
             per_city_D=per_city_D, per_city_alpha=per_city_alpha,
             pooled_D=pooled_D, pooled_alpha=pooled_alpha,
             grid_combos=np.array(list(grid_D.keys())),
             grid_D=np.array([grid_D[k] for k in grid_D]),
             grid_alpha=np.array([grid_alpha[k] for k in grid_D]),
             grid_argmin=np.array([grid_argmin[k] for k in grid_D]))
    print('saved', os.path.join(OUT, 'e1_ks_elbow.npz'))


if __name__ == '__main__':
    main()
