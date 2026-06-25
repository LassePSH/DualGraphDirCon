"""E3: Test the endpoint-bound hypothesis on the decomposed cities.

Key quantities, pooled and per city:
  - distribution of per-end distinct-neighbour counts (predict: <=2 for ~98%)
  - distribution of k_end (predict: bounded at ~4)
  - composition by k: fraction endpoint-only (k_int=0) streets (predict: ~1
    below 4, ~0 above)
  - causal KS test: elbow for all / through-only / endpoint-only sequences
  - junction valence distribution
"""
import os
import sys
import numpy as np

sys.path.insert(0, '/home/lpsha/s154446/fractality/the_answer/scripts')
from e1_ks_elbow import value_counts, ks_curve, XMINS

DEC = '/home/lpsha/s154446/fractality/the_answer/results/decomp_t10_a20'
OUT = '/home/lpsha/s154446/fractality/the_answer/results'


def elbow(D):
    r = D[:-1] / D[1:]
    if np.all(np.isnan(r)):
        return -1
    return int(XMINS[1:][np.nanargmax(r)])


def ks_of(deg):
    v, c = value_counts(deg)
    _, D, _ = ks_curve(v, c)
    return D


def main():
    files = sorted(f for f in os.listdir(DEC) if f.endswith('.npz'))
    print(f'{len(files)} cities decomposed\n')

    pool = {k: [] for k in ('k', 'k_int', 'k_end', 'k_unk', 'n_seg', 'length',
                            'n_dead', 'per_end_counts', 'valence')}
    rows = []
    for f in files:
        d = np.load(os.path.join(DEC, f))
        for key in pool:
            pool[key].append(d[key])
        k, k_int, k_end = d['k'], d['k_int'], d['k_end']
        pe = d['per_end_counts']
        val = d['valence']
        D_all = ks_of(k)
        D_thr = ks_of(k[k_int > 0])
        D_epo = ks_of(k[k_int == 0])
        rows.append(dict(
            city=f[:-4], n=len(k),
            frac_end_le2=float((pe <= 2).mean()),
            frac_kend_le4=float((k_end <= 4).mean()),
            p_kgt4_given_epo=float((k[k_int == 0] > 4).mean()),
            p_epo_given_kle4=float((k_int[(k <= 4) & (k >= 1)] == 0).mean()),
            p_int_given_kgt4=float((k_int[k > 4] > 0).mean()),
            val_ge5=float((val[val >= 3] >= 5).mean()),
            elbow_all=elbow(D_all), elbow_thr=elbow(D_thr),
            D_all=D_all, D_thr=D_thr, D_epo=D_epo,
        ))
        r = rows[-1]
        print(f"{r['city'][:28]:28s} n={r['n']:7d} "
              f"end<=2:{r['frac_end_le2']:.3f} kend<=4:{r['frac_kend_le4']:.3f} "
              f"P(k>4|epo):{r['p_kgt4_given_epo']:.3f} "
              f"P(epo|k<=4):{r['p_epo_given_kle4']:.3f} "
              f"P(int|k>4):{r['p_int_given_kgt4']:.3f} "
              f"elbow all:{r['elbow_all']} through:{r['elbow_thr']}")

    P = {key: np.concatenate(v) for key, v in pool.items()}
    k, k_int, k_end = P['k'], P['k_int'], P['k_end']
    pe, val = P['per_end_counts'], P['valence']

    print('\n=== POOLED ===')
    print('streets:', len(k), ' ends:', len(pe))
    v_, c_ = np.unique(pe, return_counts=True)
    print('per-end distinct neighbours:',
          {int(a): round(float(b / c_.sum()), 4) for a, b in zip(v_, c_)})
    v_, c_ = np.unique(k_end, return_counts=True)
    print('k_end distribution:',
          {int(a): round(float(b / c_.sum()), 4) for a, b in zip(v_, c_)})
    jv = val[val >= 3]
    v_, c_ = np.unique(jv, return_counts=True)
    print('junction valence (>=3):',
          {int(a): round(float(b / c_.sum()), 4) for a, b in zip(v_, c_)[:7]}
          if False else
          {int(a): round(float(b / c_.sum()), 4)
           for a, b in list(zip(v_, c_))[:7]})

    print('\ncomposition by k (frac endpoint-only = k_int==0):')
    for kk in range(1, 13):
        m = k == kk
        if m.sum() == 0:
            continue
        print(f'  k={kk:2d}: n={m.sum():8d}  frac_epo={float((k_int[m]==0).mean()):.4f}'
              f'  frac_1seg={float((P["n_seg"][m]==1).mean()):.4f}'
              f'  median_len={float(np.median(P["length"][m])):8.1f}')

    print('\ncausal KS test (pooled):')
    for name, seq in (('all k', k),
                      ('k, through-only (k_int>0)', k[k_int > 0]),
                      ('k, endpoint-only (k_int=0)', k[k_int == 0]),
                      ('k_int alone (k_int>=1)', k_int[k_int >= 1]),
                      ('k_int alone (through)', k_int[k_int >= 1]),
                      ('k_end alone (k_end>=1)', k_end[k_end >= 1])):
        D = ks_of(seq)
        print(f'  {name:28s} elbow={elbow(D)}  '
              'D(1..8)=' + ' '.join(f'{x:.3f}' for x in D[:8]))

    # tail exponents: does k_int carry the scale-free part?
    from e1_ks_elbow import fit_discrete_pl
    for name, seq, xm in (('k tail (kmin=5)', k, 5),
                          ('k_int tail (kmin=1)', k_int[k_int >= 1], 1),
                          ('k_int tail (kmin=2)', k_int[k_int >= 1], 2),
                          ('k_int tail (kmin=5)', k_int[k_int >= 1], 5)):
        v_, c_ = value_counts(seq)
        a, D, n = fit_discrete_pl(v_, c_, xm)
        print(f'  {name:22s} alpha={a:.3f} D={D:.4f} n={n}')

    # length coupling: degree vs length within each class
    print('\nmedian length by k and class (epo vs through):')
    for kk in (1, 2, 3, 4, 6, 8, 12, 20):
        m = k == kk
        me, mt = m & (k_int == 0), m & (k_int > 0)
        le = float(np.median(P['length'][me])) if me.sum() > 20 else np.nan
        lt = float(np.median(P['length'][mt])) if mt.sum() > 20 else np.nan
        print(f'  k={kk:2d}: epo n={me.sum():7d} len={le:9.1f} | '
              f'through n={mt.sum():7d} len={lt:9.1f}')

    np.savez(os.path.join(OUT, 'e3_pooled.npz'), **P)
    import json
    with open(os.path.join(OUT, 'e3_per_city.json'), 'w') as fh:
        json.dump([{kk: (vv.tolist() if isinstance(vv, np.ndarray) else vv)
                    for kk, vv in r.items()} for r in rows], fh, indent=1)
    print('\nsaved e3_pooled.npz / e3_per_city.json')


if __name__ == '__main__':
    main()
