"""E4 analysis: compare model decompositions (noise levels) with the city pool."""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, '/home/lpsha/s154446/fractality/the_answer/scripts')
from e1_ks_elbow import value_counts, ks_curve, XMINS
from e3_analyze import elbow

RES = '/home/lpsha/s154446/fractality/the_answer/results'
MOD = os.path.join(RES, 'model_decomp')
FIG = '/home/lpsha/s154446/fractality/the_answer/figures'
BLUE, RED, GREY, GREEN = '#3b6fb6', '#c0392b', '0.45', '#27ae60'


def pk(seq):
    seq = seq[seq >= 1]
    v, c = np.unique(seq, return_counts=True)
    return v, c / c.sum()


def report(tag, d):
    k, ki, ke = d['k'], d['k_int'], d['k_end']
    pe, val = d['per_end_counts'], d['valence']
    v, D, _ = ks_curve(*value_counts(k))
    jv = val[val >= 3]
    out = dict(
        n=len(k), elbow=elbow(D), D=D,
        frac_epo=float((ki == 0).mean()),
        frac_kle4=float((k[k >= 1] <= 4).mean()),
        p_kgt4_epo=float((k[ki == 0] > 4).mean()) if (ki == 0).sum() else np.nan,
        end_le2=float((pe <= 2).mean()),
        val_ge5=float((jv >= 5).mean()),
        kend_le4=float((ke <= 4).mean()),
    )
    print(f"{tag:18s} n={out['n']:6d} elbow={out['elbow']} "
          f"frac_epo={out['frac_epo']:.3f} P(k<=4)={out['frac_kle4']:.3f} "
          f"P(k>4|epo)={out['p_kgt4_epo']:.3f} end<=2={out['end_le2']:.3f} "
          f"val>=5={out['val_ge5']:.3f}")
    if 'is_noise' in d:
        isn = d['is_noise']
        kn = k[isn]
        v_, c_ = np.unique(kn, return_counts=True)
        print(f"   noise streets: {isn.sum()} "
              f"({100*isn.mean():.0f}%), k<=4: {(kn<=4).mean():.4f}, hist:",
              {int(a): round(float(b/c_.sum()), 3) for a, b in zip(v_, c_)[:8]}
              if False else
              {int(a): round(float(b/c_.sum()), 3)
               for a, b in list(zip(v_, c_))[:8]})
        print(f"   epo among noise: {(ki[isn]==0).mean():.3f}; "
              f"share of k<=4 streets that are noise: "
              f"{float(isn[(k>=1)&(k<=4)].mean()):.3f}")
    return out


def main():
    files = sorted(f for f in os.listdir(MOD) if f.startswith('decomp_'))
    P = np.load(os.path.join(RES, 'e3_pooled.npz'))
    res = {}
    for f in files:
        tag = f[len('decomp_'):-4]
        res[tag] = (report(tag, np.load(os.path.join(MOD, f))),
                    np.load(os.path.join(MOD, f)))

    # city pool reference
    ck, cki = P['k'], P['k_int']
    cv, cD, _ = ks_curve(*value_counts(ck))
    print(f"cities pooled      elbow={elbow(cD)}")

    fig, axs = plt.subplots(2, 2, figsize=(10, 7.6))
    order = [t for t in ('m2_s42_nf0', 'm2_s42_nf0.5', 'm2_s42_nf1',
                         'm1_s42_nf0') if t in res]
    cols = {'m2_s42_nf0': RED, 'm2_s42_nf0.5': GREEN, 'm2_s42_nf1': BLUE,
            'm1_s42_nf0': 'orange'}
    labs = {'m2_s42_nf0': 'model 2, no noise',
            'm2_s42_nf0.5': 'model 2 + noise 0.5',
            'm2_s42_nf1': 'model 2 + noise 1.0',
            'm1_s42_nf0': 'model 1, no noise'}

    ax = axs[0, 0]
    v, p = pk(ck)
    ax.plot(v, p, '-', color='k', lw=2, alpha=0.5, label='cities (pooled)')
    for t in order:
        v, p = pk(res[t][1]['k'])
        ax.plot(v, p, 'o-', color=cols[t], ms=3, lw=1, label=labs[t])
    ax.axvline(4, color='k', ls=':', lw=1)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('$k$'); ax.set_ylabel('$p(k)$')
    ax.set_title('(a) degree distributions', loc='left', fontsize=11)
    ax.legend(frameon=False, fontsize=8)

    ax = axs[0, 1]
    ks = np.arange(1, 16)
    cf = [float((cki[ck == kk] == 0).mean()) for kk in ks]
    ax.plot(ks, cf, '-', color='k', lw=2, alpha=0.5, label='cities')
    for t in order:
        k, ki = res[t][1]['k'], res[t][1]['k_int']
        f = [float((ki[k == kk] == 0).mean()) if (k == kk).sum() > 10
             else np.nan for kk in ks]
        ax.plot(ks, f, 'o-', color=cols[t], ms=3, lw=1, label=labs[t])
    ax.axvline(4.5, color='k', ls=':', lw=1)
    ax.set_xlabel('$k$'); ax.set_ylabel('fraction endpoint-only')
    ax.set_title('(b) endpoint-only composition', loc='left', fontsize=11)
    ax.legend(frameon=False, fontsize=8)

    ax = axs[1, 0]
    ax.plot(XMINS, cD, '-', color='k', lw=2, alpha=0.5, label='cities')
    for t in order:
        ax.plot(XMINS, res[t][0]['D'], 'o-', color=cols[t], ms=3, lw=1,
                label=labs[t])
    ax.axvline(4, color='k', ls=':', lw=1)
    ax.set_yscale('log')
    ax.set_xlabel(r'$k_{\min}$'); ax.set_ylabel('KS $D$')
    ax.set_title('(c) KS scans', loc='left', fontsize=11)
    ax.legend(frameon=False, fontsize=8)

    ax = axs[1, 1]
    cjv = P['valence']; cjv = cjv[cjv >= 3]
    v, c = np.unique(cjv, return_counts=True)
    ax.bar(v - 0.2, c / c.sum(), width=0.38, color='k', alpha=0.5, ec='w',
           label='cities')
    t0 = order[0]
    jv = res[t0][1]['valence']; jv = jv[jv >= 3]
    v, c = np.unique(jv, return_counts=True)
    ax.bar(v[v <= 9] + 0.2, (c / c.sum())[v <= 9], width=0.38, color=RED,
           ec='w', label=labs[t0])
    ax.set_xlabel('junction valence'); ax.set_ylabel('fraction')
    ax.set_title('(d) junction valence', loc='left', fontsize=11)
    ax.legend(frameon=False, fontsize=8)

    fig.savefig(os.path.join(FIG, 'fig4_model.png'), bbox_inches='tight')
    print('saved fig4_model.png')


if __name__ == '__main__':
    main()
