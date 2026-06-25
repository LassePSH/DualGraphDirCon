"""Figures for the phi=4 report (E1-E3). Model figure made separately."""
import os
import sys
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, '/home/lpsha/s154446/fractality/the_answer/scripts')
from e1_ks_elbow import value_counts, ks_curve, XMINS

RES = '/home/lpsha/s154446/fractality/the_answer/results'
FIG = '/home/lpsha/s154446/fractality/the_answer/figures'

plt.rcParams.update({'font.size': 11, 'axes.labelsize': 12,
                     'figure.dpi': 150, 'savefig.bbox': 'tight',
                     'xtick.direction': 'in', 'ytick.direction': 'in'})
BLUE, RED, GREY = '#3b6fb6', '#c0392b', '0.45'


def pk(seq, kmax=None):
    seq = seq[seq >= 1]
    v, c = np.unique(seq, return_counts=True)
    return v, c / c.sum()


def fig1_ks():
    d = np.load(os.path.join(RES, 'e1_ks_elbow.npz'), allow_pickle=True)
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.4))

    ax = axs[0]
    ax.plot(d['xmins'], d['pooled_D'], 'o-', color=BLUE, ms=4)
    ax.axvline(4, color='k', ls=':', lw=1)
    ax.set_yscale('log')
    ax.set_xlabel(r'candidate cutoff $k_{\min}$')
    ax.set_ylabel('KS distance $D$')
    ax.set_title('(a) pooled KS scan, $T_a$=20°, $R$=10 m', loc='left', fontsize=11)
    ax.annotate(r'$\phi=4$', (4.2, d['pooled_D'][1]), fontsize=12)

    ax = axs[1]
    per_D = d['per_city_D']
    r = per_D[:, :-1] / per_D[:, 1:]
    elbows = []
    for i in range(per_D.shape[0]):
        ri = r[i]
        elbows.append(XMINS[1:][np.nanargmax(ri)] if not np.all(np.isnan(ri))
                      else -1)
    elbows = np.array(elbows)
    ax.hist(elbows, bins=np.arange(0.5, 12.5), color=BLUE, ec='w')
    ax.axvline(4, color='k', ls=':', lw=1)
    ax.set_xlabel(r'per-city elbow of $D(k_{\min})$')
    ax.set_ylabel('# cities')
    ax.set_title(f'(b) 97 cities: {np.sum(elbows==4)} have elbow at 4',
                 loc='left', fontsize=11)

    ax = axs[2]
    for i, cb in enumerate(d['grid_combos']):
        ax.plot(d['xmins'], d['grid_D'][i], '-', color=GREY, alpha=0.5, lw=1)
    ax.plot(d['xmins'], d['pooled_D'], 'o-', color=BLUE, ms=4, zorder=5,
            label='$T_a$=20°, $R$=10 m')
    ax.axvline(4, color='k', ls=':', lw=1)
    ax.set_yscale('log')
    ax.set_xlabel(r'candidate cutoff $k_{\min}$')
    ax.set_ylabel('KS distance $D$')
    ax.set_title('(c) all 20 $(T_a, R)$ combinations', loc='left', fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    fig.savefig(os.path.join(FIG, 'fig1_ks_elbow.png'))
    plt.close(fig)


def fig2_decomposition():
    P = np.load(os.path.join(RES, 'e3_pooled.npz'))
    k, k_int, k_end = P['k'], P['k_int'], P['k_end']
    pe, val = P['per_end_counts'], P['valence']

    fig, axs = plt.subplots(2, 2, figsize=(10, 7.6))

    # (a) mixture decomposition of p(k)
    ax = axs[0, 0]
    n_tot = (k >= 1).sum()
    v, p = pk(k)
    ax.plot(v, p, 'o-', color='k', ms=4, lw=1, label='all streets')
    epo, thr = k[k_int == 0], k[k_int > 0]
    v, c = np.unique(epo[epo >= 1], return_counts=True)
    ax.plot(v, c / n_tot, 's-', color=BLUE, ms=4, lw=1,
            label='endpoint-only ($k_{int}=0$)')
    v, c = np.unique(thr[thr >= 1], return_counts=True)
    ax.plot(v, c / n_tot, '^-', color=RED, ms=4, lw=1,
            label='through ($k_{int}\\geq 1$)')
    ax.axvline(4, color='k', ls=':', lw=1)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('$k$'); ax.set_ylabel('$p(k)$')
    ax.set_title('(a) $p(k)$ is a mixture of two classes', loc='left',
                 fontsize=11)
    ax.legend(frameon=False, fontsize=9)

    # (b) composition fraction vs k
    ax = axs[0, 1]
    ks = np.arange(1, 16)
    frac = [float((k_int[k == kk] == 0).mean()) if (k == kk).sum() else np.nan
            for kk in ks]
    ax.bar(ks, frac, color=BLUE, ec='w')
    ax.axvline(4.5, color='k', ls=':', lw=1)
    ax.set_xlabel('$k$')
    ax.set_ylabel('fraction endpoint-only')
    ax.set_title('(b) bounded class dies at $k=4$', loc='left', fontsize=11)

    # (c) k_end vs k_int distributions
    ax = axs[1, 0]
    v, c = np.unique(k_end, return_counts=True)
    ax.plot(v[v >= 1], (c / c.sum())[v >= 1], 's-', color=BLUE, ms=4, lw=1,
            label='$k_{end}$ (endpoint-attached)')
    v, c = np.unique(k_int[k_int >= 1], return_counts=True)
    ax.plot(v, c / c.sum(), '^-', color=RED, ms=4, lw=1,
            label='$k_{int}$ (interior-attached)')
    ax.axvline(4, color='k', ls=':', lw=1)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('contribution to $k$'); ax.set_ylabel('$p$')
    ax.set_title('(c) endpoint part bounded, interior part scale-free',
                 loc='left', fontsize=11)
    ax.legend(frameon=False, fontsize=9)

    # (d) per-end neighbours + junction valence
    ax = axs[1, 1]
    v, c = np.unique(pe, return_counts=True)
    ax.bar(v - 0.2, c / c.sum(), width=0.38, color=BLUE, ec='w',
           label='distinct streets per street-end')
    jv = val[val >= 3]
    v, c = np.unique(jv, return_counts=True)
    ax.bar(v + 0.2, c / c.sum(), width=0.38, color=GREY, ec='w',
           label='junction valence ($\\geq$3)')
    ax.set_xlim(-0.8, 8)
    ax.set_xlabel('count'); ax.set_ylabel('fraction')
    ax.set_title('(d) ends see $\\leq 2$ streets; junctions are 3/4-valent',
                 loc='left', fontsize=11)
    ax.legend(frameon=False, fontsize=9)

    fig.savefig(os.path.join(FIG, 'fig2_decomposition.png'))
    plt.close(fig)


def fig3_length():
    P = np.load(os.path.join(RES, 'e3_pooled.npz'))
    k, k_int, L = P['k'], P['k_int'], P['length']
    fig, ax = plt.subplots(figsize=(5, 3.6))
    ks = np.unique(k[(k >= 1) & (k <= 200)])
    for m, col, lab, mk in ((k_int == 0, BLUE, 'endpoint-only', 's'),
                            (k_int > 0, RED, 'through', '^')):
        med = [np.median(L[m & (k == kk)]) if (m & (k == kk)).sum() > 30
               else np.nan for kk in ks]
        ax.plot(ks, med, mk, color=col, ms=4, label=lab)
    ax.axvline(4, color='k', ls=':', lw=1)
    ax.set_xscale('log'); ax.set_yscale('log')
    ax.set_xlabel('$k$'); ax.set_ylabel('median street length [m]')
    ax.set_title('length-degree coupling only for through streets',
                 fontsize=11)
    ax.legend(frameon=False, fontsize=9)
    fig.savefig(os.path.join(FIG, 'fig3_length.png'))
    plt.close(fig)


def fig5_schematic():
    """Why phi=4: a street has 2 ends; each end touches <=2 distinct streets."""
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.2))

    ax = axs[0]
    # endpoint-only street between a T (left) and an X (right)
    ax.plot([0, 4], [0, 0], color=BLUE, lw=3, solid_capstyle='round')
    ax.plot([0, 0], [-1.2, 1.2], color=RED, lw=2)            # T through street
    ax.plot([4, 4], [-1.2, 1.2], color=RED, lw=2)            # X cross street
    ax.plot([4, 5.2], [0, 0], color='orange', lw=2)          # opposite arm
    ax.annotate('1 neighbour\n(arms pair into\none through street)',
                (-0.1, 1.3), ha='center', fontsize=8)
    ax.annotate('2 neighbours\n(cross street + \nopposite arm)',
                (4.4, 1.3), ha='center', fontsize=8)
    ax.text(2, -1.7, r'endpoint-only street: $k = k_{end} \leq 2+2 = 4$',
            ha='center', fontsize=10)
    ax.set_xlim(-1.5, 6); ax.set_ylim(-2.2, 2.6); ax.axis('off')
    ax.set_title('(a) ends contribute at most $2+2 = 4$', loc='left',
                 fontsize=11)

    ax = axs[1]
    # through street crossing junctions
    ax.plot([0, 6], [0, 0], color=BLUE, lw=3, solid_capstyle='round')
    for x in (1, 2.5, 4, 5):
        ax.plot([x, x], [0, 1.1], color=RED, lw=2)
    ax.plot([2.5, 2.5], [-1.1, 0], color=RED, lw=2)
    ax.text(3, -1.7, r'through street: $k_{int} \propto$ junctions crossed'
            r' $\propto$ length', ha='center', fontsize=10)
    ax.set_xlim(-0.5, 6.5); ax.set_ylim(-2.2, 2.6); ax.axis('off')
    ax.set_title('(b) interior degree grows with length', loc='left',
                 fontsize=11)

    fig.savefig(os.path.join(FIG, 'fig5_schematic.png'))
    plt.close(fig)


if __name__ == '__main__':
    os.makedirs(FIG, exist_ok=True)
    fig1_ks()
    fig2_decomposition()
    fig3_length()
    fig5_schematic()
    print('figures written to', FIG)
