"""E5: sufficiency test — rebuild the degree sequence from the decomposition.

Synthetic degrees are drawn using ONLY:
  (i)  the mixture weight w = P(endpoint-only),
  (ii) the empirical k_end distributions of the two classes,
  (iii) the empirical k_int distribution of through streets,
assuming k_end and k_int independent within a street. If the KS elbow at 4
re-emerges, the two-component mixture is sufficient to explain the transition.
"""
import sys
import numpy as np

sys.path.insert(0, '/home/lpsha/s154446/fractality/the_answer/scripts')
from e1_ks_elbow import value_counts, ks_curve
from e3_analyze import elbow

P = np.load('/home/lpsha/s154446/fractality/the_answer/results/e3_pooled.npz')
k, ki, ke = P['k'], P['k_int'], P['k_end']
rng = np.random.default_rng(0)
n = len(k)
w = float((ki == 0).mean())
n_epo = int(w * n)
syn_epo = rng.choice(ke[ki == 0], n_epo)
syn_thr = rng.choice(ke[ki > 0], n - n_epo) + rng.choice(ki[ki > 0], n - n_epo)
syn = np.concatenate([syn_epo, syn_thr])

for name, seq in (('real pooled k', k), ('synthetic mixture', syn)):
    v, c = value_counts(seq)
    _, D, _ = ks_curve(v, c)
    print(f'{name:18s} elbow={elbow(D)}  D(1..8): '
          + ' '.join(f'{x:.3f}' for x in D[:8]))
