"""Render results/e3_per_city.json as a markdown appendix table."""
import json

ROWS = json.load(open(
    '/home/lpsha/s154446/fractality/the_answer/results/e3_per_city.json'))

hdr = ('| city | streets | ends ≤2 | k_end ≤4 | P(k>4 \\| epo) | '
       'P(epo \\| k≤4) | P(int \\| k>4) | elbow |')
sep = '|---|---|---|---|---|---|---|---|'
lines = [hdr, sep]
for r in sorted(ROWS, key=lambda r: -r['n']):
    lines.append(
        f"| {r['city']} | {r['n']:,} | {r['frac_end_le2']:.3f} | "
        f"{r['frac_kend_le4']:.3f} | {r['p_kgt4_given_epo']:.3f} | "
        f"{r['p_epo_given_kle4']:.3f} | {r['p_int_given_kgt4']:.3f} | "
        f"{r['elbow_all'] if r['elbow_all'] > 0 else '–'} |")

out = '/home/lpsha/s154446/fractality/the_answer/results/per_city_table.md'
with open(out, 'w') as fh:
    fh.write('# Per-city decomposition statistics (T_a=20°, R=10 m)\n\n')
    fh.write('\n'.join(lines) + '\n')
print('wrote', out)
