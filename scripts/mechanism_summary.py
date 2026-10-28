#!/usr/bin/env python3
"""Summarise analysis/mechanism/* (from mechanism_analysis.py) plus site/work/<name>/metrics.csv into analysis/mechanism/summary.json
and analysis/mechanism/SUMMARY.md. Window: t >= 20 ns unless stated. 500 ps frames.

Kinetics (pyrene x chain contact, min heavy-atom distance < 0.4 nm): contact gaps of one frame (0.5 ns) are bridged before counting, so
'binding' and 'unbinding' are events that last >= 1 ns. Residence = length of a bound segment that starts and ends inside the run
(segments cut by the start or end of the run are reported separately as 'persistent').
Autocorrelation: integrated time tau_int with Sokal's self-consistent window (c = 5); n_eff = T / (2 tau_int); SE = sd / sqrt(n_eff).
"""
import glob, json, os
import numpy as np, pandas as pd

M = 'analysis/mechanism'; W = 'site/work'; T0 = 20.0
BOX = {'sys_1to1': (1, 1, ['s1']), 'sys_2to1': (1, 2, ['s1', 's2']), 'sys_3to1': (1, 3, ['s1', 's2', 's3']), 'sys_4to1': (1, 4, ['s1', 's2', 's3', 's4']),
       'sys_pyronly': (1, 0, []), **{f'chainonly_{k}': (0, 1, [f's{k}']) for k in range(1, 5)},
       **{f'many_p8_c{c}': (8, c, [f's{(k % 4) + 1}' for k in range(c)]) for c in (0, 8, 16, 24, 32)}}

def tau_int(x, dt):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if len(x) < 10 or x.std() == 0: return np.nan
    x = x - x.mean(); n = len(x); f = np.fft.rfft(x, 2 * n); acf = np.fft.irfft(f * np.conj(f))[:n]; acf /= acf[0]
    tau = 0.5
    for m in range(1, n):
        tau += acf[m]
        if m >= 5 * tau: break
    return max(tau, 0.5) * dt
def stat(x, dt):
    x = np.asarray(x, float); x = x[np.isfinite(x)]
    if not len(x): return dict(mean=None)
    ta = tau_int(x, dt); T = len(x) * dt; neff = T / (2 * ta) if ta and np.isfinite(ta) else np.nan
    return dict(mean=float(x.mean()), sd=float(x.std(ddof=1)) if len(x) > 1 else None, tau_ns=None if not np.isfinite(ta) else float(ta),
                n_eff=None if not np.isfinite(neff) else float(neff), se=None if not np.isfinite(neff) else float(x.std(ddof=1) / np.sqrt(max(neff, 1))))
def bridge(b):
    b = b.copy()
    for i in range(1, len(b) - 1):
        if not b[i] and b[i - 1] and b[i + 1]: b[i] = True
    return b
def segments(b):
    out, i = [], 0
    while i < len(b):
        if b[i]:
            j = i
            while j + 1 < len(b) and b[j + 1]: j += 1
            out.append((i, j)); i = j + 1
        else: i += 1
    return out
J = lambda v: None if v is None or (isinstance(v, float) and not np.isfinite(v)) else (float(v) if isinstance(v, (float, np.floating)) else v)

S = {'window_ns': T0, 'boxes': {}, 'pattern': {}}
free_hyd = {}
for name, (NP, NC, pat) in BOX.items():
    ts_f = f'{M}/{name}_ts.csv'
    if not os.path.exists(ts_f): continue
    ts = pd.read_csv(ts_f); t = ts.time_ns.values; dt = float(np.median(np.diff(t))); w = t >= T0
    met = pd.read_csv(f'{W}/{name}/metrics.csv'); assert np.allclose(met.time_ns.values, t, atol=0.01)
    B = dict(npyr=NP, nchains=NC, frames=int(len(t)), t_end=float(t[-1]))
    # time series means with autocorrelation-corrected SE
    for col in ['agg_size', 'agg_rg', 'agg_kappa2', 'n_aggregates', 'hb_inter', 'hb_intra', 'hyd_core', 'hyd_glc', 'pp_stacked']:
        if col in ts: B[col] = stat(ts[col].values[w], dt)
    if NC and 'hb_inter' in ts: B['hb_inter_per_chain'] = stat(ts.hb_inter.values[w] / NC, dt); B['hb_intra_per_chain'] = stat(ts.hb_intra.values[w] / NC, dt)
    if NP:
        dcols = [c for c in ts if c.startswith('depth_pyr')]; scols = [c for c in ts if c.startswith('stack_imide_pyr')]
        if NC: B['depth'] = stat(np.nanmean(ts[dcols].values[w], axis=1), dt); B['depth_per_pyrene'] = [J(np.nanmean(ts[c].values[w])) for c in dcols]
        if NC: B['stack_imide_per_pyrene'] = stat(ts[scols].values[w].mean(1), dt); B['frac_pyrene_with_stacked_imide'] = J(float((ts[scols].values[w] > 0).mean()))
        cov = met['covpoly_mean'] if 'covpoly_mean' in met else met.get('covtot_mean')
        if cov is not None: B['coverage'] = stat(cov.values[w], dt)
        wp = met[[c for c in met if c.startswith('waters_pyr')]].values.mean(1); B['waters_per_pyrene'] = stat(wp[w], dt)
    # contact chemistry
    cf = f'{M}/{name}_contacts.csv'
    if os.path.exists(cf): B['contacts'] = pd.read_csv(cf).to_dict('list')
    pf = f'{M}/{name}_profiles.csv'
    if os.path.exists(pf): B['profiles'] = {k: [J(x) for x in v] for k, v in pd.read_csv(pf).to_dict('list').items()}
    # kinetics + pairs
    z = np.load(f'{M}/{name}_pairs.npz', allow_pickle=True)
    if NP and NC:
        pr = z['pairs']                                                    # frames x NP x NC
        res, pers, nb, nu, first = [], 0, 0, 0, []
        for i in range(NP):
            any_b = pr[:, i, :].any(1); first.append(float(t[np.argmax(any_b)]) if any_b.any() else None)
            for k in range(NC):
                b = bridge(pr[:, i, k]); sg = segments(b)
                for (s0, s1) in sg:
                    if s0 > 0: nb += 1
                    if s1 < len(b) - 1: nu += 1
                    if s0 > 0 and s1 < len(b) - 1: res.append((s1 - s0 + 1) * dt)
                    else: pers += 1
        ever = pr.any(0)                                                   # NP x NC
        late = pr[w]
        B['kinetics'] = dict(binding_events=nb, unbinding_events=nu, events_per_100ns=(nb + nu) / (t[-1] - t[0]) * 100, residence_ns_median=J(float(np.median(res))) if res else None,
                             residence_ns_max=J(float(max(res))) if res else None, n_complete_segments=len(res), n_persistent_segments=pers,
                             pairs_ever_in_contact=int(ever.sum()), pairs_total=int(NP * NC), distinct_chains_per_pyrene=float(ever.sum(1).mean()),
                             chains_ever_touching_pyrene=int(ever.any(0).sum()), first_contact_ns=first,
                             contacts_at_t0_kept_to_end=J(float((pr[w][0] & pr[-1]).sum() / max(pr[w][0].sum(), 1))))
        # per-chain binding by cis/trans pattern
        rgc = [c for c in met if c.startswith('rg_chain')]
        for k in range(NC):
            p_ = pat[k]; d = S['pattern'].setdefault(p_, dict(bound_frac=[], rg=[], boxes=[]))
            d['bound_frac'].append(float(late[:, :, k].any(1).mean())); d['rg'].append(float(met[rgc[k]].values[w].mean())); d['boxes'].append(name)
    if NP > 1:
        pp = z['pp'][w]; ag = z['agg_of_pyr'][w] if len(z['agg_of_pyr']) else None
        cp = met[[f'covpoly_pyr{i + 1}' for i in range(NP)]].values[w] if NC else None
        n = same = enc = 0
        for f in range(len(pp)):
            for i in range(NP):
                for j in range(i):
                    if not pp[f, i, j]: continue
                    n += 1
                    if ag is not None and ag[f, i] >= 0 and ag[f, i] == ag[f, j]: same += 1
                    if cp is not None and cp[f, i] > 70 and cp[f, j] > 70: enc += 1
        B['pyrene_pairs'] = dict(pair_contacts_per_frame=n / len(pp), frac_in_same_aggregate=J(same / n) if n and NC else None,
                                 frac_both_covered_gt70=J(enc / n) if n and NC else None, stacked_pairs_per_frame=J(float(z['pp'][w].shape[0] and pd.read_csv(ts_f).pp_stacked.values[w].mean())))
    if name.startswith('chainonly') and 'hyd_core' in ts: free_hyd[name] = (ts.hyd_core.values[w].mean(), ts.hyd_glc.values[w].mean())
    S['boxes'][name] = B
if free_hyd:
    S['free_chain_hydration'] = dict(core=float(np.mean([v[0] for v in free_hyd.values()])), glucose=float(np.mean([v[1] for v in free_hyd.values()])))
for p_, d in S['pattern'].items():
    d.update(n_chains=len(d['rg']), bound_frac_mean=float(np.mean(d['bound_frac'])), bound_frac_sd=float(np.std(d['bound_frac'], ddof=1)) if len(d['rg']) > 1 else None,
             rg_mean=float(np.mean(d['rg'])), rg_sd=float(np.std(d['rg'], ddof=1)) if len(d['rg']) > 1 else None)
json.dump(S, open(f'{M}/summary.json', 'w'), default=lambda o: J(o) if not isinstance(o, np.ndarray) else o.tolist(), indent=0)

# ---- markdown summary
f = lambda d, k=2: '–' if not d or d.get('mean') is None else f"{d['mean']:.{k}f} ± {d['se']:.{k}f}" if d.get('se') is not None else f"{d['mean']:.{k}f}"
L = ['# Mechanism analyses (t >= 20 ns; ± = SE corrected for autocorrelation)', '',
     '| box | coverage % (tau ns, n_eff) | pyrene depth r/R | imide stacked / pyrene | chains in contact ever (of N) | bind / unbind events | median residence ns | agg Rg nm | kappa2 | inter-chain HB / chain | core waters / atom |',
     '|---|---|---|---|---|---|---|---|---|---|---|']
for n, B in S['boxes'].items():
    k = B.get('kinetics', {}); c = B.get('coverage') or {}
    L.append(f"| {n} | {f(c, 1)} ({c.get('tau_ns') or 0:.1f}, {c.get('n_eff') or 0:.0f}) | {f(B.get('depth'))} | {f(B.get('stack_imide_per_pyrene'))} | "
             f"{k.get('chains_ever_touching_pyrene', '–')} / {B['nchains']} | {k.get('binding_events', '–')} / {k.get('unbinding_events', '–')} | {k.get('residence_ns_median') or '–'} | "
             f"{f(B.get('agg_rg'))} | {f(B.get('agg_kappa2'))} | {f(B.get('hb_inter_per_chain'))} | {f(B.get('hyd_core'))} |")
L += ['', '## Contact chemistry (enrichment = fraction of contacting chain atoms / fraction of chain heavy atoms)', '', '| box | core | linker | glucose | backbone |', '|---|---|---|---|---|']
for n, B in S['boxes'].items():
    if 'contacts' in B: L.append(f"| {n} | " + ' | '.join(f"{e:.2f}" for e in B['contacts']['enrichment']) + ' |')
L += ['', '## cis/trans pattern (chains in pyrene boxes, t >= 20 ns)', '', '| pattern | chains | fraction of frames touching a pyrene | Rg nm |', '|---|---|---|---|']
for p_ in sorted(S['pattern']):
    d = S['pattern'][p_]; L.append(f"| {p_} | {d['n_chains']} | {d['bound_frac_mean']:.2f} ± {d['bound_frac_sd'] or 0:.2f} (sd) | {d['rg_mean']:.3f} ± {d['rg_sd'] or 0:.3f} (sd) |")
L += ['', '## Pyrene pairs (boxes with 8 pyrene)', '', '| box | pyrene-pyrene contacts / frame | stacked pairs / frame | in same polymer aggregate | both >70% covered by polymer |', '|---|---|---|---|---|']
for n, B in S['boxes'].items():
    q = B.get('pyrene_pairs')
    if q: L.append(f"| {n} | {q['pair_contacts_per_frame']:.2f} | {q['stacked_pairs_per_frame'] or 0:.2f} | {('%.0f%%' % (100 * q['frac_in_same_aggregate'])) if q['frac_in_same_aggregate'] is not None else '–'} | {('%.0f%%' % (100 * q['frac_both_covered_gt70'])) if q['frac_both_covered_gt70'] is not None else '–'} |")
if 'free_chain_hydration' in S: L += ['', f"Free-chain hydration (chain-only controls): core {S['free_chain_hydration']['core']:.2f}, glucose {S['free_chain_hydration']['glucose']:.2f} waters per heavy atom."]
open(f'{M}/SUMMARY.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
