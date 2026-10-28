#!/usr/bin/env python3
"""Build / refresh the trajectory dashboard in site/ from whatever trajectories are on disk.

  python scripts/build_site.py                 # process every box whose md_<name>.xtc changed since the last build
  python scripts/build_site.py --watch 300     # keep going: re-check every 300 s (for while rsync is still bringing files in)
  python scripts/build_site.py --only many_p8_c32 --force

Open site/index.html in a browser (works from file://; or `python -m http.server -d site 8000`).

Per box, when boxes_*/md_<name>.xtc is new or its size/mtime changed (incremental; state in site/work/<name>/key.json):
  1. gmx trjconv -pbc whole -dt 500 (full system, water kept)          -> site/work/<name>/whole500.xtc
  2. scripts/site_metrics.py (waters at pyrene, Rg per chain, contacts, clusters)  -> metrics.csv
     with --sasa: pyrene SASA coverage (Shrake-Rupley on pyrene atoms, 0.14 nm probe, 960 dots, Bondi radii; matches gmx sasa
     via analyze_many.py to <0.2 percentage points on many_p8_c8)
  3. optional --gmx-sasa: scripts/analyze_many.py --dt 500 (gmx sasa) -> cov.csv, as a cross-check
  4. movie: reuse movies/<set>/<name>.mp4 if newer than the xtc and frame count matches, else render it
     (render_box_movie.py for C1, render_many_movie.py for many-pyrene), then transcode to site/videos/<name>.mp4 (720p, GOP 5 for fast seeking)
  5. contacts by chain part (scripts/site_contacts.py, boxes with pyrene and chains) -> contacts.csv, written to site/contacts.js
  6. run status from md_<name>.log (last time, Finished, ns/day) and any live rsync temp file (.md_<name>.xtc.XXXX)
Then site/data.js is rewritten (atomic) and site/index.html is copied from scripts/site_template/index.html.
CPU: one box at a time, OMP_NUM_THREADS=2 (the Mac is unstable under full load).
"""
import argparse, glob, json, os, re, shutil, subprocess, sys, time, traceback
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore', category=pd.errors.PerformanceWarning)

GMX_SASA = False
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ENVS = os.path.expanduser('~/.claude-science/conda/envs')
GMX = os.environ.get('GMX', f'{ENVS}/gmx-local/bin/gmx')
PY = os.environ.get('SITE_PY', f'{ENVS}/chem-rdkit/bin/python')
SITE = os.path.join(ROOT, 'site'); WORKD = os.path.join(SITE, 'work'); VID = os.path.join(SITE, 'videos')
MV = 2   # metrics version: bump to recompute metrics.csv for every box without redoing trjconv/movies
CV = 1   # contacts version (scripts/site_contacts.py -> contacts.csv, site/contacts.js); bump to recompute
ENV = dict(os.environ, OMP_NUM_THREADS='2', MKL_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', GMX=GMX, GMX_MAXBACKUP='-1')

def B(name, set_, folder, npyr, nchains, label, movie, pattern):
    return dict(name=name, set=set_, folder=folder, npyr=npyr, nchains=nchains, label=label, movie=movie, patterns=pattern)
BOXES = []
for n in range(1, 5):
    BOXES.append(B(f'sys_{n}to1', 'C1 stoichiometry (1 pyrene)', 'boxes_c1', 1, n, f'{n}:1 polymer:pyrene', f'movies/c1/c1_{n}to1.mp4', [f's{k}' for k in range(1, n + 1)]))
BOXES.append(B('sys_pyronly', 'C1 stoichiometry (1 pyrene)', 'boxes_c1', 1, 0, 'pyrene only (free reference)', 'movies/c1/c1_pyronly.mp4', []))
for n in range(1, 5):
    BOXES.append(B(f'chainonly_{n}', 'C1 chain-only controls', 'boxes_c1', 0, 1, f'chain s{n} alone', f'movies/c1/c1_chainonly_{n}.mp4', [f's{n}']))
for c in (0, 8, 16, 24, 32):
    BOXES.append(B(f'many_p8_c{c}', 'Many-pyrene (8 pyrene)', 'boxes_many', 8, c, f'8 pyrene + {c} chains', f'movies/many/many_p8_c{c}.mp4',
                   [f's{(k % 4) + 1}' for k in range(c)]))

def log(*m): print(time.strftime('%H:%M:%S'), *m, flush=True)
def run(cmd, cwd=ROOT, stdin=None, logf=None):
    r = subprocess.run(cmd, shell=True, cwd=cwd, env=ENV, input=stdin, capture_output=True, text=True)
    if logf: open(logf, 'a').write(f'$ {cmd}\n{r.stdout[-4000:]}\n{r.stderr[-4000:]}\n')
    if r.returncode: raise RuntimeError(f'{cmd[:120]} -> exit {r.returncode}: {r.stderr[-600:]}')
    return r

def paths(b):
    f = os.path.join(ROOT, b['folder']); n = b['name']
    gro = next((p for p in (f'{f}/md_{n}.gro', f'{f}/npt_{n}.gro', f'{f}/{n}.gro') if os.path.exists(p)), None)
    return dict(xtc=f'{f}/md_{n}.xtc', tpr=f'{f}/md_{n}.tpr', log=f'{f}/md_{n}.log', gro=gro, folder=f)

def run_status(b):
    p = paths(b); s = dict(log_time_ns=None, finished=False, ns_per_day=None, target_ns=None, syncing=False, sync_bytes=None)
    if os.path.exists(p['log']):
        txt = open(p['log'], errors='ignore').read()
        m = re.findall(r'^\s+Step\s+Time\s*\n\s+\d+\s+([\d.]+)', txt, re.M)
        if m: s['log_time_ns'] = float(m[-1]) / 1000
        s['finished'] = 'Finished mdrun' in txt
        m = re.findall(r'^Performance:\s+([\d.]+)', txt, re.M); s['ns_per_day'] = float(m[-1]) if m else None
        ns = re.search(r'^\s+nsteps\s+=\s+(\d+)', txt, re.M); dt = re.search(r'^\s+dt\s+=\s+([\d.]+)', txt, re.M)
        if ns and dt: s['target_ns'] = int(ns.group(1)) * float(dt.group(1)) / 1000
    tmps = glob.glob(os.path.join(p['folder'], f'.md_{b["name"]}.xtc.*'))
    live = [t for t in tmps if time.time() - os.path.getmtime(t) < 300]
    if live: s['syncing'] = True; s['sync_bytes'] = max(os.path.getsize(t) for t in live)
    if os.path.exists(p['xtc']): s['xtc_bytes'] = os.path.getsize(p['xtc'])
    return s

def xtc_key(b):
    x = paths(b)['xtc']
    return None if not os.path.exists(x) else dict(size=os.path.getsize(x), mtime=int(os.path.getmtime(x)))

def process(b, force=False):
    n = b['name']; p = paths(b); w = os.path.join(WORKD, n); os.makedirs(w, exist_ok=True); lg = f'{w}/build.log'
    key = xtc_key(b); kf = f'{w}/key.json'
    old = json.load(open(kf)) if os.path.exists(kf) else None
    if key is None: return 'no xtc'
    same_xtc = old is not None and all(old.get(k) == key[k] for k in ('size', 'mtime'))
    if same_xtc and old.get('mv') == MV and not force and os.path.exists(f'{w}/metrics.csv') and os.path.exists(f'{VID}/{n}.mp4'): return 'up to date'
    if same_xtc and old.get('mv') != MV and os.path.exists(f'{w}/metrics.csv'): os.remove(f'{w}/metrics.csv')
    if not same_xtc or force:
        for f in ('whole500.xtc', 'metrics.csv', 'cov.csv', 'video.json'):
            if os.path.exists(f'{w}/{f}'): os.remove(f'{w}/{f}')
        shutil.rmtree(f'{w}/cov_work', ignore_errors=True)
    t0 = time.time(); log(n, 'processing', key)
    if not os.path.exists(f'{w}/whole500.xtc'):
        run(f'"{GMX}" trjconv -s "{p["tpr"]}" -f "{p["xtc"]}" -pbc whole -dt 500 -o "{w}/whole500.xtc"', stdin='System\n', logf=lg)
    if not os.path.exists(f'{w}/metrics.csv'):
        run(f'"{PY}" scripts/site_metrics.py --sasa --gro "{p["gro"]}" --xtc "{w}/whole500.xtc" --npyr {b["npyr"]} --nchains {b["nchains"]} --out "{w}/metrics.csv"', logf=lg)
    if GMX_SASA and b['npyr'] and not os.path.exists(f'{w}/cov.csv'):   # optional cross-check with gmx sasa (slow for big boxes)
        run(f'"{PY}" scripts/analyze_many.py --tpr "{p["tpr"]}" --xtc "{w}/whole500.xtc" --npyr {b["npyr"]} --nchains {b["nchains"]} --dt 500 --out "{w}/cov"', logf=lg)
        for f in ('whole.xtc',):
            if os.path.exists(f'{w}/cov_work/{f}'): os.remove(f'{w}/cov_work/{f}')
    nfr = len(pd.read_csv(f'{w}/metrics.csv'))
    # movie: reuse if newer than the trajectory and frame-matched, else render
    import imageio_ffmpeg
    mv = os.path.join(ROOT, b['movie']); ok = False
    if os.path.exists(mv) and os.path.getmtime(mv) > key['mtime']:
        ok = imageio_ffmpeg.count_frames_and_secs(mv)[0] == nfr
    if not ok:
        log(n, f'rendering movie ({nfr} frames)')
        if b['folder'] == 'boxes_many':
            csv = f'--csv "{w}/metrics.csv"' if b['npyr'] else ''
            run(f'"{PY}" scripts/render_many_movie.py --gro "{p["gro"]}" --xtc "{w}/whole500.xtc" --npyr {b["npyr"]} --nchains {b["nchains"]} '
                f'--label "{b["label"]}" {csv} {"--fixed-zoom" if b["nchains"] == 0 else ""} --out "{mv}"', logf=lg)
        else:
            c1csv = os.path.join(ROOT, 'analysis/c1', f'{n}.csv')
            csv = f'--csv "{c1csv}"' if os.path.exists(c1csv) else ''
            run(f'"{PY}" scripts/render_box_movie.py --gro "{p["gro"]}" --xtc "{w}/whole500.xtc" --nchains {b["nchains"]} {"--no-pyrene" if not b["npyr"] else ""} '
                f'--label "{b["label"]}" {csv} --out "{mv}"', logf=lg)
    os.makedirs(VID, exist_ok=True); tmp = f'{VID}/.{n}.tmp.mp4'
    run(f'"{imageio_ffmpeg.get_ffmpeg_exe()}" -y -loglevel error -i "{mv}" -vf "scale=1280:-2" -c:v libx264 -preset veryfast -crf 26 -g 5 -keyint_min 5 '
        f'-pix_fmt yuv420p -movflags +faststart -an -threads 2 "{tmp}"', logf=lg)
    os.replace(tmp, f'{VID}/{n}.mp4')
    nv, secs = imageio_ffmpeg.count_frames_and_secs(f'{VID}/{n}.mp4')
    json.dump(dict(frames=nv, secs=secs, fps=25, source=os.path.relpath(mv, ROOT)), open(f'{w}/video.json', 'w'))
    json.dump(dict(key, mv=MV), open(kf, 'w')); log(n, f'done in {time.time() - t0:.0f} s ({nfr} frames, video {nv})')
    return 'processed'

def contacts_step(b):
    """Which part of the chain touches pyrene (site_contacts.py). Separate from process() so it never re-transcodes movies."""
    if not (b['npyr'] and b['nchains']): return None
    n = b['name']; p = paths(b); w = os.path.join(WORKD, n); kf = f'{w}/key.json'
    if not (os.path.exists(kf) and os.path.exists(f'{w}/whole500.xtc')): return None
    key = json.load(open(kf)); ck = f'{w}/contacts_key.json'
    want = dict(size=key.get('size'), mtime=key.get('mtime'), cv=CV)
    if os.path.exists(ck) and json.load(open(ck)) == want and os.path.exists(f'{w}/contacts.csv'): return 'up to date'
    t0 = time.time(); log(n, 'contacts by chain part')
    run(f'"{PY}" scripts/site_contacts.py --gro "{p["gro"]}" --xtc "{w}/whole500.xtc" --npyr {b["npyr"]} --nchains {b["nchains"]} '
        f'--out "{w}/contacts.csv" --meta "{w}/contacts_meta.json"', logf=f'{w}/build.log')
    json.dump(want, open(ck, 'w')); log(n, f'contacts done in {time.time() - t0:.0f} s'); return 'processed'

def write_contacts():
    out = {}
    for b in BOXES:
        w = os.path.join(WORKD, b['name'])
        if not (os.path.exists(f'{w}/contacts.csv') and os.path.exists(f'{w}/contacts_meta.json')): continue
        m = pd.read_csv(f'{w}/contacts.csv'); meta = json.load(open(f'{w}/contacts_meta.json'))
        meta['series'] = {k: [round(float(x), 2) for x in m[k]] for k in m.columns if k != 'covpoly_check'}
        out[b['name']] = meta
    js = 'window.CONTACTS = ' + json.dumps(out, separators=(',', ':')) + ';\n'
    tmp = f'{SITE}/.contacts.js.tmp'; open(tmp, 'w').write(js); os.replace(tmp, f'{SITE}/contacts.js')
    log('wrote site/contacts.js', f'{len(js) // 1024} KB, {len(out)} boxes')

def series_for(b):
    w = os.path.join(WORKD, b['name'])
    if not os.path.exists(f'{w}/metrics.csv'): return None
    m = pd.read_csv(f'{w}/metrics.csv'); m['tk'] = (m.time_ns * 100).round().astype(int)
    if os.path.exists(f'{w}/cov.csv') and 'covtot_mean' not in m.columns:
        c = pd.read_csv(f'{w}/cov.csv'); c['tk'] = (c.time_ns * 100).round().astype(int)
        keep = ['tk'] + [k for k in c.columns if k.startswith('covpoly') or k.startswith('covtot')]
        m = m.merge(c[keep], on='tk', how='left')
    m = m.drop(columns='tk').copy()
    out = {}
    for k in m.columns:
        v = m[k].to_numpy(dtype=float)
        if np.all(np.isnan(v)): continue
        out[k] = [None if np.isnan(x) else round(float(x), 4) for x in v]
    return out

def write_site():
    boxes = []
    for b in BOXES:
        w = os.path.join(WORKD, b['name']); s = series_for(b)
        vj = json.load(open(f'{w}/video.json')) if os.path.exists(f'{w}/video.json') else None
        key = json.load(open(f'{w}/key.json')) if os.path.exists(f'{w}/key.json') else None
        geo = {}
        g = paths(b)['gro']
        if g:
            with open(g, 'rb') as fh:
                fh.seek(-200, 2); bx = [float(x) for x in fh.read().decode().strip().splitlines()[-1].split()[:3]]
            V = bx[0] * bx[1] * bx[2] * 1e-24          # nm^3 -> L
            geo = dict(box_nm=[round(x, 2) for x in bx], pyr_mM=round(b['npyr'] / 6.02214076e23 / V * 1e3, 2) if b['npyr'] else None,
                       chain_mM=round(b['nchains'] / 6.02214076e23 / V * 1e3, 2) if b['nchains'] else None)
        st = run_status(b); st['processed_key'] = key; xk = xtc_key(b); st['stale'] = bool(key and xk and any(key.get(k) != xk[k] for k in ('size', 'mtime')))
        boxes.append(dict({k: b[k] for k in ('name', 'set', 'npyr', 'nchains', 'label', 'patterns')}, status=st, geo=geo, video=(f'videos/{b["name"]}.mp4' if vj else None),
                          video_info=vj, series=s))
    js = 'window.SITE = ' + json.dumps(dict(built=time.strftime('%Y-%m-%d %H:%M'), boxes=boxes), separators=(',', ':')) + ';\n'
    tmp = f'{SITE}/.data.js.tmp'; open(tmp, 'w').write(js); os.replace(tmp, f'{SITE}/data.js')
    tpl = os.path.join(ROOT, 'scripts/site_template/index.html')
    if os.path.exists(tpl): shutil.copy(tpl, f'{SITE}/index.html')
    open(f'{SITE}/.gitignore', 'w').write('# dashboard build output (regenerate with scripts/build_site.py)\n*\n')
    log('wrote site/data.js', f'{len(js) // 1024} KB')

def once(a):
    sel = [b for b in BOXES if not a.only or b['name'] in a.only]
    # boxes whose xtc is being replaced by a live rsync go last (they will likely be complete by then)
    sel.sort(key=lambda b: (run_status(b)['syncing'], b['folder'] != 'boxes_many'))
    for b in sel:
        st = run_status(b)
        if st['syncing'] and not a.include_syncing:
            log(b['name'], 'rsync in progress, skipped this pass'); continue
        try:
            r = process(b, a.force)
            if r == 'processed': write_site()
            elif r != 'up to date': log(b['name'], r)
        except Exception as e:
            log(b['name'], 'FAILED', e); traceback.print_exc()
    for b in sel:
        if run_status(b)['syncing'] and not a.include_syncing: continue
        try:
            if contacts_step(b) == 'processed': write_contacts()
        except Exception as e:
            log(b['name'], 'contacts FAILED', e); traceback.print_exc()
    write_contacts()
    write_site()

if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--watch', type=float, default=0); ap.add_argument('--only', nargs='*')
    ap.add_argument('--gmx-sasa', action='store_true', help='also run analyze_many.py (gmx sasa) as a cross-check; slow for big boxes')
    ap.add_argument('--force', action='store_true'); ap.add_argument('--include-syncing', action='store_true')
    a = ap.parse_args(); os.makedirs(WORKD, exist_ok=True); GMX_SASA = a.gmx_sasa
    while True:
        once(a)
        if not a.watch: break
        log(f'sleeping {a.watch:.0f} s'); time.sleep(a.watch)
