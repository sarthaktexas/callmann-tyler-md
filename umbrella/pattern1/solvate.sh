#!/bin/bash
# Solvate solute_unboxed.gro in a cubic box and write sys_pull.gro / sys_pull.top (TIP3P WAT, same ff as boxes/sys_1to1.top)
set -eo pipefail
BOX=${1:-7.0}
gmx editconf -f solute_unboxed.gro -o solute_box.gro -box $BOX $BOX $BOX -noc >/dev/null 2>&1
gmx solvate -cp solute_box.gro -cs spc216.gro -o solv_raw.gro -maxsol 100000 2>&1 | grep -E "Number of SOL|Generated"
# rename SOL -> WAT with the atom names used by the tleap topology (O, H1, H2)
python3 - <<'PY'
lines=open('solv_raw.gro').read().split('\n')
out=[lines[0]]; n=int(lines[1]); body=lines[2:2+n]; box=lines[2+n]
new=[]
for l in body:
    if l[5:10].strip()=='SOL':
        nm=l[10:15].strip(); nm={'OW':'O','HW1':'H1','HW2':'H2'}[nm]
        l=l[:5]+'WAT  '+f'{nm:>5}'+l[15:]
    new.append(l)
nw=sum(1 for l in new if l[5:10].strip()=='WAT' and l[10:15].strip()=='O')
open('sys_pull.gro','w').write('\n'.join([lines[0],lines[1]]+new+[box,''])); print('waters',nw)
open('nwat.txt','w').write(str(nw))
PY
NW=$(cat nwat.txt)
python3 - "$NW" <<'PY'
import sys
nw=sys.argv[1]; s=open('../../boxes/sys_1to1.top').read()
i=s.index('[ molecules ]')
open('sys_pull.top','w').write(s[:i]+'[ molecules ]\n; Compound       #mols\nPYR                  1\nsystem1              1\nWAT               %s\n'%nw)
PY
