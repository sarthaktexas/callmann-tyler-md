#!/bin/bash
# usage: bash run_seed.sh <system: sys_1to1|sys_pyronly> <seed> <xtc group: Other|PYR>
set -eo pipefail
export PATH=/opt/conda/bin:$PATH
SYS=$1; SEED=$2; GRP=$3
mkdir -p out
gmx --version | grep -E "GROMACS version|GPU support"; nvidia-smi -L
for f in nvt npt md; do sed "s/@GRP@/$GRP/" $f.mdp > ${f}_run.mdp; done
sed -i "s/^gen-seed.*/gen-seed = $SEED/" nvt_run.mdp
printf 'q\n' | gmx make_ndx -f $SYS.gro -o index.ndx >/dev/null 2>&1
G="-ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -update gpu"
gmx grompp -f em.mdp -c $SYS.gro -p $SYS.top -o em.tpr -maxwarn 1
gmx mdrun -deffnm em -ntmpi 1 -ntomp 8 -nb gpu
gmx grompp -f nvt_run.mdp -c em.gro -p $SYS.top -n index.ndx -o nvt.tpr -maxwarn 1
gmx mdrun -deffnm nvt $G
gmx grompp -f npt_run.mdp -c nvt.gro -t nvt.cpt -p $SYS.top -n index.ndx -o npt.tpr -maxwarn 1
gmx mdrun -deffnm npt $G
gmx grompp -f md_run.mdp -c npt.gro -t npt.cpt -p $SYS.top -n index.ndx -o md.tpr -maxwarn 1
gmx mdrun -deffnm md $G -cpt 5 -maxh ${MAXH:-1.6}
echo $GRP | gmx trjconv -s md.tpr -f $SYS.gro -n index.ndx -o solute_ref.gro >/dev/null 2>&1
cp nvt.xtc npt.xtc md.xtc solute_ref.gro em.gro npt.gro md.gro md.cpt nvt.edr npt.edr md.edr md.log index.ndx out/
ls -la out
