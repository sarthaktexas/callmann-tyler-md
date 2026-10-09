#!/bin/bash
#SBATCH --job-name=us_wham
#SBATCH --partition=compute1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=/work/wrz135/callmann-tyler-md/umbrella/pattern1/logs/%x.%j.out
#SBATCH --error=/work/wrz135/callmann-tyler-md/umbrella/pattern1/logs/%x.%j.err
# WHAM does not need a GPU, so it runs on the CPU partition compute1. It can also be run on the login node
# or locally (bash 03_wham.sh with WORKDIR edited): it only reads windows/us_XX.tpr and us_XX_pullx.xvg.
set -eo pipefail
source /work/wrz135/callmann-tyler-md/umbrella/pattern1/arc_env.sh
cd windows; mkdir -p ../wham
ls us_??.tpr | sort > ../wham/tpr-files.dat
ls us_??_pullx.xvg | sort > ../wham/pullx-files.dat
cd ../wham
# pass the file lists relative to windows/
sed -i 's#^#../windows/#' tpr-files.dat pullx-files.dat
EQ_PS=${EQ_PS:-2000}     # discard first 2 ns of each window
$GMX wham -it tpr-files.dat -ix pullx-files.dat -o profile.xvg -hist histo.xvg -unit kJ -b $EQ_PS -bsres bsres.xvg -nBootstrap 50 -bs-method b-hist -zprof0 3.0 -temp 300
echo "WHAM DONE"
