#!/usr/bin/bash

#SBATCH -J variational
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --cpus-per-task 4
#SBATCH -p skx-normal
#SBATCH -o variational.%j.out
#SBATCH -e variational.%j.err
#SBATCH -t 4:29:00
#SBATCH -A A-ib1


#------------------------------------------------------

module load launcher

wdir=/work2/05863/mgarciat/frontera/covid_timing
cd $wdir
export LAUNCHER_WORKDIR=$wdir
export LAUNCHER_JOB_FILE=tacc/jobs_variational.sh

$LAUNCHER_DIR/paramrun