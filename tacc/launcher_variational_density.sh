#!/usr/bin/bash

#SBATCH -J var_popdens
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --cpus-per-task 4
#SBATCH -p small
#SBATCH -o var_popdens.%j.out
#SBATCH -e var_popdens.%j.err
#SBATCH -t 1:59:00
#SBATCH -A A-ib1


#------------------------------------------------------

module load launcher

wdir=/work2/05863/mgarciat/frontera/covid_timing
cd $wdir
export LAUNCHER_WORKDIR=$wdir
export LAUNCHER_JOB_FILE=tacc/jobs_variational_density.sh
export PRE_VARS="college age_65_plus black hispanic popdensity"

$LAUNCHER_DIR/paramrun