#!/usr/bin/bash

#SBATCH -J mcmc
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --cpus-per-task 4
#SBATCH -p small
#SBATCH -o mcmc.%j.out
#SBATCH -e mcmc.%j.err
#SBATCH -t 47:59:00
#SBATCH -A A-ib1


#------------------------------------------------------

module load gcc/9
module load launcher

cd $WORK2/covid_timing/
export LAUNCHER_WORKDIR=$WORK2/covid_timing/
export LAUNCHER_JOB_FILE=$WORK2/covid_timing/tacc/jobs_mcmc.sh

$LAUNCHER_DIR/paramrun