#!/usr/bin/bash

#SBATCH -J launcher
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --cpus-per-task 8
#SBATCH -p small
#SBATCH -o decrease.%j.out
#SBATCH -e decrease.%j.err
#SBATCH -t 47:59:00
#SBATCH -A A-ib1


#------------------------------------------------------

module load gcc/9
module load launcher

cd $WORK2/covid_timing/
export LAUNCHER_WORKDIR=$WORK2/covid_timing/
export LAUNCHER_JOB_FILE=$WORK2/covid_timing/tacc/jobs_variational.sh

$LAUNCHER_DIR/paramrun