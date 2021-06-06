#!/usr/bin/bash

#SBATCH -J variational
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --cpus-per-task 4
#SBATCH -p small
#SBATCH -o variational.%j.out
#SBATCH -e variational.%j.err
#SBATCH -t 2:59:00
#SBATCH -A A-ib1


#------------------------------------------------------

module load gcc/9
module load launcher

cd $WORK2/covid_timing/
export LAUNCHER_WORKDIR=$WORK2/covid_timing/
export LAUNCHER_JOB_FILE=$WORK2/covid_timing/tacc/jobs_variational.sh

$LAUNCHER_DIR/paramrun