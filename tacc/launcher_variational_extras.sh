#!/usr/bin/bash

#SBATCH -J extras
#SBATCH -N 1
#SBATCH -n 10
#SBATCH --cpus-per-task 4
#SBATCH -p small
#SBATCH -o extras.%j.out
#SBATCH -e extras.%j.err
#SBATCH -t 2:59:00
#SBATCH -A A-ib1


#------------------------------------------------------

module load gcc/9
module load launcher

cd $WORK2/covid_timing/
export LAUNCHER_WORKDIR=$WORK2/covid_timing/
export LAUNCHER_JOB_FILE=$WORK2/covid_timing/tacc/jobs_variational_extras.sh

$LAUNCHER_DIR/paramrun