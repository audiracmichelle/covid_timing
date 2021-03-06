#!/usr/bin/bash

#SBATCH -J var_relax
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --cpus-per-task 4
#SBATCH -p small
#SBATCH -o var_relax.%j.out
#SBATCH -e var_relax.%j.err
#SBATCH -t 2:59:00
#SBATCH -A A-ib1


#------------------------------------------------------

module load gcc/9
module load launcher

cd $WORK2/covid_timing/
export LAUNCHER_WORKDIR=$WORK2/covid_timing/
export LAUNCHER_JOB_FILE=$WORK2/covid_timing/tacc/jobs_variational_relax.sh

$LAUNCHER_DIR/paramrun