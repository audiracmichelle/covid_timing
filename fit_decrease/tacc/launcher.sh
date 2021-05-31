!/bin/bash

#SBATCH -J launcher
#SBATCH -N 1
#SBATCH -n 3
#SBATCH --cpus-per-task 6
#SBATCH -p development
#SBATCH -o decrease.%j.out
#SBATCH -e decrease.%j.err
#SBATCH -t 48:00:00
#SBATCH -A A-ib1


#------------------------------------------------------

module load gcc/9
module load launcher

cd $WORK2/covid_timing/fit_decrease
export LAUNCHER_WORKDIR=$WORK2/covid_timing/fit_decrease
export LAUNCHER_JOB_FILE=$WORK2/covid_timing/fit_decrease/tacc/jobs.txt

$LAUNCHER_DIR/paramrun