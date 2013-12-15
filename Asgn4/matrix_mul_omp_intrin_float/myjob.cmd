#!/bin/bash

#SBATCH -o /home/hpc/t1221/lu26xum/workspace/hpc-lab-ws13/Asgn2/matrix_mul_intrinsic/output_myjob.%j.%N.out
#SBATCH -D /home/hpc/t1221/lu26xum/workspace/hpc-lab-ws13/Asgn2/matrix_mul_intrinsic
#SBATCH -J mm_int_par
BATCH --get-user-env
#SBATCH --partition=snb
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
BATCH --mail-type=end
#SBATCH --mail-user=gmkukreja@gmail.com
#SBATCH --export=NONE
#SBATCH --time=08:00:00

source /etc/profile.d/modules.sh
cd /home/hpc/t1221/lu26xum/workspace/hpc-lab-ws13/Asgn2/matrix_mul_intrinsic

make run_test
