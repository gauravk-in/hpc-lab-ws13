#!/bin/bash

#SBATCH -o /home/hpc/t1221/lu26xum/workspace/hpc-lab-ws13/Asgn2/matrix_mul/output_myjob.%j.%N.out
#SBATCH -D /home/hpc/t1221/lu26xum/workspace/hpc-lab-ws13/Asgn2/matrix_mul
#SBATCH -J BPar_8
#SBATCH --get-user-env
#SBATCH --partition=snb
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
BATCH --mail-type=end
#SBATCH --mail-user=gmkukreja@gmail.com
#SBATCH --export=NONE
#SBATCH --time=08:00:00

source /etc/profile.d/modules.sh
cd /home/hpc/t1221/lu26xum/workspace/hpc-lab-ws13/Asgn2/matrix_mul
export OMP_NUM_THREADS=8

./dgemm 2048 output_2048.txt > ouput_flops_2048.txt
./dgemm 8192 output_8192.txt > ouput_flops_8192.txt
