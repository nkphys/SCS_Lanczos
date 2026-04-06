#!/bin/bash
#SBATCH --job-name=Run_w_l
#SBATCH -A st-ianaffle-1
#SBATCH -t 48:00:00
#SBATCH --output=output.txt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=96000M
cd $SLURM_SUBMIT_DIR
module load gcc/9.4.0
module load intel-mkl/2020.4.304
date
./lanczos input_run.inp > Run_out.txt
date
