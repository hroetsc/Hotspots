#!/bin/bash

module purge
module load python/3.8.2 gcc/8.2.0 openmpi/gcc/64/3.1.4

source hot/bin/activate

module load cuda10.1/toolkit/10.1.105
module load cuda10.1/blas/10.1.105
module load cuda10.1/fft/10.1.105
module load cuda10.1/nsight/10.1.105
module load cuda10.1/profiler/10.1.105
module load cudnn/10.1v7.6.5

salloc -p gpu -C scratch2 -N 16 -n 16 --tasks-per-node 1 --gpus-per-task=1 --mem-per-gpu=20G -t 04:00:00 --job-name='validation'
scontrol show hostnames $SLURM_JOB_NODELIST

srun --mpi=pmix -o validation-%J.out python src/4_prediction.py

cp -rf /scratch2/hroetsc/Hotspots/results/best_model_prediction_*aa_*_rank*.csv results/
cp -rf /scratch2/hroetsc/Hotspots/results/last_model_prediction_*aa_*_rank*.csv results/
