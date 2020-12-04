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

salloc -p gpu -N 8 -n 8 --gpus-per-task=1 --mem-per-gpu=20G -t 12:00:00 --job-name='hotspots'
scontrol show hostnames $SLURM_JOB_NODELIST
sinfo show job $SLURM_JOB_ID

#zip -r /scratch2/hroetsc/Hotspots/data/cnt.zip /scratch2/hroetsc/Hotspots/data/cnt/*
sbcast /scratch2/hroetsc/Hotspots/data/cnt.zip /local/cnt.zip
srun --mpi=pmix -o hotspots-%J.out python src/6_fitModel.py

cp -rf /scratch2/hroetsc/Hotspots/results/model_metrics.txt results/
cp -rf /scratch2/hroetsc/Hotspots/results/model/* results/
cp -rf /scratch2/hroetsc/Hotspots/results/*_prediction_rank*.csv results/


#srun --mpi=pmix -o hotspots-prediction-%J.out python src/5_makePrediction.py
#mpirun --mca mpi_warn_on_fork 0 --output-filename hotspots-%J-%N.out python C_fitModel.py

