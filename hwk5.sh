#!/bin/sh
#$ -V
#$ -cwd
#$ -S /bin/bash
#$ -N CS4379Hwk5
#$ -o $JOB_NAME.o$JOB_ID
#$ -e $JOB_NAME.e$JOB_ID
#$ -q omni
#$ -pe fill 8
#$ -P quanah

module load gnu openmpi
mpirun --machinefile machinefile.$JOB_ID -np $NSLOTS ./hwk5
