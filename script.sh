#!/bin/bash
#$ -S /bin/bash
#$ -N dnsbox_test
#$  -pe openmpi 8
#$ -l h_vmem=2G
#$ -l h_rt=72:00:00
#$ -j y
#$ -q hybrid.q
#$ -cwd
ulimit -a
ulimit -s unlimited
TASKS=8
module load openmpi/1.4.5
# run
mpirun  --mca btl openib,self -np $TASKS ./dns.x
#mpirun --mca btl openib,self hello_c

