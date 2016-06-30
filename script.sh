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
export LD_LIBRARY_PATH=/clusterhome/nbudanur/local2/lib/
# run
/clusterhome/nbudanur/local2/bin/mpirun  --mca btl openib,self -np $TASKS ./dns.x
#mpirun --mca btl openib,self hello_c

