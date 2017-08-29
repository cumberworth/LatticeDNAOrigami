#!/bin/sh

# Name of job
#SBATCH -J %OUTPUTFILEBASE

# Walltime limit (hours:mins:secs)
#SBATCH -t %WALLTIME:00:00

# Nodes and procs
#SBATCH -N %NODES
#SBATCH -n %PROCS

# Standard error and out files
#SBATCH -o %OUTPUTFILEBASE.e
#SBATCH -e %OUTPUTFILEBASE.o

module unload gcc
module load gcc/6.2.0

echo "Starting job $SLURM_JOB_ID"
echo
echo "SLURM assigned me this node:"
cat $SLURM_JOB_NODELIST
echo

export LD_LIBRARY_PATH=~/lib:$LD_LIBRARY_PATH
export PATH=~/bin/$PATH

# Main job
latticeDNAOrigami -i %OUTPUTFILEBASE > %OUTPUTDIRFILE/%OUTPUTFILEBASE.out

# Copy results to slowscratch mirror
cp %OUTPUTDIRFILE/%OUTPUTFILEBASE.* /sharedscratch/amc226/projects/origami/calculations/17-03-29_us-validation

echo
echo "Job finished. SLURM details are:"
echo
qstat -f ${SLURM_JOB_ID}
echo
echo Finished at `date`
