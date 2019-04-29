#!/bin/sh

# Name of job
#SBATCH -J %OUTPUTFILEBASE

# Walltime limit (hours:mins:secs)
#SBATCH -t %WALLTIME:00:00

# Nodes and procs
#SBATCH -n %PROCS

# Standard error and out files
#SBATCH -o outs/%OUTPUTFILEBASE.o
#SBATCH -e outs/%OUTPUTFILEBASE.e

module unload gcc
module load gcc/6.2.0

echo "Starting job $SLURM_JOB_ID"

export LD_LIBRARY_PATH=~/lib:$LD_LIBRARY_PATH
export PATH=~/bin:$PATH
mpirun -n %PROCS mkdir -p %OUTPUTFILEDIR

# Main job
mpirun -n %NUMREPS latticeDNAOrigami -i %INPDIR/%OUTPUTFILEBASE.inp > %OUTPUTFILEDIR/%OUTPUTFILEBASE.out

# Copy results to sharedscratch
targetdir=$(pwd | sed "s:home:sharedscratch:")/outs/
mkdir -p $targetdir
cp %OUTPUTFILEDIR/%OUTPUTFILEBASE* $targetdir/

echo
echo "Job finished. SLURM details are:"
echo
qstat -f ${SLURM_JOB_ID}
echo
echo Finished at `date`
