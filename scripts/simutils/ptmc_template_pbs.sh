# Name of job
#PBS -N %OUTPUTFILEBASE

# Queue to use
#PBS -q %QUEUE

# Walltime limit (hours:mins:secs)
#PBS -l walltime=%WALLTIME:00:00

# Nodes and procs
#PBS -l nodes=%NODES:ppn=%PROCS

# Standard error and out files
#PBS -o outs/%OUTPUTFILEBASE.o
#PBS -e outs/%OUTPUTFILEBASE.e

# Environment setup
cd ${PBS_O_WORKDIR}
module unload gcc
module load gcc/6.2.0

echo "Starting job $PBS_JOBID"
echo
echo "PBS assigned me this node:"
cat $PBS_NODEFILE
echo

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
echo "Job finished. PBS details are:"
echo
qstat -f ${PBS_JOBID}
echo
echo Finished at `date`
