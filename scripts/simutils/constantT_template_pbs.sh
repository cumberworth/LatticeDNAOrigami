# Name of job
#PBS -N %OUTPUTFILEBASE

# Queue to use
#PBS -q %QUEUE

# Nodes and procs
#PBS -l nodes=1:ppn=1

# Walltime limit (hours:mins:secs)
#PBS -l walltime=%WALLTIME:00:00

# Standard error and out files
#PBS -o %OUTPUTFILEBASE.o
#PBS -e %OUTPUTFILEBASE.e

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
export PATH=~/bin/$PATH

# Main job
latticeDNAOrigami -i %OUTPUTFILEBASE.inp > %OUTPUTFILEDIR/%OUTPUTFILEBASE.out

# Copy results to slowscratch mirror
cp %OUTPUTFILEDIR/%OUTPUTFILEBASE.* /sharedscratch/amc226/projects/origami/calculations/17-03-29_us-validation

echo
echo "Job finished. PBS details are:"
echo
qstat -f ${PBS_JOBID}
echo
echo Finished at `date`
