# Name of job
#PBS -N %OUTPUTFILEBASE

# Queue to use
#PBS -q %QUEUE

# Nodes and procs
#PBS -l nodes=%NODENAMES:ppn=1

# Walltime limit (hours:mins:secs)
#PBS -l walltime=%WALLTIME:00:00

# Standard error and out files
#PBS -o outs/%OUTPUTFILEBASE.o
#PBS -e outs/%OUTPUTFILEBASE.e

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
latticeDNAOrigami -i %INPDIR/%OUTPUTFILEBASE.inp > %OUTPUTFILEDIR/%OUTPUTFILEBASE.out

# Copy results to slowscratch mirror
targetdir=$(pwd | sed "s:home:sharedscratch:")/outs/
mkdir -p $targetdir
cp %OUTPUTFILEDIR/%OUTPUTFILEBASE.* $targetdir/

echo
echo "Job finished. PBS details are:"
echo
qstat -f ${PBS_JOBID}
echo
echo Finished at `date`
