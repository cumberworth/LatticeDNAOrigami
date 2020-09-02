# Environment setup
cd ${PBS_O_WORKDIR}
module unload gcc
module load gcc/6.2.0
nodelist=$(sort $PBS_NODEFILE | uniq)

system=%SYSTEM
var=%VARIANT
wins=%WINS
outputfiledir=%OUTPUTFILEDIR
sharedfiledir=%SHAREDFILEDIR
rep=${PBS_ARRAYID}

filebase=${system}-${var}_run-${run}_rep-${rep}

echo "Starting job $PBS_JOBID"
echo "PBS assigned me these nodes:"
echo $nodelist
echo

export LD_LIBRARY_PATH=~/lib:$LD_LIBRARY_PATH
export PATH=~/bin:~/projects/origami/bin:$PATH

# Create output directories on each node
for node in $nodelist
do
    ssh $node mkdir -p $outputfiledir
done

# Run simulation
inputfile=inps/${filebase}.inp
outputfile=${outputfiledir}/${filebase}.out
mpirun -n $wins latticeDNAOrigami -i $inputfile > $outputfile

# Create central output directory for collection
mkdir -p ${sharedfiledir}

# Copy results central directory
for node in $nodelist
do
    ssh $node rsync -avz ${outputfiledir}/${filebase}* ${sharedfiledir}
done

echo
echo "Job finished. PBS details are"
echo
qstat -f ${PBS_JOBID}
echo
echo Finished at `date`
