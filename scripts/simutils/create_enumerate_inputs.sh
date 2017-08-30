#!/bin/bash

# Script to create input files for enumeration runs from template

# Note the variables being added to the array must be in quotes,
# otherwise it will add nothing

fields=()
vars=()

echo "Input directory:"; read inpdir
#inpdir=
fields+=(INPDIR)
vars+=("$inpdir")

echo "System:"; read system
#system=
fields+=(SYSTEM)
vars+=("$system")

echo "State (unbound or assembled):"; read state
#state=
fields+=(STATE)
vars+=("$state")

sysfile=$inpdir/${system}_${state}.json
fields+=(SYSFILE)
vars+=("$sysfile")

echo "Variant:"; read variant
#variant=
fields+=(VARIANT)
vars+=("$variant")

echo "Temps (seperate with spaces:"; read temps
#temps
temps=($temps)

echo "Queue:"; read queue
#queue=
fields+=(QUEUE)
vars+=("$queue")

echo "Walltime hours:"; read walltime
#walltime=
fields+=(WALLTIME)
vars+=("$walltime")

echo "Output file directory:"; read outputfiledir
#outputfiledir=/scratch/amc226
fields+=(OUTPUTFILEDIR)
vars+=("$outputfiledir")

echo "Sceduler (pbs or slurm):"; read ; ched
#sched=

echo "Staple concentration (mol/L):"; read staplem
#staplem=1e-6
fields+=(STAPLEM)
vars+=("$staplem")

echo "Cation concentration (mol/L):"; read cationm
#cationm=1
fields+=(CATIONM)
vars+=("$cationm")

echo "Maximum number of total staples:"; read maxstaples
#maxstaples=
fields+=(MAXSTAPLES)
vars+=("$maxstaples")

echo "Maximum number of staples of given type:"; read maxtypestaples
#maxtypestaples=
fields+=(MAXTYPESTAPLES)
vars+=("$maxtypestaples")

echo "Domain level biases present?"; read domainbiases
#domainbiases=
fields+=(DOMAINBIASES)
vars+=("$domainbiases")

echo "Order parameter file:"; read opfile
#opfile=
fields+=(OPFILE)
if [ -z $opfile ]
then
    vars+=("")
else
    vars+=("$inpdir/$opfile")
fi

echo "Bias function file:"; read biasfile
#biasfile=
fields+=(BIASFILE)
if [ -z $biasfile ]
then
    vars+=("")
else
    vars+=("$inpdir/$biasfile")
fi

echo "Tags of order parameters to output"; read ops
#ops=
fields+=(OPS)
vars+=("$ops")

numfields=${#fields[@]}
fields+=(OUTPUTFILEBASE)
fields+=(TEMP)
for temp in ${temps[@]}
do
    outputfilebase=${system}-${variant}_temp-${temp}
    vars[$numfields]=$outputfilebase
    vars[$(numfields + 1)]=$temp

    sedcommand=""
    for i in ${!fields[@]}
    do
        sedcommand+="s:%${fields[i]}:${vars[i]}:g;"
    done
    sed "$sedcommand" serial_template_${sched}.sh > $inpdir/$outputfilebase.sh
    sed "$sedcommand" enumerate_template.inp > $inpdir/$outputfilebase.inp
done
