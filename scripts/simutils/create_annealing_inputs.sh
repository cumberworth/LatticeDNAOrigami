#!/bin/bash

# Script to create input files for constant temperature runs from template

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

echo "Run number:"; read run
#run=
fields+=(RUN)
vars+=("$run")

echo "Reps:"; read reps
#reps=

echo "Max temp:"; read maxt
#maxt=
fields+=(MAXT)
vars=("$maxt")

echo "Min temp:"; read mint
#mint=
fields+=(MINT)
vars=("$mint")

echo "Temperature interval:"; read tint
#tint=
fields+=(TINT)
vars=("$tint")

echo "Steps per interval:"; read steps
#steps=
fields+=(STEPS)
vars+=("$steps")

echo "Max time per interval: (s)"; read maxdur
#maxdur=
fields+=(MAXDUR)
vars+=("$maxdur")

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

echo "Sceduler (pbs or slurm):"; read sched
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

echo "Movetype file:"; read movetypefile
#movetypefile=
fields+=(MOVETYPEFILE)
vars+=("$inpdir/$movetypefile")

echo "Restart file:"; read restartfile
#restartfile=
fields+=(RESTARTFILE)
if [ -z $restartfile ]
then
    vars+=("")
else
    vars+=("$inpdir/$restartfile")
fi

echo "Restart step:"; read restartstep
#restartstep=
fields+=(RESTARTSTEP)
vars+=("$restartstep")

echo "Default check/center/write/etc freq:"; read defaultint
#defaultint=1000000

echo "Centering freq:"; read centeringfreq
#centeringfreq=$defaultint
fields+=(CENTERINGFREQ)
vars+=("$centeringfreq")

echo "Constraint check freq:"; read concheckfreq
#concheckfreq=$defaultint
fields+=(CONCHECKFREQ)
vars+=("$concheckfreq")

echo "Logging freq:"; read loggingfreq
#loggingfreq=$defaultint
fields+=(LOGGINGFREQ)
vars+=("$loggingfreq")

echo "Config write freq (all formats):"; read configsfreq
#configsfreq=$defaultint
fields+=(CONFIGSFREQ)
vars+=("$configsfreq")

echo "Counts write freq:"; read countsfreq
#countsfreq=$defaultint
fields+=(COUNTSFREQ)
vars+=("$countsfreq")

echo "Tags of order parameters to output"; read ops
#ops=
fields+=(OPS)
vars+=("$ops")

echo "Order parameter write freq:"; read opfreq
#opfreq=$defaultint
fields+=(OPFREQ)
vars+=("$opfreq")

numfields=${#fields[@]}
fields+=(OUTPUTFILEBASE)
for ((rep=0; $rep<$reps; rep += 1))
do
    outputfilebase=${system}-${variant}_run-${run}_rep-${rep}
    vars[$numfields]=$outputfilebase

    sedcommand=""
    for i in ${!fields[@]}
    do
        sedcommand+="s:%${fields[i]}:${vars[i]}:g;"
    done
    sed "$sedcommand" serial_template_${sched}.sh > $inpdir/$outputfilebase.sh
    sed "$sedcommand" annealing_template.inp > $inpdir/$outputfilebase.inp
done
