#!/bin/bash

fields=()
vars=()

echo "System:"
read system
#system=
fields+=(SYSTEM)
vars+=($system)

echo "State (unbound or assembled):"
read state
#state=
fields+=(STATE)
vars+=($state)

echo "Variant:"
read variant
#variant=
fields+=(VARIANT)
vars+=($variant)

echo "Run number:"
read run
#run=
fields+=(RUN)
vars+=($run)

echo "Reps:"
read reps
#reps=

echo "Temps (seperate with spaces:"
read temps
#temps
temps=($temps)

echo "Queue:"
read queue
#queue=
fields+=(QUEUE)
vars+=($queue)

echo "Walltime hours:"
read walltime
#walltime=
fields+=(WALLTIME)
vars+=($walltime)

echo "Output file directory:"
read outputfiledir
#outputfiledir=/scratch/amc226
fields+=(OUTPUTFILEDIR)
vars+=($outputfiledir)

echo "Scheduler (pbs or slurm):"
read sched
#sched=

echo "Staple concentration (mol/L):"
read staplem
#staplem=1e-6
fields+=(STAPLEM)
vars+=($staplem)

echo "Cation concentration (mol/L):"
read cationm
#cationm=1
fields+=(CATIONM)
vars+=($cationm)

echo "Maximum number of total staples:"
read maxstaples
#maxstaples=
fields+=(MAXSTAPLES)
vars+=($maxstaples)

echo "Maximum number of staples of given type:"
read maxstaplestype
#maxstaplestype=
fields+=(MAXSTAPLESTYPE)
vars+=($maxstaplestype)

echo "Domain level biases present?"
read domainbiases
#domainbiases=
fields+=(DOMAINBIASES)
vars+=($domainbiases)

echo "Order parameter file:"
read opfile
#opfile=
fields+=(OPFILE)
vars+=($opfile)

echo "Bias function file:"
read biasfile
#biasfile=
fields+=(BIASFILE)
vars+=($biasfile)

echo "Movetype file:"
read movetypefile
#movetypefile=
fields+=(MOVETYPEFILE)
vars+=($movetypefile)

echo "Restart file:"
read restartfile
#restartfile=
fields+=(RESTARTFILE)
vars+=($restartfile)

echo "Restart step:"
read restartstep
#restartstep=
fields+=(RESTARTSTEP)
vars+=($restartstep)

echo "Default check/center/write/etc freq:"
read defaultint
#defaultint=1000000

echo "Centering freq:"
read centeringfreq
#centeringfreq=$defaultint
fields+=(CENTERINGFREQ)
vars+=($centeringfreq)

echo "Constraint check freq:"
read concheckfreq
#concheckfreq=$defaultint
fields+=(CONCHECKFREQ)
vars+=($concheckfreq)

echo "Steps:"
read steps
#steps=
fields+=(STEPS)
vars+=($steps)

echo "Logging freq:"
read loggingfreq
#loggingfreq=$defaultint
fields+=(LOGGINGFREQ)
vars+=($loggingfreq)

echo "Config write freq (all formats):"
read configsfreq
#configsfreq=$defaultint
fields+=(CONFIGSFREQ)
vars+=($configsfreq)

echo "Counts write freq:"
read countsfreq
#countsfreq=$defaultint
fields+=(COUNTSFREQ)
vars+=($countsfreq)

echo "Tags of order parameters to output"
read ops
#ops=
fields+=(OPS)
vars+=($ops)

echo "Order parameter write freq:"
read opsfreq
#opsfreq=$defaultint
fields+=(OPSFREQ)
vars+=($opsfreq)

outputfilebase=%SYSTEM-%VARIANT_run-%RUN_rep-%REP-%TEMP
fields+=(OUTPUTFILEBASE)
vars+=($outputfilebase)

for ((rep=0; $rep<$reps; rep += 1))
do
    for temp in {$temps[@]}
    do
        sedcommand=""
        for i in ${!fields[@]}
        do
            sedcommand+="s/%${fields[i]}/${vars[i]}/g;"
        done
        sed "$sedcommand" constantT_template_${sched}.sh > $outputfilebase.sh
        sed "$sedcommand" constantT_template.inp > $outputfilebase.inp
    done
done
