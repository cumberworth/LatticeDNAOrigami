#!/bin/bash

# Script to create input files for PTMC runs from template

# Note the variables being added to the array must be in quotes,
# otherwise it will add nothing

fields=()
vars=()

#echo "Input directory:"; read inpdir
inpdir=inps
fields+=(INPDIR)
vars+=("$inpdir")

echo "System:"; read system
#system=
fields+=(SYSTEM)
vars+=("$system")

#echo "State (unbound or assembled):"; read state
state=unbound
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
#run=0
fields+=(RUN)
vars+=("$run")

echo "Reps:"; read reps
#reps=3

echo "Number of replicas:"; read numreps
#numreps = 
fields+=(NUMREPS)
vars+=("$numreps")

echo "Exchange frequency:"; read exchangefreq
#exchangefreq=
fields+=(EXCHANGEFREQ)
vars+=("$exchangefreq")

echo "Temps (seperate with spaces:"; read temps
#temps="330 335 340 345 350 355 360 365 370"
fields+=(TEMPS)
vars+=("$temps")

#echo "Chemical potential multipliers:"; read chempotmults
chempotmults=""
for ((i=0; $i<$numreps; i+=1))
do
    chempotmults+=" 1"
done
fields+=(CHEMPOTMULTS)
vars+=("$chempotmults")

#echo "Bias potential multipliers:"; read biasmults
biasmults=""
for ((i=0; $i<$numreps; i+=1))
do
    biasmults+=" 1"
done
fields+=(BIASMULTS)
vars+=("$biasmults")

echo "Stacking energy multipliers:"; read stackingmults
#stackingmults="1.0 0.5 0.0"
fields+=(STACKINGMULTS)
vars+=("$stackingmults")

#echo "Queue:"; read queue
queue=l1
fields+=(QUEUE)
vars+=("$queue")

echo "Walltime hours:"; read walltime
#walltime=30
fields+=(WALLTIME)
vars+=("$walltime")

#echo "Number of nodes:"; read nodes
nodes=1
fields+=(NODES)
vars+=("$nodes")

#echo "Number of procs per node:"; read procs
procs=$numreps
fields+=(PROCS)
vars+=("$procs")

#echo "Node list:"; read nodenames
#procs=$nodenames
#fields+=(NODENAMES)
#vars+=("$nodenames")

#echo "Output file directory:"; read outputfiledir
outputfiledir=/scratch/amc226
fields+=(OUTPUTFILEDIR)
vars+=("$outputfiledir")

echo "Sceduler (pbs or slurm):"; read sched
#sched=pbs

echo "Hybridization potential:"; read hybridpot
#hybridpot=NearestNeighbour
fields+=(HYBRIDPOT)
vars+=("$hybridpot")

#echo "Domain type:"; read domaintype
domaintype=ThreeQuarterTurn
fields+=(DOMAINTYPE)
vars+=("$domaintype")

#echo "Binding potential:"; read bindpot
bindpot=ConKinkLinearFlexible
fields+=(BINDPOT)
vars+=("$bindpot")

#echo "Misbinding potential:"; read misbindpot
misbindpot=Opposing
fields+=(MISBINDPOT)
vars+=("$misbindpot")

#echo "Stacking potential:"; read stackpot
stackpot=Constant
fields+=(STACKPOT)
vars+=("$stackpot")

echo "Staple concentration (mol/L):"; read staplem
#staplem=1e-7
fields+=(STAPLEM)
vars+=("$staplem")

echo "Cation concentration (mol/L):"; read cationm
#cationm=0.5
fields+=(CATIONM)
vars+=("$cationm")

echo "Maximum number of total staples:"; read maxstaples
#maxstaples=4
fields+=(MAXSTAPLES)
vars+=("$maxstaples")

echo "Maximum number of staples of given type:"; read maxtypestaples
#maxtypestaples=$maxstaples
fields+=(MAXTYPESTAPLES)
vars+=("$maxtypestaples")

echo "Maximum number of domains per staple:"; read maxsizestaple
#maxstaplesize=$maxstaplesize
fields+=(MAXSIZESTAPLE)
vars+=("$maxsizestaple")

#echo "Domain level biases present?"; read domainbiases
domainbiases=false
fields+=(DOMAINBIASES)
vars+=("$domainbiases")

echo "Binding enthalpy"; read bindh
#bindh=-60059.27167992321
#bindh=0
fields+=(BINDH)
vars+=("$bindh")

echo "Binding entropy"; read binds
#binds=-164.0135846334504
#binds=0
fields+=(BINDS)
vars+=("$binds")

echo "Misbinding enthalpy"; read misbindh
#misbindh=-8989.354123837773
#misbindh=0
fields+=(MISBINDH)
vars+=("$misbindh")

echo "Misbinding entropy"; read misbinds
#misbinds=-26.857736488186973
#misbinds=0
fields+=(MISBINDS)
vars+=("$misbinds")

echo "Stacking energy"; read stackene
#stackene=
fields+=(STACKENE)
vars+=("$stackene")

echo "Order parameter file:"; read opfile
#opfile=ops_default.json
fields+=(OPFILE)
if [ -z $opfile ]
then
    vars+=("")
else
    vars+=("$inpdir/$opfile")
fi

#echo "Bias function file:"; read biasfile
biasfile=
fields+=(BIASFILE)
if [ -z $biasfile ]
then
    vars+=("")
else
    vars+=("$inpdir/$biasfile")
fi

#echo "Restart step:"; read restartstep
restartstep=0
fields+=(RESTARTSTEP)
vars+=("$restartstep")

#echo "Restart postfix:"; read restartpostfix
restartpostfix=.trj.restart
fields+=(RESTARTPOSTFIX)
vars+=("$restartpostfix")

echo "Default check/center/write/etc freq:"; read defaultint
#defaultint=100000

#echo "Centering freq:"; read centeringfreq
centeringfreq=$defaultint
fields+=(CENTERINGFREQ)
vars+=("$centeringfreq")

#echo "Constraint check freq:"; read concheckfreq
concheckfreq=1000000
fields+=(CONCHECKFREQ)
vars+=("$concheckfreq")

#echo "Swaps:"; read swaps
swaps=1000000000000
fields+=(SWAPS)
vars+=("$swaps")

echo "Maximum swap duration (s):"; read maxdur
#maxdur=100000
fields+=(MAXDUR)
vars+=("$maxdur")

#echo "Logging freq:"; read loggingfreq
loggingfreq=$defaultint
fields+=(LOGGINGFREQ)
vars+=("$loggingfreq")

#echo "Config write freq (all formats):"; read configsfreq
configsfreq=$defaultint
fields+=(CONFIGSFREQ)
vars+=("$configsfreq")

#echo "Counts write freq:"; read countsfreq
countsfreq=$defaultint
fields+=(COUNTSFREQ)
vars+=("$countsfreq")

echo "Tags of order parameters to output"; read ops
#ops="numfulldomains nummisdomains numstackedpairs numstaples numlinearhelices numstackedjuncts"
fields+=(OPS)
vars+=("$ops")

#echo "Order parameter write freq:"; read opfreq
opfreq=$defaultint
fields+=(OPFREQ)
vars+=("$opfreq")

#echo "Timing write freq:"; read timefreq
timefreq=$defaultint
fields+=(TIMEFREQ)
vars+=("$timefreq")

#echo "Energy write freq:"; read energyfreq
energyfreq=$defaultint
fields+=(ENERGYFREQ)
vars+=("$energyfreq")

#movetypefile=${inpdir}/moveset_default.json
movetypefile=${inpdir}/moveset_${system}-${variant}.json
fields+=(MOVETYPEFILE)
vars+=("$movetypefile")

numfields=${#fields[@]}
fields+=(OUTPUTFILEBASE)
fields+=(REP)
fields+=(RESTARTFILEBASE)
fields+=(RESTARTSWAPFILE)
for ((rep=0; $rep<$reps; rep += 1))
do
    outputfilebase=${system}-${variant}_run-${run}_rep-${rep}
    restartfilebase=${outputfilebase}
    restartswapfile=${inpdir}/${restartfilebase}.swp.restart
    vars[$numfields]=$outputfilebase
    vars[$numfields + 1]=$rep
    vars[$numfields + 2]=$restartfilebase
    vars[$numfields + 3]=$restartswapfile

    sedcommand=""
    for i in ${!fields[@]}
    do
        sedcommand+="s:%${fields[i]}:${vars[i]}:g;"
    done

    if (( $run > 0 ))
    then
        input_template=2d_ptmc_restart_template.inp
    else
        input_template=2d_ptmc_template.inp
    fi

    sed "$sedcommand" ptmc_template_${sched}.sh > $inpdir/$outputfilebase.sh
    sed "$sedcommand" $input_template > $inpdir/$outputfilebase.inp
done
