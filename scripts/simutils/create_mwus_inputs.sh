#!/bin/bash

# Script to create input files for PTMC runs from template

# Note the variables being added to the array must be in quotes,
# otherwise it will add nothing

source $1

fields=()
vars=()

#echo "Input directory:"; read inpdir
fields+=(INPDIR)
vars+=("$inpdir")

#echo "System:"; read system
fields+=(SYSTEM)
vars+=("$system")

fields+=(DOMAINPAIRS)
vars+=("$domainpairs")

fields+=(STAPLETYPES)
vars+=("$stapletypes")

fields+=(FULLYSTACKEDPAIRS)
vars+=("$fullystackedpairs")

fields+=(SCAFFOLDDOMAINS)
vars+=("$scaffolddomains")

#echo "State (unbound or assembled):"; read state
fields+=(STATE)
vars+=("$state")

sysfile=$inpdir/${system}_${state}.json
fields+=(SYSFILE)
vars+=("$sysfile")

#echo "Variant:"; read variant
fields+=(VARIANT)
vars+=("$var")

#echo "Reps:"; read reps
fields+=(REPS)
vars+=("$reps")

#echo "Temp:"; read temp
fields+=(TEMP)
vars+=("$temp")

#echo "Windows file:"; read winfile
fields+=(WINFILE)
vars+=("$winfile")

#echo "Max bias change:"; read maxdbias
fields+=(MAXDBIAS)
vars+=("$maxdbias")

#echo "Equilibrium steps:"; read esteps
fields+=(ESTEPS)
vars+=("$esteps")

#echo "Max time for equilibration (s):"; read emaxdur
fields+=(EMAXDUR)
vars+=("$emaxdur")

#echo "Steps per iteration:"; read isteps
fields+=(ISTEPS)
vars+=("$isteps")

#echo "Max time per iteration (s):"; read imaxdur
fields+=(IMAXDUR)
vars+=("$imaxdur")

#echo "Production steps:"; read psteps
fields+=(PSTEPS)
vars+=("$psteps")

#echo "Max time for production (s):"; read pmaxdur
fields+=(PMAXDUR)
vars+=("$pmaxdur")

#echo "Iterations:"; read iters
fields+=(ITERS)
vars+=("$iters")

#echo "Number of windows:"; read wins
fields+=(WINS)
vars+=("$wins")

#echo "Queue:"; read queue
fields+=(QUEUE)
vars+=("$queue")

#echo "Walltime hours:"; read walltime
fields+=(WALLTIME)
vars+=("$walltime")

#echo "Number of nodes:"; read nodes
fields+=(NODES)
vars+=("$nodes")

#echo "Number of procs per node:"; read procspernode
fields+=(PROCSPERNODE)
vars+=("$procspernode")

#echo "Grid bias tag:"; read gridtag
fields+=(GRIDTAG)
vars+=("$gridtag")

#echo "Node list:"; read nodenames
#fields+=(NODENAMES)
#vars+=("$nodenames")

#echo "Output file directory:"; read outputfiledir
fields+=(OUTPUTFILEDIR)
vars+=("$outputfiledir")

fields+=(SHAREDFILEDIR)
vars+=("$sharedfiledir")

#echo "Sceduler (pbs or slurm):"; read sched

#echo "Hybridization potential:"; read hybridpot
fields+=(HYBRIDPOT)
vars+=("$hybridpot")

#echo "Domain type:"; read domaintype
fields+=(DOMAINTYPE)
vars+=("$domaintype")

#echo "Binding potential:"; read bindpot
fields+=(BINDPOT)
vars+=("$bindpot")

#echo "Misbinding potential:"; read misbindpot
fields+=(MISBINDPOT)
vars+=("$misbindpot")

#echo "Stacking potential:"; read stackpot
fields+=(STACKPOT)
vars+=("$stackpot")

#echo "Staple concentration (mol/L):"; read staplem
fields+=(STAPLEM)
vars+=("$staplem")

#echo "Cation concentration (mol/L):"; read cationm
fields+=(CATIONM)
vars+=("$cationm")

fields+=(TEMPERROR)
vars+=("$temp_error")

#echo "Maximum number of total staples:"; read maxstaples
fields+=(MAXSTAPLES)
vars+=("$maxstaples")

#echo "Maximum number of staples of given type:"; read maxtypestaples
fields+=(MAXTYPESTAPLES)
vars+=("$maxtypestaples")

#echo "Maximum number of domains per staple:"; read maxsizestaple
fields+=(MAXSIZESTAPLE)
vars+=("$maxsizestaple")

#echo "Domain level biases present?"; read domainbiases
fields+=(DOMAINBIASES)
vars+=("$domainbiases")

#echo "Binding enthalpy for uniform potential"; read bindh
fields+=(BINDH)
vars+=("$bindh")

#echo "Binding entropy for uniform potential"; read binds
fields+=(BINDS)
vars+=("$binds")

#echo "Misbinding enthalpy for uniform potential"; read misbindh
fields+=(MISBINDH)
vars+=("$misbindh")

#echo "Misbinding entropy for uniform potential"; read misbinds
fields+=(MISBINDS)
vars+=("$misbinds")

#echo "Stacking energy"; read stackene
fields+=(STACKENE)
vars+=("$stackene")

#echo "Order parameter file:"; read opfile
fields+=(OPFILE)
vars+=("$opfile")

#echo "Bias function file:"; read biasfile
fields+=(BIASFILE)
vars+=("$biasfile")

#echo "Restart step:"; read restartstep
fields+=(RESTARTSTEP)
vars+=("$restartstep")

#echo "Restart postfix:"; read restartpostfix
fields+=(RESTARTPOSTFIX)
vars+=("$restartpostfix")

#echo "Default check/center/write/etc freq:"; read defaultint

fields+=(SKIP)
vars+=("$skip")

#echo "Centering freq:"; read centeringfreq
fields+=(CENTERINGFREQ)
vars+=("$centeringfreq")

#echo "Constraint check freq:"; read concheckfreq
fields+=(CONCHECKFREQ)
vars+=("$concheckfreq")

#echo "Logging freq:"; read loggingfreq
fields+=(LOGGINGFREQ)
vars+=("$loggingfreq")

#echo "Config write freq (all formats):"; read configsfreq
fields+=(CONFIGSFREQ)
vars+=("$configsfreq")

#echo "Counts write freq:"; read countsfreq
fields+=(COUNTSFREQ)
vars+=("$countsfreq")

#echo "Tags of order parameters to output"; read tags
fields+=(TAGS)
vars+=("$tags")

fields+=(TAGPAIRS)
vars+=("$tagpairs")

#echo "Order parameter write freq:"; read opfreq
fields+=(OPFREQ)
vars+=("$opfreq")

#echo "Timing write freq:"; read timefreq
fields+=(TIMEFREQ)
vars+=("$timefreq")

#echo "Energy write freq:"; read energyfreq
fields+=(ENERGYFREQ)
vars+=("$energyfreq")

fields+=(MOVETYPEFILE)
vars+=("$movetypefile")

sedcommand=""
for i in ${!fields[@]}
do
    sedcommand+="s:%${fields[i]}:${vars[i]}:g;"
done

sed "$sedcommand" mwus_simulation_template_${sched}.sh > $inpdir/${system}-${var}_simulation.sh
sed "$sedcommand" mwus_analysis_template_${sched}.sh > $inpdir/${system}-${var}_analysis.sh
sed "$sedcommand" mwus_template.inp > $inpdir/${system}-${var}_template.inp

# Create inputs for first run
run=0

numfields=${#fields[@]}
fields+=(OUTPUTFILEBASE)
fields+=(REP)
fields+=(RESTARTFILEBASE)
for ((rep=0; $rep<$reps; rep += 1))
do
    select_starting_configs.py inps/${system}_unbound.json $winfile $biasfile ${seedrunprerep}rep-${rep}${seedrunposrep} inps/${system}-${var}_run-${run}_rep-${rep} .trj.restart
    filebase=${system}-${var}_run-${run}_rep-${rep}

    vars[$numfields]=$filebase
    vars[$numfields + 1]=$rep
    vars[$numfields + 2]=$filebase

    sedcommand=""
    for i in ${!fields[@]}
    do
        sedcommand+="s:%${fields[i]}:${vars[i]}:g;"
    done

    sed "$sedcommand" mwus_init_template.inp > $inpdir/${filebase}.inp
done
runfilebase=${system}-${var}_run-${run}
echo "arrayid=\$(qsub inps/${system}-${var}_simulation.sh\
    -N ${runfilebase}\
    -q ${queue}\
    -l walltime=${walltime}:00:00\
    -l nodes=${nodes}:ppn=${procspernode}\
    -o outs/${runfilebase}.o\
    -e outs/${runfilebase}.e\
    -t 0-$((${reps} - 1))\
    -v run=${run})" >\
    inps/${system}-${var}_init.sh
echo "qsub inps/${system}-${var}_analysis.sh\
    -N ${runfilebase}_analysis\
    -l walltime=3:00:00\
    -l nodes=1:ppn=1\
    -o outs/${runfilebase}_analysis.o\
    -e outs/${runfilebase}_analysis.e\
    -v run=${run}\
    -W depend=afterokarray:\${arrayid}" >>\
    inps/${system}-${var}_init.sh 

echo -e "run=${run}" >\
    inps/${system}-${var}_current.sh
