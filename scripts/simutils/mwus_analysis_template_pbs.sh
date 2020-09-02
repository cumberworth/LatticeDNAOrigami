# Environment setup
cd ${PBS_O_WORKDIR}

system=%SYSTEM
var=%VARIANT
stapletypes=%STAPLETYPES
scaffolddomains=%SCAFFOLDDOMAINS
fullystackedpairs=%FULLYSTACKEDPAIRS
domainpairs="%DOMAINPAIRS"
temp=%TEMP
staplem=%STAPLEM
stackene=%STACKENE
tags="%TAGS"
tagpairs="%TAGPAIRS"
queue=%QUEUE
nodes=%NODES
procspernode=%PROCSPERNODE
prod_steps=%PRODSTEPS
walltime=%WALLTIME
sharedfiledir=%SHAREDFILEDIR
threads=%THREADS
reps=%REPS
skip=%SKIP
winfile=%WINFILE
biasfile=%BIASFILE

export PATH=~/bin:~/projects/origami/bin:~/.local/bin/:$PATH

#cp ${sharedfiledir}/${system}-${var}_run-${run}_rep-${rep}-0.vsf ${sharedfiledir}/${system}-${var}_run-${run}_rep-${rep}-${mult}.vsf
calc_mwus_numfullyboundstaples.py ${system}-${var} ${sharedfiledir} ${sharedfiledir} $winfile $biasfile $temp $staplem $reps $run
calc_mwus_scaffold_distances.py inps/${system}_unbound.json ${system}-${var} ${sharedfiledir} ${sharedfiledir} $winfile $biasfile $temp $staplem $reps $run --domain_pairs $domainpairs
calc_mwus_domain_occupancies.py ${system}-${var} ${sharedfiledir} ${sharedfiledir} $winfile $biasfile $temp $staplem $scaffolddomains $reps $run
perform_mwus_decorrelation.py ${system}-${var} ${sharedfiledir} analysis $winfile $biasfile $temp $staplem $stackene $skip $reps $run
perform_mwus_mbar.py ${system}-${var} ${sharedfiledir} analysis $winfile $biasfile $temp $staplem $stackene $reps $run --tags $tags
perform_mwus_slice_mbar.py ${system}-${var} ${sharedfiledir} analysis $winfile $biasfile $temp $staplem $stackene $reps $run numstaples $stapletypes $stapletypes $scaffolddomains --tags $tags
perform_mwus_slice_mbar.py ${system}-${var} ${sharedfiledir} analysis $winfile $biasfile $temp $staplem $stackene $reps $run numfullyboundstaples $stapletypes $stapletypes $scaffolddomains --tags $tags
perform_mwus_slice_mbar.py ${system}-${var} ${sharedfiledir} analysis $winfile $biasfile $temp $staplem $stackene $reps $run numfulldomains $stapletypes $stapletypes $scaffolddomains --tags $tags

# Run the plotting scripts
plot_mwus_ops_series.py ${sharedfiledir} plots $winfile ${system} ${var} $reps $run --assembled_values $stapletypes $scaffolddomains 0 $fullystackedpairs
plot_single-staple-curves.py analysis plots ${system} ${var} $stapletypes --xtag bias
plot_melting-lfes.py analysis plots/${system}-${var} numstaples --systems ${system} --varis ${var}
plot_melting-lfes.py analysis plots/${system}-${var} numfullyboundstaples --systems ${system} --varis ${var}
plot_melting-lfes.py analysis plots/${system}-${var} numfulldomains --systems ${system} --varis ${var}
plot_melting-lfes.py analysis plots/${system}-${var} numstackedpairs --systems ${system} --varis ${var}
plot_frequencies.py analysis plots ${system} ${var} ${stapletypes} inps/${system}_stapletype-map.dat
plot_domain-frequencies.py analysis plots ${system} ${var} ${scaffolddomains} inps/${system}_domaintype-map.dat

# Check if sampling converged
#sampling_converged=$(check_mwus_converged.py ${system}-${var} ${sharedfiledir} $staplem $stackene $fullystackedpairs $prod_steps --temps $smults)
#if [[ $sampling_converged == 1 ]]
#then
#    echo Sampling has converged
#else
#    echo Sampling has not converged
#fi

# Check bias converegence
#biases_converged=$(check_mwus_biases_converged.py $bias_error --old_biases $biases --new_biases $new_biases)
#new_biases=biases
#if [[ $biases_converged == 0 ]]
#then
#    echo Biases have not converged
#    biases=$new_biases
#    new_biases_i=$((biases_i + 1))
#    new_run=0
#else
#    echo Biases have converged
#    new_biases_i=${biases_i}
#    new_run=$((run + 1))
#fi

fields=()
vars=()
#fields+=(SMULTS)
#vars+=("$smults")

numfields=${#fields[@]}
fields+=(OUTPUTFILEBASE)
fields+=(RESTARTFILEBASE)
new_run=$((run + 1))
for ((rep=0; $rep<$reps; rep += 1))
do
    old_filebase=${system}-${var}_run-${run}_rep-${rep}
    new_filebase=${system}-${var}_run-${new_run}_rep-${rep}
    extract_mwus_restart_configs.py ${sharedfiledir}/${old_filebase} inps/${new_filebase} $winfile prod

    for file in $(ls ${sharedfiledir}| grep $old_filebase | grep biases)
    do
        newfile=inps/$(echo $file | sed "s/run-${run}/run-${new_run}/")
        cp -v ${sharedfiledir}/$file $newfile
    done

    vars[$numfields]=$new_filebase
    vars[$numfields + 1]=$new_filebase

    sedcommand=""
    for i in ${!fields[@]}
    do
        sedcommand+="s:%${fields[$i]}:${vars[$i]}:g;"
    done
    input_template=inps/${system}-${var}_template.inp
    sed "$sedcommand" $input_template > inps/${new_filebase}.inp
done

new_runfilebase=${system}-${var}_run-${new_run}

## Write current variables
#echo -e "smults_i=${new_smults_i}\nrun=${new_run}\nsmults=\"${smults}\"" >\
#    inps/${system}-${var}_current.sh
#
# Write restart script
echo "arrayid=\$(qsub inps/${system}-${var}_simulation.sh\
        -N ${new_runfilebase}\
        -q ${queue}\
        -l walltime=${walltime}:00:00\
        -l nodes=${nodes}:ppn=${procspernode}\
        -o outs/${new_runfilebase}.o\
        -e outs/${new_runfilebase}.e\
        -t 0-$((${reps} - 1))\
        -v run=${new_run})" >\
        inps/${system}-${var}_restart.sh

# Write analysis script
echo "qsub inps/${system}-${var}_analysis.sh\
        -N ${new_runfilebase}_analysis\
        -q ${queue}\
        -l walltime=3:00:00\
        -l nodes=1:ppn=1\
        -o outs/${new_runfilebase}_analysis.o\
        -e outs/${new_runfilebase}_analysis.e\
        -v run=${new_run}\
        -W depend=afterokarray:\${arrayid}" >>\
        inps/${system}-${var}_restart.sh

#if [[ $sampling_converged == 0 || $smults_converged == 0 ]]
#then
#    bash inps/${system}-${var}_restart.sh
#fi
