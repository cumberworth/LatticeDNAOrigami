# Shitty pipe

set vmd_file_dir [lindex $argv 0]
set filebase [lindex $argv 1]
set staplelength [lindex $argv 2]

source $vmd_file_dir/liborigami.tcl

set origami [mol new $filebase.vsf]
set num_scaffold_domains [calc_num_scaffold_domains]
set num_staple_domains [calc_num_staple_domains]

color Display Background white
display shadows on
display ambientocclusion on
display aoambient 1.0
display aodirect 0.4
display resize 1000 1000
mol material AOChalky

set radius 0.20
set arrowheadlength 0.1
set cylinderradius 0.03

mol delrep 0 0
create_domain_reps
animate read vcf $filebase.vcf waitfor all $origami
animate read vcf $filebase.vcf waitfor all $origami
#create_legend
axes location off
display projection orthographic
mol top $origami

trace variable vmd_frame(0) w update_frame_trace

animate speed 0.1
animate forward
