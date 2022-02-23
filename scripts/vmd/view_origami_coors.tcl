# Requires the variables libdir, filebase, system, and staplelength have been set

source $libdir/liborigami.tcl

set origami [mol new $system.vsf]
set ores_raw [load_matrix_as_lists $filebase.ores]
set ores [unpack_ores $ores_raw]
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
mol addfile $filebase.vcf type vcf waitfor all
#create_legend
axes location off
display projection orthographic
mol top $origami
set states [load_matrix_as_lists $filebase.states]
trace variable vmd_frame(0) w update_colors_trace
trace variable vmd_frame(0) w draw_graphics_trace
trace variable vmd_frame(0) w update_radii_trace

animate goto start
