# Functions for visualizing lattice DNA origami simulations
#
# origami must be set to the relevant mol id
# Note this assumes throughout all staples are two domains

# Color scheme and labels
set labels(scaffold_domain) "Unbound scaffold domain"
set colors(scaffold_domain) 0
# Dark-A
color change rgb 0 0.10588235294117647 0.6196078431372549 0.4666666666666667

set labels(staple_domain) "Unbound staple domain"
set colors(staple_domain) 1
# Dark2-B
color change rgb 1 0.8509803921568627 0.37254901960784315 0.00784313725490196

set labels(bound_domain) "Bound domain"
set colors(bound_domain) 2
# Dark2-A + Dark2-B
color change rgb 2 0.21568627450980393 0.3176470588235294 0.0784313725490196

set labels(misbound_domain) "Misbound domain"
set colors(misbound_domain) 3
# Complement of Dark2-A + Dark2-B
color change rgb 3 0.3176470588235294 0.0784313725490196 0.20784313725490197

set labels(scaffold_next_domain) "Scaffold next domain vector"
set colors(scaffold_next_domain) 4
# Dark2-F
color change rgb 4 0.9019607843137255 0.6705882352941176 0.00784313725490196

set labels(staple_next_domain) "Staple next domain vector"
set colors(staple_next_domain) 5
# Dark2-C
color change rgb 5 0.4588235294117647 0.4392156862745098 0.7019607843137254

set labels(scaffold_ore) "Scaffold orientation vector"
set colors(scaffold_ore) 6
# Dark2-D
color change rgb 6 0.9058823529411765 0.1607843137254902 0.5411764705882353

set labels(staple_ore) "Staple orientation vector"
set colors(staple_ore) 7
# Dark2-E
color change rgb 7 0.4 0.6509803921568628 0.11764705882352941

proc create_legend {} {
    # Create legend
    global colors
    global labels
    global legend
    if {[info exists legend]} {
        mol delrep all $legend
    }
    set legend [mol new]
    mol fix $legend
    set x 0.5
    set cur_y 1.3
    set incr 0.1
    foreach index {scaffold_domain staple_domain bound_domain misbound_domain \
            scaffold_next_domain staple_next_domain scaffold_ore staple_ore} {
        graphics $legend color $colors($index)
        graphics $legend text "$x $cur_y 0" $labels($index) size 1
        set cur_y [expr $cur_y - $incr]
    }
}

proc load_matrix_as_lists {filename} {
    # Load a matrix as a list of lists
    set f [open $filename r]
    set raw [read $f]
    close $f
    set lines [split $raw "\n"]
    set last_index [expr [llength $lines] - 1]
    set lines [lreplace $lines $last_index $last_index]

    return $lines
}

proc unpack_ores {raw_mat} {
    # Take a list of lists matrix and package the ores into a further list
    set ores {}
    foreach step $raw_mat {
        set step_ores {}
        for {set i 0} {$i != [llength $step]} {incr i 3} {
            set ore {}
            for {set j $i} {$j != [expr $i + 3]} {incr j} {
                lappend ore [lindex $step $j]
            }
            lappend step_ores $ore
        }
        lappend ores $step_ores
    }

    return $ores
}

proc calc_num_scaffold_domains {} {
    # Calculate number of scaffold domains
    global origami
    set numatoms [molinfo $origami get numatoms]
    set num_scaffolds 0
    for {set atom_i 0} {$atom_i != $numatoms} {incr atom_i} {
        set atom [atomselect $origami "index $atom_i"]
        if {[$atom get type] != "staple"} {
            incr num_scaffolds
        } else {
            break
        }
    }

    return $num_scaffolds
}

proc calc_num_staple_domains {} {
    # Calculate number of staples for given frame
    global origami
    global num_scaffold_domains
    global staplelength
    set numatoms [molinfo $origami get numatoms]

    return [expr $numatoms - $num_scaffold_domains]
}

proc create_domain_reps {} {
    global origami
    mol rep vdw 1.0 25
    mol rep vdw
    set numatoms [molinfo $origami get numatoms]
    for {set i 0} {$i != $numatoms} {incr i} {
        mol addrep $origami
        mol modselect $i $origami "index $i"
    }
}

proc update_colors {} {
    # Color the domains based on their binding states (and type)
    global colors
    global states
    global origami
    set numatoms [molinfo $origami get numatoms]
    set frame [molinfo $origami get frame]
    for {set i 0} {$i != $numatoms} {incr i} {
        set state [lindex [lindex $states $frame] $i]
        set atom [atomselect $origami "index $i" frame $frame]

        # Unbound domains
        if {$state == 1} {
            set type [$atom get type]
            if {$type == "scaffold"} {
                mol modcolor $i $origami ColorID $colors(scaffold_domain)
            } elseif {$type == "staple"} {
                mol modcolor $i $origami ColorID $colors(staple_domain)
            }

        # Fully bound domains
        } elseif {$state == 2} {
            mol modcolor $i $origami ColorID $colors(bound_domain)

        # Misbound domains
        } elseif {$state == 3} {
            mol modcolor $i $origami ColorID $colors(misbound_domain)
        }
    }
}

# Not sure how to pass arguments properly with callbacks, so just use globals
# It's what they do in the examples in the VMD docs

proc update_radii {} {
    # Set undefined domains' radii to 0 for current frame
    global states
    global origami
    global radius
    set frame [molinfo $origami get frame]
    set numatoms [molinfo $origami get numatoms]
    for {set atom_i 0} {$atom_i != $numatoms} {incr atom_i} {
        set atom [atomselect $origami "index $atom_i" frame $frame]
        set state [lindex [lindex $states $frame] $atom_i]
        if {$state == 0 || $state == -1} {
            $atom set radius 0
        } else {
            $atom set radius $radius
        }
    }
}

proc draw_3d_vector {origin vector color cylradius} {
    # Draw vector from origin
    global origami
    global arrowheadlength
    graphics $origami color $color
    set end [vecadd $origin $vector]
    set middle [vecadd $origin [vecsub $vector [vecscale $arrowheadlength [vecnorm $vector]]]]
    graphics $origami cylinder $origin $middle radius $cylradius resolution 25
    graphics $origami cone $middle $end radius 0.10 resolution 25
}

proc draw_next_domain_vectors {} {
    # Calculate and draw next domain vectors for current frame
    # Must clear previous first
    global colors
    global num_scaffold_domains
    global origami
    global states
    global radius
    global cylinderradius
    global staplelength
    global num_staple_domains
    set frame [molinfo $origami get frame]

    # Draw scaffold vectors
    set d1 [atomselect $origami "index 0" frame $frame]
    set d1_coors [lindex [$d1 get {x y z}] 0]
    for {set i 0} {$i != $num_scaffold_domains - 1} {incr i} {
        set d2_i [expr $i + 1]
        set d2 [atomselect $origami "index $d2_i" frame $frame]
        set d2_coors [lindex [$d2 get {x y z}] 0]
        set diff [vecsub $d2_coors $d1_coors]
        set vector [vecscale [expr 1 - $radius] $diff]
        set state1 [lindex [lindex $states $frame] $i]
        set state2 [lindex [lindex $states $frame] $d2_i]
        if {$state1 != 0 && $state2 != 0} {
            draw_3d_vector $d1_coors $vector $colors(scaffold_next_domain) \
                [vecscale 1.00 $cylinderradius]
        }
        set d1_coors $d2_coors
    }

    # Draw staple vectors
    for {set i 0} {$i != $num_staple_domains} {incr i $staplelength} {
        for {set j 0} {$j != [expr $staplelength - 1]} {incr j} {
            set d1_i [expr $i + $num_scaffold_domains + $j]
            set d2_i [expr $d1_i + 1]
            set d1 [atomselect $origami "index $d1_i" frame $frame]
            set d1_coors [lindex [$d1 get {x y z}] 0]
            set d2 [atomselect $origami "index $d2_i" frame $frame]
            set d2_coors [lindex [$d2 get {x y z}] 0]
            set diff [vecsub $d2_coors $d1_coors]
            set vector [vecscale $diff [expr 1 - $radius]]
            set state1 [lindex [lindex $states $frame] $d1_i]
            set state2 [lindex [lindex $states $frame] $d2_i]
            if {$state1 != 0 && $state1 != -1 && $state2 != 0 && $state2 != -1 && [veclength $vector] <= 1} {
                draw_3d_vector $d1_coors $vector $colors(staple_next_domain) \
                    [vecscale 0.99 $cylinderradius]
            }
        }
    }
}

proc draw_ore_vectors {} {
    # Draw orientation vectors for current frame
    # Must clear previous first
    global colors
    global ores
    global num_scaffold_domains
    global origami
    global states
    global cylinderradius
    global num_staple_domains
    set frame [molinfo $origami get frame]

    # Draw scaffold vectors
    for {set i 0} {$i != $num_scaffold_domains} {incr i} {
        set state [lindex [lindex $states $frame] $i]
        if {$state == 0} {
            continue
        }
        set d [atomselect $origami "index $i" frame $frame]
        set d_coors [lindex [$d get {x y z}] 0]
        set d_ore [lindex [lindex $ores $frame] $i]
        set vector [vecscale $d_ore 0.5]
        draw_3d_vector $d_coors $vector $colors(scaffold_ore) \
                [vecscale 0.98 $cylinderradius]
    }

    # Draw staple vectors
    set num_domains [expr $num_scaffold_domains + $num_staple_domains]
    for {set i $num_scaffold_domains} {$i != $num_domains} {incr i} {
        set state [lindex [lindex $states $frame] $i]
        if {$state == 0 || $state == -1} {
            continue
        }
        set d [atomselect $origami "index $i" frame $frame]
        set d_coors [lindex [$d get {x y z}] 0]
        set d_ore [lindex [lindex $ores $frame] $i]
        set vector [vecscale $d_ore 0.5]
        draw_3d_vector $d_coors $vector $colors(staple_ore) \
                [vecscale 0.97 $cylinderradius]
    }
}

proc update_frame {} {
    # Load new configuration and delete previous
    global origami
    global filebase
    global states
    global ores
    animate delete beg 0 $origami

    # Save visulation state
    foreach mol [molinfo list] {
        set viewpoints($mol) [molinfo $mol get {
            center_matrix rotate_matrix scale_matrix global_matrix}]
    }
    animate read vcf $filebase.vcf waitfor all $origami
    animate read vcf $filebase.vcf waitfor all $origami

    # Return to previous visulation state
    foreach mol [molinfo list] {
        molinfo $mol set {center_matrix rotate_matrix scale_matrix
            global_matrix} $viewpoints($mol)
    }

    set states [load_matrix_as_lists $filebase.states]
    lappend states [lindex $states 0]
    set ores_raw [load_matrix_as_lists $filebase.ores]
    set ores [unpack_ores $ores_raw]
    lappend ores [lindex $ores 0]

    graphics $origami delete all
    update_colors
    draw_next_domain_vectors
    draw_ore_vectors
    update_radii
}

proc update_frame_trace {args} {
    update_frame
}

proc draw_all_vectors {} {
    # Clear all previous graphics and draw vectors for current frame
    global origami
    graphics $origami delete all
    draw_next_domain_vectors
    draw_ore_vectors
}

proc update_colors_trace {args} {
    update_colors
}

proc update_radii_trace {args} {
    update_radii
}

proc draw_graphics_trace {args} {
    # Have to have all graphics stuff in one callback because of deletion
    draw_all_vectors
}
