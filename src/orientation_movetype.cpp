// orientation_movetype.cpp

#include "orientation_movetype.h"

namespace movetypes {

    using utility::Occupancy;

	OrientationRotationMCMovetype::OrientationRotationMCMovetype(
                OrigamiSystem& origami_system,
                RandomGens& random_gens,
                IdealRandomWalks& ideal_random_walks,
                vector<OrigamiOutputFile*> config_files,
                string label,
                SystemOrderParams& ops,
                SystemBiases& biases,
                InputParameters& params) :
        MCMovetype(origami_system, random_gens, ideal_random_walks,
                config_files, label, ops, biases, params) {
	}

    void OrientationRotationMCMovetype::write_log_summary(ostream*) {
    }

    bool OrientationRotationMCMovetype::internal_attempt_move() {
        bool accepted {false};
        
        // Select random chain, domain, and orientation
        Domain* domain {select_random_domain()};
        VectorThree o_new {select_random_orientation()};

        if (domain->m_state == Occupancy::bound or domain->m_state == Occupancy::misbound) {
            VectorThree o_old {domain->m_ore};
            Domain* bound_domain {domain->m_bound_domain};
            m_origami_system.unassign_domain(*bound_domain);
            m_origami_system.set_domain_orientation(*domain, o_new);
            VectorThree pos {domain->m_pos};
            m_origami_system.set_domain_config(*bound_domain, pos, -o_new);
            if (m_origami_system.m_constraints_violated) {
                accepted = false;
                m_origami_system.unassign_domain(*bound_domain);
                m_origami_system.set_domain_orientation(*domain, o_old);
                m_origami_system.set_checked_domain_config(*bound_domain, pos, -o_old);
            }
            else {
                accepted = true;
            }
        }
        else {
            m_origami_system.set_domain_orientation(*domain, o_new);
            accepted = true;
        }
        write_config();

        return accepted;
    }

    void OrientationRotationMCMovetype::add_tracker(bool) {
    }
}
