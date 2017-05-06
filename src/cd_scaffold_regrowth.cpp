#include "cd_scaffold_regrowth.h"

using namespace CDScaffoldRegrowth;

CDScaffoldRegrowthMCMovetype::CDScaffoldRegrowthMCMovetype(
        OrigamiSystem& origami_system,
        RandomGens& random_gens,
        IdealRandomWalks& ideal_random_walks,
        InputParameters& params) :
        MCMovetype(
                origami_system,
                random_gens,
                ideal_random_walks,
                params) {

    // Setup staple simulation
    m_staple_sim = new ConstantTGCMCSimulation {m_origami_system, params};
    m_staple_sim->m_op_freq = 10000;
}

CDScaffoldRegrowthMCMovetype::~CDScaffoldRegrowthMCMovetype() {
    delete m_staple_sim;
}

bool CDScaffoldRegrowthMCMovetype::attempt_move() {
    bool accepted;

    // Select a scaffold segment to regrow
    vector<Domain*> scaffold {m_origami_system.get_chain(m_origami_system.c_scaffold)};
    int starting_d_i {m_random_gens.uniform_int(0, scaffold.size() - 1)};
    int dir {m_random_gens.uniform_int(0, 1)};
    int end_d_i;
    if (dir == 0) {
        dir = -1;
        end_d_i = 0;
    }
    else {
        end_d_i = scaffold.size() - 1;
    }

    vector<Domain*> domains {};
    for (int d_i {starting_d_i}; d_i != end_d_i + dir; d_i += dir) {
        Domain* cur_domain {scaffold[d_i]};
        domains.push_back(cur_domain);
    }

    for (size_t i {1}; i != domains.size(); i++) {
        Domain* domain {domains[i]};
        pair<int, int> key {domain->m_c, domain->m_d};
        m_prev_pos[key] = domain->m_pos;
        m_prev_ore[key] = domain->m_ore;
        m_modified_domains.push_back(key);
        m_delta_e += m_origami_system.unassign_domain(*domain);
    }

    // Grow chain
    for (size_t i {1}; i != domains.size(); i++) {
        Domain* domain {domains[i]};
        Domain* prev_domain {domains[i - 1]};
        VectorThree new_p {select_random_position(prev_domain->m_pos)};
        VectorThree new_o {select_random_orientation()};
        m_delta_e += m_origami_system.set_domain_config(*domain, new_p, new_o);
        if (m_origami_system.m_constraints_violated) {
            m_rejected = true;
            break;
        }
        else {
            pair<int, int> key {domain->m_c, domain->m_d};
            m_assigned_domains.push_back(key);
        }
    }

    if (m_rejected) {
        accepted = false;
        return accepted;
    }

    // Run staple simulation
    m_staple_sim->run();

    // Strip staples back off
    vector<vector<Domain*>> all_chains {m_origami_system.get_chains()};
    for (size_t i {1}; i != all_chains.size(); i++) {
        vector<Domain*> chain {all_chains[i]};
        int c_i {chain[0]->m_c};
        for (auto domain: chain) {
            m_origami_system.unassign_domain(*domain);
        }
        m_origami_system.delete_chain(c_i);
    }

    // Calculate averages and variances
    vector<double> enes {m_staple_sim->get_energies()};
    double enes_sum {std::accumulate(enes.begin() + m_burnin, enes.end(), 0.0)};
    m_enes_mean = enes_sum / (enes.size() - m_burnin);
    double enes_sq_sum {std::inner_product(enes.begin() + m_burnin, enes.end(), enes.begin(), 0.0)};
    double enes_var {enes_sq_sum / (enes.size() - m_burnin) - m_enes_mean * m_enes_mean};

    vector<int> staples {m_staple_sim->get_staples()};
    double staples_sum {std::accumulate(staples.begin() + m_burnin, staples.end(), 0.0)};
    m_staples_mean = staples_sum / (staples.size() - m_burnin);

    vector<int> domains_set {m_staple_sim->get_domains()};
    double domains_sum {std::accumulate(domains_set.begin() + m_burnin, domains_set.end(), 0.0)};
    m_domains_mean = domains_sum / (domains_set.size() - m_burnin);

    // Calculate CD factor
    double delta_staple_mean {m_enes_mean - m_prev_enes_mean};
    double cd_factor {exp(-(m_delta_e + delta_staple_mean + 2*enes_var))};
    accepted = test_acceptance(cd_factor);

    return accepted;
}
