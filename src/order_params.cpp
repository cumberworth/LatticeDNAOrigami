// order_params.cpp

#include <numeric>
#include <utility>

#include "files.h"
#include "order_params.h"

namespace orderParams {

    using std::abs;
    using std::cout;
    using std::set;

    using files::OrigamiOrderParamsFile;

    string OrderParam::get_label() {
        return m_label;
    }

    int OrderParam::get_param() {
        return m_param;
    }

    int OrderParam::get_checked_param() {
        return m_checked_param;
    }

    bool OrderParam::defined() {
        return m_defined;
    }

    DistOrderParam::DistOrderParam(
            Domain& domain_1,
            Domain& domain_2,
            string label) :
            m_domain_1 {domain_1},
            m_domain_2 {domain_2} {

        m_label = label;
        calc_param();
    }

    int DistOrderParam::calc_param() {
        if (m_domain_1.m_state != Occupancy::unassigned and
                m_domain_2.m_state != Occupancy::unassigned) {

            m_defined = true;
            VectorThree diff_vec {m_domain_2.m_pos - m_domain_1.m_pos};
            m_param = diff_vec.abssum();
            m_checked_param = m_param;
        }
        else {
            m_defined = false;
        }

        return m_param;
    }

    int DistOrderParam::check_param(Domain& domain, VectorThree new_pos, VectorThree,
            Occupancy state) {
        Domain* unmodded_domain {&m_domain_1};
        if (domain.m_d == m_domain_1.m_d) {
            unmodded_domain = &m_domain_2;
        }
        if (unmodded_domain->m_state != Occupancy::unassigned and
                state != Occupancy::unassigned) {
            m_checked_param = (new_pos - unmodded_domain->m_pos).abssum();
            // Will be defined for configurations where that domain is unassigned
            // Consider using a check defined variable instead
            m_defined = true;
        }
        else {
            m_defined = false;
        }

        return m_checked_param;
    }

    AdjacentSiteOrderParam::AdjacentSiteOrderParam(
            Domain& domain_1,
            Domain& domain_2,
            string label) :
            m_domain_1 {domain_1},
            m_domain_2 {domain_2} {

        m_label = label;
        calc_param();
    }

    int AdjacentSiteOrderParam::calc_param() {
        if (m_domain_1.m_state != Occupancy::unassigned and
                m_domain_2.m_state != Occupancy::unassigned) {

            m_defined = true;
            VectorThree diff_vec {m_domain_2.m_pos - m_domain_1.m_pos};
            auto dist = diff_vec.abssum();
            if (dist == 1) {
                m_param = 1;
            }
            else {
                m_param = 0;
            }
            m_checked_param = m_param;
        }
        else {
            m_defined = false;
        }

        return m_param;
    }

    int AdjacentSiteOrderParam::check_param(Domain& domain, VectorThree new_pos, VectorThree,
            Occupancy state) {
        Domain* unmodded_domain {&m_domain_1};
        if (domain.m_d == m_domain_1.m_d) {
            unmodded_domain = &m_domain_2;
        }
        if (unmodded_domain->m_state != Occupancy::unassigned and
                state != Occupancy::unassigned) {
            // Will be defined for configurations where that domain is unassigned
            // Consider using a check defined variable instead
            auto dist = (new_pos - unmodded_domain->m_pos).abssum();
            if (dist == 1) {
            m_checked_param = 1;
            }
            else {
            m_checked_param = 0;
            }
            m_checked_param = 
            m_defined = true;
        }
        else {
            m_defined = false;
        }

        return m_checked_param;
    }

    SumOrderParam::SumOrderParam(vector<reference_wrapper<OrderParam>> ops, string label) :
            m_ops {ops} {
        m_label = label;
        calc_param();
    }

    int SumOrderParam::calc_param() {
        int op_sum {0};
        m_defined = true;
        for (auto param: m_ops) {

            // All constituent parameters must be defined
            if (param.get().defined()) {
                op_sum += param.get().get_param();
            }
            else {
                m_defined = false;
                op_sum = m_param;
                break;
            }
        }
        m_param = op_sum;
        m_checked_param = op_sum;

        return m_param;
    }

    int SumOrderParam::check_param(Domain&, VectorThree, VectorThree,
            Occupancy) {
        int op_sum {0};
        m_defined = true;
        for (auto param: m_ops) {

            // All constituent parameters must be defined
            if (param.get().defined()) {
                op_sum += param.get().get_checked_param();
            }
            else {
                m_defined = false;
                op_sum = m_param;
                break;
            }
        }
        m_checked_param = op_sum;

        return op_sum;
    }

    NumStaplesOrderParam::NumStaplesOrderParam(
            OrigamiSystem& origami,
            string label) :
            m_origami {origami} {

        m_label = label;
        calc_param();
        m_defined = true;
    }

    int NumStaplesOrderParam::calc_param() {
        m_param = m_origami.num_staples();

        return m_param;
    }

    int NumStaplesOrderParam::check_param(Domain&, VectorThree, VectorThree,
            Occupancy) {

        m_checked_param = m_origami.num_staples();

        return m_checked_param;
    }

    NumStaplesTypeOrderParam::NumStaplesTypeOrderParam(
            OrigamiSystem& origami,
            int c_ident,
            string label) :
            m_origami {origami} {

        m_c_ident = c_ident;
        m_label = label;
        calc_param();
        m_defined = true;
    }

    int NumStaplesTypeOrderParam::calc_param() {
        m_param = m_origami.staples_of_ident(m_c_ident).size();

        return m_param;
    }

    int NumStaplesTypeOrderParam::check_param(Domain&, VectorThree, VectorThree,
            Occupancy) {

        cout << "CHECK STACKING NOT IMPLEMENTED\n";
        return 0;
    }

    StapleTypeFullyBoundOrderParam::StapleTypeFullyBoundOrderParam(
            OrigamiSystem& origami,
            int c_ident,
            string label) :
            m_origami {origami} {

        m_c_ident = c_ident;
        m_label = label;
        calc_param();
        m_defined = true;
    }

    int StapleTypeFullyBoundOrderParam::calc_param() {
        m_param = 0;
        for (auto staple_i: m_origami.staples_of_ident(m_c_ident)) {
            bool fully_bound {true};
            for (auto domain: m_origami.get_chain(staple_i)) {
                if (domain->m_state != Occupancy::bound) {
                    fully_bound = false;
                    break;
                }
            }
            if (fully_bound) {
                m_param = 1;
                break;
            }
        }

        return m_param;
    }

    int StapleTypeFullyBoundOrderParam::check_param(Domain&, VectorThree, VectorThree,
            Occupancy) {

        cout << "CHECK STACKING NOT IMPLEMENTED\n";
        return 0;
    }

    NumBoundDomainPairsOrderParam::NumBoundDomainPairsOrderParam(
            OrigamiSystem& origami,
            string label) :
            m_origami {origami} {

        m_label = label;
        calc_param();
        m_defined = true;
    }

    int NumBoundDomainPairsOrderParam::calc_param() {
        m_param = m_origami.num_fully_bound_domain_pairs();

        return m_param;
    }

    int NumBoundDomainPairsOrderParam::check_param(Domain&, VectorThree,
            VectorThree, Occupancy state) {
        
        if (state != Occupancy::unassigned) {
            m_defined = true;
            int num_domain_pairs = m_origami.num_fully_bound_domain_pairs();
            if (state == Occupancy::bound) {
                m_checked_param = num_domain_pairs + 1;
            }
            else {
                m_checked_param = num_domain_pairs;
            }
        }

        return m_checked_param;
    }

    NumMisboundDomainPairsOrderParam::NumMisboundDomainPairsOrderParam(
            OrigamiSystem& origami,
            string label) :
            m_origami {origami} {

        m_label = label;
        calc_param();
        m_defined = true;
    }

    int NumMisboundDomainPairsOrderParam::calc_param() {
        m_param = m_origami.num_misbound_domain_pairs();

        return m_param;
    }

    int NumMisboundDomainPairsOrderParam::check_param(Domain&, VectorThree,
            VectorThree, Occupancy state) {
        
        if (state != Occupancy::unassigned) {
            m_defined = true;
            int num_domain_pairs = m_origami.num_misbound_domain_pairs();
            if (state == Occupancy::misbound) {
                m_checked_param = num_domain_pairs + 1;
            }
            else {
                m_checked_param = num_domain_pairs;
            }
        }

        return m_checked_param;
    }

    NumStackedPairsOrderParam::NumStackedPairsOrderParam(
            OrigamiSystem& origami,
            string label) :
            m_origami {origami} {

        m_label = label;
        calc_param();
        m_defined = true;
    }

    int NumStackedPairsOrderParam::calc_param() {
        m_param = m_origami.num_stacked_domain_pairs();

        return m_param;
    }

    int NumStackedPairsOrderParam::check_param(Domain&, VectorThree,
            VectorThree, Occupancy) {
        
/*        if (state != Occupancy::unassigned) {
            m_defined = true;
            m_checked_param = m_origami.num_stacked_domain_pairs();
            m_checked_param += m_origami.check_stacking(domain);
        }

        return m_checked_param;
        */
        cout << "CHECK STACKING NOT IMPLEMENTED\n";
        return 0;
    }

    NumLinearHelicesOrderParam::NumLinearHelicesOrderParam(
            OrigamiSystem& origami,
            string label) :
            m_origami {origami} {

        m_label = label;
        calc_param();
        m_defined = true;
    }

    int NumLinearHelicesOrderParam::calc_param() {
        m_param = m_origami.num_linear_helix_trips();

        return m_param;
    }

    int NumLinearHelicesOrderParam::check_param(Domain&, VectorThree,
            VectorThree, Occupancy) {
        
/*        if (state != Occupancy::unassigned) {
            m_defined = true;
            m_checked_param = m_origami.num_stacked_domain_pairs();
            m_checked_param += m_origami.check_stacking(domain);
        }

        return m_checked_param;
        */
        cout << "CHECK LINEAR HELICES NOT IMPLEMENTED\n";
        return 0;
    }

    NumStackedJunctsOrderParam::NumStackedJunctsOrderParam(
            OrigamiSystem& origami,
            string label) :
            m_origami {origami} {

        m_label = label;
        calc_param();
        m_defined = true;
    }

    int NumStackedJunctsOrderParam::calc_param() {
        m_param = m_origami.num_stacked_junct_quads();

        return m_param;
    }

    int NumStackedJunctsOrderParam::check_param(Domain&, VectorThree,
            VectorThree, Occupancy) {
        
/*        if (state != Occupancy::unassigned) {
            m_defined = true;
            m_checked_param = m_origami.num_stacked_domain_pairs();
            m_checked_param += m_origami.check_stacking(domain);
        }

        return m_checked_param;
        */
        cout << "CHECK STACKED JUNCTIONS NOT IMPLEMENTED\n";
        return 0;
    }

    SystemOrderParams::SystemOrderParams(InputParameters& params,
            OrigamiSystem& origami) :
            m_origami {origami} {

        // Setup dependency table
        vector<pair<int, int>> keys {};
        for (auto chain: origami.get_chains()) {
            for (auto domain: chain) {
                pair<int, int> key {domain->m_c, domain->m_d};
                m_domain_update_ops[key] = {};
                keys.push_back(key);
            }
        }

        if (params.m_ops_filename != "") {
            setup_ops(params.m_ops_filename, keys);
        }
    }

    void SystemOrderParams::setup_ops(
            string ops_filename,
            vector<pair<int, int>> keys) {

        OrigamiOrderParamsFile ops_file {ops_filename};
        vector<vector<string>> level_to_types {ops_file.get_types_by_level()};
        vector<vector<string>> level_to_tags {ops_file.get_tags_by_level()};
        vector<vector<string>> level_to_labels {ops_file.get_labels_by_level()};
        for (size_t i {0}; i != level_to_types.size(); i++) {
            m_level_to_ops.push_back({});
            m_move_update_ops.push_back({});
            m_levels++;
            for (auto key: keys) {
                m_domain_update_ops[key].push_back({});
            }

            bool update_per_domain {false};
            for (size_t j {0}; j != level_to_types[i].size(); j++) {
                OrderParam* op;
            	set<pair<int, int>> d_domains {};
                string type {level_to_types[i][j]};
                string tag {level_to_tags[i][j]};
                string label {level_to_labels[i][j]};
                m_tag_to_domains[tag] = {};
                if (type == "Dist" or type == "AdjacentSite") {
                    int c1i {ops_file.get_int_option(i, j, "chain1")};
                    int d1i {ops_file.get_int_option(i, j, "domain1")};
                    int c2i {ops_file.get_int_option(i, j, "chain2")};
                    int d2i {ops_file.get_int_option(i, j, "domain2")};
					d_domains.insert({c1i, d1i});
					d_domains.insert({c2i, d2i});
                    Domain& d1 {*m_origami.get_domain(c1i, d1i)};
                    Domain& d2 {*m_origami.get_domain(c2i, d2i)};
                    if (type == "Dist") {
                        op = new DistOrderParam {d1, d2, label};
                    }
                    else if (type == "AdjacentSite") {
                        op = new AdjacentSiteOrderParam {d1, d2, label};
                    }

                    update_per_domain = ops_file.get_bool_option(i, j,
                            "update_per_domain");
                }
                else if (type == "Sum") {
                    vector<string> op_tags_to_sum {
                        ops_file.get_vector_string_option(i, j, "ops")};
                    vector<reference_wrapper<OrderParam>> d_ops {};
                    update_per_domain = true;
                    for (auto d_tag: op_tags_to_sum) {
                        vector<pair<int, int>> vd_domains {m_tag_to_domains[d_tag]};
                        if (vd_domains.size() == 0) {
                            update_per_domain = false;
                        }
                        for (auto d_domain: vd_domains) {
                            d_domains.insert(d_domain);
                        }
                        reference_wrapper<OrderParam> d_op {m_tag_to_op.at(d_tag)};
                        d_ops.push_back(d_op);
                    }
                    op = new SumOrderParam {d_ops, label};
                }
                else if (type == "NumStaples")  {
                    op = new NumStaplesOrderParam {m_origami, label};
                }
                else if (type == "NumStaplesType")  {
                    int staple {ops_file.get_int_option(i, j, "staple")};
                    op = new NumStaplesTypeOrderParam {m_origami, staple,
                            label};
                }
                else if (type == "StapleTypeFullyBound")  {
                    int staple {ops_file.get_int_option(i, j, "staple")};
                    op = new StapleTypeFullyBoundOrderParam {m_origami, staple,
                            label};
                }
                else if (type == "NumBoundDomainPairs") {
                    op = new NumBoundDomainPairsOrderParam {
                            m_origami, label};
                }
                else if (type == "NumMisboundDomainPairs") {
                    op = new NumMisboundDomainPairsOrderParam {
                            m_origami, label};
                }
                else if (type == "NumStackedPairs") {
                    op = new NumStackedPairsOrderParam {
                            m_origami, label};
                }
                else if (type == "NumLinearHelices") {
                    op = new NumLinearHelicesOrderParam {
                            m_origami, label};
                }
                else if (type == "NumStackedJuncts") {
                    op = new NumStackedJunctsOrderParam {
                            m_origami, label};
                }
                else {
                    cout << "Order parameter type does not exist";
                    throw utility::SimulationMisuse {};
                }
                m_level_to_ops[i].emplace_back(op);
				OrderParam& op_ref {*m_level_to_ops[i].back()};
                m_tag_to_op.emplace(tag, op_ref);
                if (not update_per_domain) {
                    m_move_update_ops[i].emplace_back(op_ref);
                }
                else {
                    for (auto d_domain: d_domains) {
                        m_domain_update_ops[d_domain][i].emplace_back(op_ref);
                        m_tag_to_domains[tag].push_back(d_domain);
                    }
                }
            }
        }
    }

    OrderParam& SystemOrderParams::get_order_param(string tag) {
        return m_tag_to_op.at(tag);
    }

    vector<pair<int, int>> SystemOrderParams::get_dependent_domains(
            string tag) {
        return m_tag_to_domains[tag];
    }

    void SystemOrderParams::update_all_params() {
        for (auto& level: m_level_to_ops) {
            for (auto& op: level) {
                op->calc_param();
            }
        }
    }

    void SystemOrderParams::update_move_params() {
        for (auto level: m_move_update_ops) {
            for (auto op: level) {
                op.get().calc_param();
            }
        }
    }

    void SystemOrderParams::update_one_domain(Domain& domain) {

        pair<int, int> key {domain.m_c, domain.m_d};
        for (auto level: m_domain_update_ops[key]) {
            for (auto op: level) {
                op.get().calc_param();
            }
        }
    }

    void SystemOrderParams::check_one_domain(Domain& domain, VectorThree pos,
            VectorThree ore, Occupancy state) {

        pair<int, int> key {domain.m_c, domain.m_d};
        for (auto level: m_domain_update_ops[key]) {
            for (auto op: level) {
                op.get().check_param(domain, pos, ore, state);
            }
        }
    }
}
