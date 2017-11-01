// order_params.h

#ifndef ORDER_PARAMS_H
#define ORDER_PARAMS_H

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "domain.h"
#include "hash.h"
#include "origami_system.h"
#include "parser.h"
#include "utility.h"

namespace orderParams {

    using std::pair;
    using std::reference_wrapper;
    using std::string;
    using std::unique_ptr;
    using std::unordered_map;
    using std::vector;

    using domainContainer::Domain;
    using origami::OrigamiSystem;
    using parser::InputParameters;
    using utility::VectorThree;
    using utility::Occupancy;

    class OrderParam {
        // System parameters calculated from current config
        public:
            virtual ~OrderParam() {}
            virtual int calc_param() = 0;
            virtual int check_param(
                    Domain& domain,
                    VectorThree new_pos,
                    VectorThree new_ore,
                    Occupancy new_state) = 0;
            int get_param();
            int get_checked_param();
            string get_label();
            bool defined();

        protected:
            string m_label {};
            int m_param;
            int m_checked_param;
            bool m_defined;
    };

    class DistOrderParam: public OrderParam {
        // Distance between two given domains
        public:
            DistOrderParam(Domain& domain_1, Domain& domain_2, string label);

            int calc_param() override final;
            int check_param(
                    Domain& domain,
                    VectorThree new_pos,
                    VectorThree,
                    Occupancy) override final;

        private:
            Domain& m_domain_1;
            Domain& m_domain_2;
    };

    class AdjacentSiteOrderParam: public OrderParam {
        // If given domain pair occupies adjacent lattice sites
        public:
            AdjacentSiteOrderParam(
                    Domain& domain_1,
                    Domain& domain_2,
                    string label);

            int calc_param() override final;
            int check_param(
                    Domain& domain,
                    VectorThree new_pos,
                    VectorThree,
                    Occupancy) override final;

        private:
            Domain& m_domain_1;
            Domain& m_domain_2;
    };

    class SumOrderParam: public OrderParam {
        public:
            SumOrderParam(
                    vector<reference_wrapper<OrderParam>> ops,
                    string label);

            int calc_param() override final;
            int check_param(
                    Domain&,
                    VectorThree,
                    VectorThree,
                    Occupancy) override final;

        private:
            vector<reference_wrapper<OrderParam>> m_ops {};
    };

    class NumStaplesOrderParam: public OrderParam {
        public:
            NumStaplesOrderParam(OrigamiSystem& origami, string label);
            int calc_param() override final;
            int check_param(
                    Domain& domain,
                    VectorThree new_pos,
                    VectorThree,
                    Occupancy) override final;

        private:
            OrigamiSystem& m_origami;
    };

    class NumBoundDomainPairsOrderParam: public OrderParam {
        public:
            NumBoundDomainPairsOrderParam(OrigamiSystem& origami, string label);
            int calc_param() override final;
            int check_param(
                    Domain& domain,
                    VectorThree new_pos,
                    VectorThree,
                    Occupancy) override final;

        private:
            OrigamiSystem& m_origami;
    };

    class NumMisBoundDomainPairsOrderParam: public OrderParam {
        public:
            NumMisBoundDomainPairsOrderParam(OrigamiSystem& origami, string label);
            int calc_param() override final;
            int check_param(
                    Domain& domain,
                    VectorThree new_pos,
                    VectorThree,
                    Occupancy) override final;

        private:
            OrigamiSystem& m_origami;
    };

    class SystemOrderParams {
        public:
            SystemOrderParams(InputParameters& params, OrigamiSystem& origami);
            SystemOrderParams(const SystemOrderParams&) = delete;
            SystemOrderParams& operator=(const SystemOrderParams&) = delete;

            OrderParam& get_order_param(string tag);
            vector<pair<int, int>> get_dependent_domains(string tag);

            void update_all_params();
            void update_move_params();
            void update_one_domain(Domain& domain);
            void check_one_domain(
                    Domain& domain,
                    VectorThree pos,
                    VectorThree ore,
                    Occupancy state);

        private:
            void setup_ops(
                    string ops_filename,
                    vector<pair<int, int>> keys);

            OrigamiSystem& m_origami;

            // Vectors of each type of order param
            vector<vector<unique_ptr<OrderParam>>> m_level_to_ops {};
            unordered_map<string, reference_wrapper<OrderParam>> m_tag_to_op;

            // Store distance order parameters that are dependent on given domain
            unordered_map<pair<int, int>, vector<vector<reference_wrapper<OrderParam>>>>
                    m_domain_update_ops {};
            vector<vector<reference_wrapper<OrderParam>>> m_move_update_ops {};
            unordered_map<string, vector<pair<int, int>>> m_tag_to_domains {};

            int m_levels;
    };
}

#endif // ORDER_PARAMS_H
