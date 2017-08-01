// bias_functions.h

#ifndef BIAS_FUNCTIONS_H
#define BIAS_FUNCTIONS_H

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "order_params.h"
#include "domain.h"
#include "hash.h"
#include "origami_system.h"
#include "parser.h"
#include "utility.h"

namespace biasFunctions {

    using std::pair;
    using std::reference_wrapper;
    using std::string;
	using std::unique_ptr;
    using std::unordered_map;
    using std::vector;

    using parser::InputParameters;
    using domainContainer::Domain;
    using orderParams::OrderParam;
    using orderParams::SystemOrderParams;
    using origami::OrigamiSystem;
    using utility::VectorThree;
    using utility::Occupancy;

    class BiasFunction {
        // Functions of order parameters that calculate an energy bias
        public:
            virtual ~BiasFunction() {};
            virtual double update_bias() = 0;
            virtual double check_bias() = 0;
            virtual double get_bias() = 0;

        protected:
            double m_bias {0};
    };

    class LinearStepBiasFunction: public BiasFunction {
        // Constant below a min and above a max value, and linear between
        public:
            LinearStepBiasFunction(
                    OrderParam& order_param,
                    int min_param,
                    int max_param,
                    double max_bias);
            ~LinearStepBiasFunction() {};
            double update_bias();
            double check_bias();
            double get_bias();
        private:
            OrderParam& m_order_param;
            int m_min_param;
            int m_max_param;
            double m_max_bias;
            double m_slope;

            double calc_bias(int param);
    };

    class SquareWellBiasFunction: public BiasFunction {
        // a between and including min and max values of given parameter, b otherwise
        public:
            SquareWellBiasFunction(
                    OrderParam& order_param,
                    int min_param,
                    int max_param,
                    double well_bias,
                    double oustide_bias);
            ~SquareWellBiasFunction() {};
                double update_bias();
            double check_bias();
            double get_bias();

        private:
            OrderParam& m_order_param;
            int m_min_param;
            int m_max_param;
            double m_well_bias;
            double m_outside_bias;

            double calc_bias(int param);
    };

    using GridPoint = vector<int>;
    class GridBiasFunction: public BiasFunction {
        public:
            GridBiasFunction(vector<reference_wrapper<OrderParam>> ops);
            ~GridBiasFunction() {};

            int get_dim();
            vector<int> get_point();

            void replace_biases(unordered_map<GridPoint, double> bias_grid);

            double update_bias();
            double check_bias();
            double get_bias();
            double calc_bias(vector<int> params);

        private:
            vector<reference_wrapper<OrderParam>> m_ops;
            unordered_map<vector<int>, double> m_bias_grid {};
            double m_off_grid_bias {0};
    };

    // Maybe have two classes here; one for bias functions that should be
    // updated per domain change, and those that are only valid at a system
    // level
    class SystemBiases {
        // All bias functions for the system
        public:
            SystemBiases(
                    OrigamiSystem& origami,
                    SystemOrderParams& system_order_params,
                    InputParameters& params);

            double get_bias();

            // Bias can be tuned with a multiplier
            void update_bias_mult(double bias_mult);

            // Recalculate everything
            double calc_total();
            double calc_move();

            // Calculate change in bias from the setting of one domain
            /* At this point, the bias object requires this to be called
                every time the origami system object changes a domain (assigning
                or unassing); if not done, the change in bias will be wrong,
                as it will recalculate the changed biases with the current
                configuration of the origami system, but will have an out-of-date
                value for the last bias. The burden of calling each
                time can be on either the movetypes or the origami object.
                Consider having a way of checking whether other domains have
                changed and throwing an exception.
            */
            double calc_one_domain(Domain& domain);

            // Check change in bias from changing one domain but don't update total bias
            double check_one_domain(Domain& domain);

            GridBiasFunction& get_grid_bias(string tag);

        private:
            //void add_chain(vector<Domain*> chain);
            //void remove_chain(vector<Domain*> chain);

            OrigamiSystem& m_origami;
            SystemOrderParams& m_system_order_params;
            double m_bias_mult {1};

            // Vectors of each type of order param
            vector<vector<unique_ptr<BiasFunction>>> m_level_to_biases {};
            unordered_map<string, reference_wrapper<BiasFunction>> m_tag_to_biases;

            // Store distance order parameters that are dependent on given domain
            unordered_map<pair<int, int>, vector<vector<reference_wrapper<BiasFunction>>>>
                    m_domain_update_biases {};
            vector<vector<reference_wrapper<BiasFunction>>> m_move_update_biases {};

            double m_total_bias {0};
            int m_levels;
    };
}

#endif // BIAS_FUNCTIONS_H
