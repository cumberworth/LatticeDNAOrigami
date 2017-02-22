// order_params.h

#ifndef ORDER_PARAMS_H
#define ORDER_PARAMS_H

#include <vector>

#include "hash.h"
#include "parser.h"
#include "domain.h"
#include "origami_system.h"

using std::vector;

using namespace Parser;
using namespace DomainContainer;
using namespace Origami;

namespace OrderParams {

    class OrderParam {
        // System parameters calculated from current config
        public:
            virtual ~OrderParam() {};
            virtual int calc_param() = 0;
            virtual int calc_param(
                    Domain& domain,
                    VectorThree new_pos,
                    VectorThree new_ore,
                    Occupancy new_state) = 0;
            virtual int get_param() = 0;
            //virtual void update_param() = 0;
            //virtual int get_param() = 0;

            // Check whether order parameters dependent on given domain
            virtual bool dependent_on(Domain& domain) = 0;

            // Check whether order parameter currently defined
            virtual bool defined() = 0;
    };

    class DistOrderParam: public OrderParam {
        // Distance between two given domains
        public:
            DistOrderParam(Domain& domain_1, Domain& domain_2);
            ~DistOrderParam() {}

            int calc_param();
            int calc_param(
                    Domain& domain,
                    VectorThree new_pos,
                    VectorThree,
                    Occupancy);
            int get_param();
            bool dependent_on(Domain& domain);
            bool defined();

        private:
            Domain& m_domain_1;
            Domain& m_domain_2;
            int m_param;
            bool m_defined;
    };

    class BiasFunction {
        // Functions of order parameters that calculate an energy bias
        public:
            virtual ~BiasFunction() {};
            virtual double calc_bias() = 0;
            virtual double update_bias() = 0;
            virtual bool dependent_on(OrderParam& order_param) = 0;
            virtual double get_bias() = 0;
    };

    class LinearStepBiasFunction: public BiasFunction {
        // Constant below a min and above a max value, and linear between
        public:
            LinearStepBiasFunction(
                    OrderParam& order_param,
                    int min_param,
                    int max_param,
                    double max_bias);
            ~LinearStepBiasFunction();
            double calc_bias();
            double update_bias();
            bool dependent_on(OrderParam& order_param);
            double get_bias();
        private:
            OrderParam& m_order_param;
            int m_min_param;
            int m_max_param;
            double m_max_bias;
            double m_slope;
            double m_bias;
    };

    class SystemBias {
        // All bias functions for the system
        public:
            SystemBias(InputParameters& params, OrigamiSystem& origami);
            ~SystemBias();

            // Recalculate everything
            double calc_bias();

            // Bias can be tuned with a multiplier
            void update_bias_mult(double bias_mult);

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
            double check_one_domain(
                    Domain& domain,
                    VectorThree new_pos,
                    VectorThree new_ore,
                    Occupancy new_state);
        private:
            OrigamiSystem& m_origami;
            double m_bias_mult {1};
            vector<BiasFunction*> m_bias_fs {};
            vector<OrderParam*> m_order_params {};
            unordered_map<pair<int, int>, vector<OrderParam*>> m_domain_to_order_params {};
            unordered_map<pair<int, int>, vector<BiasFunction*>> m_domain_to_bias_fs{};

            void setup_distance_bias(InputParameters& params);
            void add_order_param_dependency(Domain& domain, OrderParam* order_param);
            void add_bias_f_dependency(Domain& domain, BiasFunction* bias_f);
            void remove_domain_dependencies(Domain& domain);
            void add_chain(vector<Domain*> chain);
            void remove_chain(vector<Domain*> chain);
            vector<OrderParam*> get_dependent_order_params(Domain& domain);
            vector<BiasFunction*> get_dependent_bias_fs(Domain& domain);
    };
}

#endif // ORDER_PARAMS_H
