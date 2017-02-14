// order_params.h

#ifndef ORDER_PARAMS_H
#define ORDER_PARAMS_H

#include <vector>

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
    };

    class DistOrderParam: public OrderParam {
        // Distance between two given domains
        public:
            DistOrderParam(Domain& domain_1, Domain& domain_2);
            int calc_param();

        private:
            Domain& m_domain_1;
            Domain& m_domain_2;
    };

    class BiasFunction {
        // Functions of order parameters that calculate an energy bias
        public:
            virtual ~BiasFunction() {};
            virtual double calc_bias() = 0;
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
        private:
            OrderParam& m_order_param;
            int m_min_param;
            int m_max_param;
            double m_max_bias;
            double m_slope;
    };

    class SystemBias {
        // All bias functions for the system
        public:
            SystemBias(InputParameters& params, OrigamiSystem& origami);
            ~SystemBias();
            double calc_bias();
            void update_bias_mult(double bias_mult);
        private:
            double m_bias_mult {1};
            vector<BiasFunction*> m_biases {};
    };
}

#endif // ORDER_PARAMS_H
