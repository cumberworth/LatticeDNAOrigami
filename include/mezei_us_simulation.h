// mezei_us_simulation.h

#ifndef MEZEI_US_SIMULATION_H
#define MEZEI_US_SIMULATION_H

#include "ceres/ceres.h"

#include "simulation.h"
#include "us_simulation.h"

namespace MezeiUS {

    class MezeiUSGCMCSimulation: public USGCMCSimulation {
        // Adaptive US method described in Mezei 1987
        public:
            MezeiUSGCMCSimulation(
                    OrigamiSystem& origami,
                    InputParameters& params);
            ~MezeiUSGCMCSimulation() {}

        private:
            ofstream m_solver_file;

            // Variable names set to be consistent with mezei1987
            ArrayOfSets m_s_n {}; // grid points visited at iteration n
            GridOfIntArrays m_f_n {}; // number of visits at each grid point at each iteration
            GridOfIntArrays m_F_n {}; // total visits at each grid point over all iterations
            GridOfFloatArrays m_r_n {}; // relative contribution of grid points across iterations
            GridOfFloatArrays m_w_n {}; // relative contribution of grid points for each iteration

            GridOfFloatArrays m_p_n {}; // locally normalized weight of grid point for each iteration
            GridFloats m_P_n {}; // unormalized weight of grid point for all iterations
            vector<double> m_N {}; // normalization of each each iteration

            bool iteration_equilibrium_step();
            void update_bias(int n);
            void update_grids(int n);
            void output_summary(int n);
            void estimate_normalizations(int n);
            double estimate_initial_normalization(int n);
            void estimate_final_normalizations();
    };

    class RDevSquareSumCostFunction: public ceres::CostFunction {
        public: 
            RDevSquareSumCostFunction(int num_residuals, vector<int>& parameter_block_sizes);
            ~RDevSquareSumCostFunction() {}
            bool Evaluate(double const* const* parameters,
                        double* residuals,
                        double** jacobians) const;
    };

    struct RDevSquareSum {
        // For solving the relative deviation square sum with ceres
        template <typename T> 
        bool operator() (
                T const* const* params,
                T* residual) const;
    };
}

#endif // MEZEI_US_SIMULATION_H
