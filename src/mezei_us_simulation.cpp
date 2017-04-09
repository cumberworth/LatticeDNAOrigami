// mezei_simulation.cpp

#include "mezei_us_simulation.h"
#include "json/json.h"
#include "ceres/ceres.h"
#include "glog/logging.h"

using ceres::DynamicAutoDiffCostFunction;
using ceres::CostFunction;
using ceres::Problem;
using ceres::Solver;
using ceres::Solve;
using ceres::LossFunction;
using ceres::ScaledLoss;

using namespace Simulation;
using namespace US;
using namespace MezeiUS;

MezeiUSGCMCSimulation::MezeiUSGCMCSimulation(OrigamiSystem& origami,
        InputParameters& params):
        USGCMCSimulation {origami, params} {

    // Initalize log file
    m_solver_file.open(params.m_output_filebase + ".solver");
}

void MezeiUSGCMCSimulation::update_grids(int n) {
    // Update tracking stuff
    m_s_n.push_back(m_s_i);

    for (auto point: m_new_points) {
        vector<int> f_k_n(n);
        f_k_n.push_back(m_f_i[point]);
        m_f_n[point] = f_k_n;

        m_F_n[point] = f_k_n;

        vector<double> r_k_n(n);
        r_k_n.push_back(1);
        m_r_n[point] = r_k_n;

        vector<double> w_k_n(n);
        w_k_n.push_back(static_cast<double>(m_f_i[point]) / m_steps);
        m_w_n[point] = w_k_n;
    }
    for (auto point: m_old_points) {
        m_f_n[point].push_back(m_f_i[point]);

        int F_k_n {std::accumulate(m_f_n[point].begin(), m_f_n[point].end(), 0)};
        m_F_n[point].push_back(F_k_n);

        m_r_n[point].push_back(0);
        for (int i {0}; i != n + 1; i++) {
            m_r_n[point][i] = static_cast<double>(m_f_n[point][i]) / F_k_n;
        }

        m_w_n[point].push_back(static_cast<double>(m_f_i[point]) / m_steps);
    }

    for (auto point: m_old_only_points) {
        m_f_n[point].push_back(0);
        m_F_n[point].push_back(m_F_n[point].back());
        m_r_n[point].push_back(0);
        m_w_n[point].push_back(0);
    }

    // Setup vectors for new points
    for (auto point: m_new_points) {
        vector<double> p_k_n(n);
        m_p_n[point] = p_k_n;
    }
    for (auto point: m_S_n) {
        m_p_n[point].push_back(m_p_i[point]);
    }
}

void MezeiUSGCMCSimulation::update_bias(int n) {
    update_grids(n);
    estimate_normalizations(n);
    for (auto point: m_S_n) {

        // No T to be consistent with biases here
        double old_bias {m_E_w[point]};
        double new_bias {std::log(m_lP_n[point])};
        double updated_bias {new_bias};
        double D_bias = new_bias - old_bias;

        // Limit how quckly bias can change
        if (abs(D_bias) > m_max_D_bias) {
            if (D_bias > 0) {
                updated_bias = old_bias + m_max_D_bias;
            }
            else {
                updated_bias = old_bias - m_max_D_bias;
            }
        }

        m_E_w[point] = updated_bias;
    }
    m_grid_bias->replace_biases(m_E_w);
}

void MezeiUSGCMCSimulation::output_summary(int n) {
    *m_us_stream << "Iteration: " << n << "\n";
    *m_us_stream << "Normalization constants:\n";
    for (auto N: m_N) {
        *m_us_stream << N << " ";
    }
    *m_us_stream << "\n";
    *m_us_stream << "Gridpoint w, r, p, P, E:\n";
    for (auto point: m_S_n) {
        for (auto coor: point) {
            *m_us_stream << coor << " ";
        }
        *m_us_stream << std::setprecision(3);
        *m_us_stream << ": ";
        *m_us_stream << std::setw(10) << m_w_n[point][n];
        *m_us_stream << std::setw(10) << m_r_n[point][n];
        *m_us_stream << std::setw(10) << m_p_n[point][n];
        *m_us_stream << std::setw(10) << m_lP_n[point];
        *m_us_stream << std::setw(10) << m_E_w[point];
        *m_us_stream << "\n";
    }
    *m_us_stream << "\n";
}

void MezeiUSGCMCSimulation::estimate_normalizations(int n) {
    // Ugly special case
    if (n == 0) {
        m_N.push_back(1);
    }
    else {
        m_N.push_back(estimate_initial_normalization(n));
        estimate_final_normalizations();
    }

    // Update probabilities
    m_old_lP_n = m_lP_n;
    double P_n_sum {0};
    for (auto point: m_S_n) {
        for (int i {0}; i != n + 1; i++) {
            P_n_sum += m_r_n[point][i] * m_N[i] * m_p_n[point][i];
        }
    }
    for (auto point: m_S_n) {
        double P_k_n {0};
        for (int i {0}; i != n + 1; i++) {
            P_k_n += m_r_n[point][i] * m_N[i] * m_p_n[point][i];
        }
        m_P_n[point] = P_k_n;
        m_lP_n[point] = P_k_n / P_n_sum;
    }
}

double MezeiUSGCMCSimulation::estimate_initial_normalization(int n) {
    // Use a linear one-step optimization
    // Uses the previous Nis and estimates the current iteration's Ni
    
    // N_n = a(c)/b(c)

    // Calculate grid of c
    GridFloats c_grid {};
    long int F_prev_sum {n * m_steps};

    for (auto point: m_S_n) {
        double w_n_k {m_w_n[point][n]};
        double r_n_k {m_r_n[point][n]};
        double c {pow(r_n_k, 2) * m_F_n[point][n - 1] / F_prev_sum +
                w_n_k*pow((1 - r_n_k), 2)};
        c_grid[point] = c;
    }

    // Calculate a(c)
    double N_n_num {0};
    for (auto point: m_S_n) {
        N_n_num += m_P_n[point] * m_p_n[point][n] * c_grid[point];
    }

    // Calculate b(c)
    double N_n_denom {0};
    for (auto point: m_S_n) {
        N_n_denom += pow(m_p_n[point][n], 2) * c_grid[point];
    }

    return N_n_num / N_n_denom;
}

void MezeiUSGCMCSimulation::estimate_final_normalizations() {
    // With relative deviation square sum
    int num_iters {static_cast<int>(m_N.size())};
    double num_iters_as_double {static_cast<double>(m_N.size())};
    Problem problem;
    vector<int> parameter_block_sizes {num_iters - 1, num_iters, num_iters, 1, 1};
    for (int n {0}; n != num_iters; n++) {
        double n_as_double {static_cast<double>(n)};
        for (auto point: m_s_n[n]) {
            CostFunction* cost_function =
                    new RDevSquareSumCostFunction {1, parameter_block_sizes};

            /* Autodiff
            //DynamicAutoDiffCostFunction<RDevSquareSum, 4>* cost_function =
            //    new DynamicAutoDiffCostFunction<RDevSquareSum, 4>(
            //            new RDevSquareSum {});

            //cost_function->AddParameterBlock(num_iters);
            //cost_function->AddParameterBlock(num_iters);
            //cost_function->AddParameterBlock(num_iters);
            //cost_function->AddParameterBlock(1);
            //cost_function->AddParameterBlock(1);
            //cost_function->SetNumResiduals(1);
            */

            LossFunction* loss_function = new ScaledLoss {NULL, m_w_n[point][n],
                    ceres::TAKE_OWNERSHIP};
            // I am passing the core arrays that are held in vectors (&nvariable[0])
            problem.AddResidualBlock(cost_function, loss_function, &m_N[1],
                    &m_r_n[point][0], &m_p_n[point][0], &n_as_double, &num_iters_as_double);
            problem.SetParameterBlockConstant(&m_r_n[point][0]);
            problem.SetParameterBlockConstant(&m_p_n[point][0]);
        }
        problem.SetParameterBlockConstant(&n_as_double);
        if (n != 0) {
            problem.SetParameterLowerBound(&m_N[1], n - 1, 0);
        }
    }
    problem.SetParameterBlockConstant(&num_iters_as_double);

    Solver::Options options;
    options.max_num_iterations = 100;
    options.minimizer_type = ceres::TRUST_REGION;
    options.minimizer_progress_to_stdout = false;

    Solver::Summary summary;
    Solve(options, &problem, &summary);
    m_solver_file << "Iteration " << (num_iters - 1) << "\n";
    m_solver_file << summary.FullReport() << "\n\n";
}

// For manual diff
RDevSquareSumCostFunction::RDevSquareSumCostFunction(int num_residuals,
        vector<int>& parameter_block_sizes) {
    set_num_residuals(num_residuals);
    vector<int>* parameter_block_sizes_ {mutable_parameter_block_sizes()};
    *parameter_block_sizes_ = parameter_block_sizes;
}

// For manual diff
bool RDevSquareSumCostFunction::Evaluate(
        double const* const* params, // 0: N_n[1:]; 1: r_k_n; 2: p_k_n; 3: iter; 4 num_iters
        double* residuals,
        double** jacobians) const {

    // Calculate cost function
    double P_k_n {params[1][0] * params[2][0]};
    int num_iters {static_cast<int>(params[4][0])};
    for (int j {0}; j != num_iters - 1; j ++) {
        P_k_n += params[1][j + 1] * params[0][j] * params[2][j + 1];
    }
    
    int iter {static_cast<int>(params[3][0])};
    if (iter == 0) {
        residuals[0] = (params[2][0] - P_k_n) / P_k_n;
    }
    else {
        residuals[0] = (params[0][iter - 1]*params[2][iter] - P_k_n) / P_k_n;
    }

    // Calculate Jacobian of cost function
    if (jacobians != NULL && jacobians[0] != NULL) {
        double P_k_n_sqr {pow(P_k_n, 2)};
        if (iter == 0) {
            for (int j {0}; j != num_iters - 1; j++) {
                double term_two = -params[2][iter]*params[1][j + 1]*
                        params[2][j + 1] / P_k_n_sqr;
                jacobians[0][j] = term_two;
            }
        }
        else {
            for (int j {0}; j != num_iters - 1; j++) {
                double term_two = -params[0][iter - 1]*params[2][iter]*
                        params[1][j + 1]*params[2][j + 1] / P_k_n_sqr;
                if (j + 1 == iter) {
                    jacobians[0][j] = params[2][iter] / P_k_n + term_two;
                }
                else {
                    jacobians[0][j] = term_two;
                }
            }
        }
    }
    return true;
}

// For autodiff
template <typename T>
bool RDevSquareSum::operator() (
        T const* const* params, // 0, N_n; 1, r_n_k; 2 p_n_k; 3, i; 4 num_iters
        T* residual) const {
    residual[0] = params[1][0] * params[0][0] * params[2][0];
    for (double j {1}; j != params[4][0]; j += 1) {
        int j_int {static_cast<int>(j)};
        residual[0] += params[1][j_int] * params[0][j_int] * params[2][j_int];
    }
    
    // Hack way of getting an index from the passed parameter
    double i;
    for (i = 0; i != params[3][0]; i++) {}
    int i_int = static_cast<int>(i);
    residual[0] = (params[0][i_int] * params[2][i_int] - residual[0]) / residual[0];

    return true;
}
