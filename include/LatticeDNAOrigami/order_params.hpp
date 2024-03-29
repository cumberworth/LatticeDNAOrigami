// order_params.h

#ifndef ORDER_PARAMS_H
#define ORDER_PARAMS_H

#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "LatticeDNAOrigami/domain.hpp"
#include "LatticeDNAOrigami/hash.hpp"
#include "LatticeDNAOrigami/origami_system.hpp"
#include "LatticeDNAOrigami/parser.hpp"
#include "LatticeDNAOrigami/utility.hpp"

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
using utility::Occupancy;
using utility::VectorThree;

// System parameters calculated from current config
class OrderParam {
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

// Distance between two given domains
class DistOrderParam: public OrderParam {
  public:
    DistOrderParam(Domain& domain_1, Domain& domain_2, string label);

    int calc_param() override final;
    int check_param(Domain& domain, VectorThree new_pos, VectorThree, Occupancy)
            override final;

  private:
    Domain& m_domain_1;
    Domain& m_domain_2;
};

// If given domain pair occupies adjacent lattice sites
class AdjacentSiteOrderParam: public OrderParam {
  public:
    AdjacentSiteOrderParam(Domain& domain_1, Domain& domain_2, string label);

    int calc_param() override final;
    int check_param(Domain& domain, VectorThree new_pos, VectorThree, Occupancy)
            override final;

  private:
    Domain& m_domain_1;
    Domain& m_domain_2;
};

class SumOrderParam: public OrderParam {
  public:
    SumOrderParam(vector<reference_wrapper<OrderParam>> ops, string label);

    int calc_param() override final;
    int check_param(Domain&, VectorThree, VectorThree, Occupancy)
            override final;

  private:
    vector<reference_wrapper<OrderParam>> m_ops {};
};

class NumStaplesOrderParam: public OrderParam {
  public:
    NumStaplesOrderParam(OrigamiSystem& origami, string label);
    int calc_param() override final;
    int check_param(Domain& domain, VectorThree new_pos, VectorThree, Occupancy)
            override final;

  private:
    OrigamiSystem& m_origami;
};

class NumStaplesTypeOrderParam: public OrderParam {
  public:
    NumStaplesTypeOrderParam(OrigamiSystem& origami, int stype, string label);
    int calc_param() override final;
    int check_param(Domain& domain, VectorThree new_pos, VectorThree, Occupancy)
            override final;

  private:
    OrigamiSystem& m_origami;
    int m_c_ident;
};

class StapleTypeFullyBoundOrderParam: public OrderParam {
  public:
    StapleTypeFullyBoundOrderParam(
            OrigamiSystem& origami,
            int stype,
            string label);
    int calc_param() override final;
    int check_param(Domain& domain, VectorThree new_pos, VectorThree, Occupancy)
            override final;

  private:
    OrigamiSystem& m_origami;
    int m_c_ident;
};

class NumBoundDomainPairsOrderParam: public OrderParam {
  public:
    NumBoundDomainPairsOrderParam(OrigamiSystem& origami, string label);
    int calc_param() override final;
    int check_param(Domain& domain, VectorThree new_pos, VectorThree, Occupancy)
            override final;

  private:
    OrigamiSystem& m_origami;
};

class NumMisboundDomainPairsOrderParam: public OrderParam {
  public:
    NumMisboundDomainPairsOrderParam(OrigamiSystem& origami, string label);
    int calc_param() override final;
    int check_param(Domain& domain, VectorThree new_pos, VectorThree, Occupancy)
            override final;

  private:
    OrigamiSystem& m_origami;
};

class NumStackedPairsOrderParam: public OrderParam {
  public:
    NumStackedPairsOrderParam(OrigamiSystem& origami, string label);
    int calc_param() override final;
    int check_param(Domain& domain, VectorThree, VectorThree, Occupancy)
            override final;

  private:
    OrigamiSystem& m_origami;
};

class NumLinearHelicesOrderParam: public OrderParam {
  public:
    NumLinearHelicesOrderParam(OrigamiSystem& origami, string label);
    int calc_param() override final;
    int check_param(Domain& domain, VectorThree, VectorThree, Occupancy)
            override final;

  private:
    OrigamiSystem& m_origami;
};

class NumStackedJunctsOrderParam: public OrderParam {
  public:
    NumStackedJunctsOrderParam(OrigamiSystem& origami, string label);
    int calc_param() override final;
    int check_param(Domain& domain, VectorThree, VectorThree, Occupancy)
            override final;

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
    void setup_ops(string ops_filename, vector<pair<int, int>> keys);

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
} // namespace orderParams

#endif // ORDER_PARAMS_H
