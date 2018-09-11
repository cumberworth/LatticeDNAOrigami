// cb_movetypes.h

#ifndef CB_MOVETYPES_H
#define CB_MOVETYPES_H

#include <iostream>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "domain.h"
#include "files.h"
#include "hash.h"
#include "ideal_random_walk.h"
#include "movetypes.h"
#include "origami_system.h"
#include "parser.h"
#include "random_gens.h"
#include "top_constraint_points.h"
#include "utility.h"

/** Movetypes involving configuration bias */
namespace movetypes {

    using std::ostream;
    using std::pair;
    using std::set;
    using std::string;
    using std::unordered_map;
    using std::vector;

    using domainContainer::Domain;
    using files::OrigamiInputFile;
    using files::OrigamiOutputFile;
    using idealRandomWalk::IdealRandomWalks;
    using movetypes::MovetypeTracking;
    using movetypes::RegrowthMCMovetype;
    using movetypes::add_tracker;
    using origami::OrigamiSystem;
    using parser::InputParameters;
    using randomGen::RandomGens;
    using topConstraintPoints::Constraintpoints;
    using utility::VectorThree;
    using utility::StapleRegrowthTracking;
    using utility::CTCBScaffoldRegrowthTracking;
    using utility::CTCBLinkerRegrowthTracking;
    using utility::VectorThree;

    typedef vector<pair<VectorThree, VectorThree>> configsT;

    /** Configurational bias base class */
    class CBMCMovetype:
        virtual public RegrowthMCMovetype {

        public:
            CBMCMovetype(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks,
                    vector<OrigamiOutputFile*> config_files,
                    string label,
                    SystemOrderParams& ops,
                    SystemBiases& biases,
                    InputParameters& params);
            CBMCMovetype(const CBMCMovetype&) = delete;
            CBMCMovetype& operator=(const CBMCMovetype&) = delete;

            virtual void reset_internal() override;

        protected:

            /** Include external bias on whole configuration in Rosenbluth */
            void add_external_bias() override;

            /** Calculate the CB trial weights and rosenbluth-type weight */
            virtual vector<double> calc_bias(
                    const vector<double> bfactors, // Boltzman weights
                    const configsT& configs, // Configs
                    Domain* domain) = 0; // Domain being regrown

            /** Calculate the Boltzmann weights for all configurations
              *
              * Fills the configuration and weights arrays with the Boltzmann
              * weights for all configurations available to the given domain.
              * If all orientations are permissible for a given position, it
              * will output zero-orientation and multiply the weight by the
              * number of possible orientations.
              */
            void calc_biases(
                    const VectorThree p_prev, // Position of growthpoint domain
                    Domain& domain, // Domain to calculate weights for
                    configsT& configs, // Configs considered
                    vector<double>& bfactors); // Weights for configs

            /** Select and set a configuration for a domain using CB */
            void select_and_set_config(
                    const int i, // Index in segment array of domain to regrow
                    vector<Domain*> domains); // Segment being regrown

            /** Select and set configuration from given configs and weights */
            void select_and_set_new_config(
                    const vector<double> weights, // Weights
                    const configsT configs, // Config
                    Domain& domain); // Domain being regrown

            /** Set configuration to old */
            void select_and_set_old_config(Domain& domain);

            /** Set growthpoint without checking constraints */
            double set_old_growth_point(
                    Domain& growth_domain_new,
                    Domain& growth_domain_old);

            /** Test move acceptance and take appropriate action
              *
              * If accepted it will revert to the new configuration (as the
              * system will currently be in the old configuration), otherwise
              * it will clear the reset arrays so that reset origami has no
              * effect.
              */
            bool test_cb_acceptance();

            /** Unassign all given domains and store old configurations */
            double unassign_domains(vector<Domain*>);

            /** Prepare interals for regrowing old configuration */
            void setup_for_regrow_old();

            /** Update the external bias without adding */
            void update_external_bias();

            long double m_bias {1}; // Rosenbluth-type weight
            long double m_new_bias {1}; // Storage for new config's Rosenbluth
            double m_new_modifier {1}; // Storage for new config's modifier
            bool m_regrow_old {false}; // Regrowing old configuration
            unordered_map<pair<int, int>, VectorThree> m_old_pos {};
            unordered_map<pair<int, int>, VectorThree> m_old_ore {};
    };

    /** CB regrowth of staples */
    class CBStapleRegrowthMCMovetype:
        public CBMCMovetype {

        public:
            CBStapleRegrowthMCMovetype(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks,
                    vector<OrigamiOutputFile*> config_files,
                    string label,
                    SystemOrderParams& ops,
                    SystemBiases& biases,
                    InputParameters& params);
            CBStapleRegrowthMCMovetype(const
                    CBStapleRegrowthMCMovetype&) = delete;
            CBStapleRegrowthMCMovetype& operator=(const
                    CBStapleRegrowthMCMovetype&) = delete;

            void write_log_summary(ostream* log_entry) override;

        private:
            bool internal_attempt_move() override;
            virtual void add_tracker(bool accepted) override;

            /** Grow out domains from first in given array */
            void grow_chain(vector<Domain*> domains) override;

            /** Calculate the CB trial weights and rosenbluth-type weight */
            vector<double> calc_bias(
                    const vector<double> bfactors, // Boltzman weights
                    const configsT& configs, // Configs
                    Domain* domain) override; // Domain being regrown

            /** Set given growthpoint and grow staple */
            void set_growthpoint_and_grow_staple(
                    domainPairT growthpoint,
                    vector<Domain*> selected_chain);

            StapleRegrowthTracking m_tracker {};
            unordered_map<StapleRegrowthTracking, MovetypeTracking>
                    m_tracking {};
    };

    /** Base for conserved topology CB regrowth moves */
    class CTCBRegrowthMCMovetype:
        virtual public CBMCMovetype,
        virtual public CTRegrowthMCMovetype {

        public:
            CTCBRegrowthMCMovetype(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks,
                    vector<OrigamiOutputFile*> config_files,
                    string label,
                    SystemOrderParams& ops,
                    SystemBiases& biases,
                    InputParameters& params,
                    int num_excluded_staples,
                    int max_regrowth);
            CTCBRegrowthMCMovetype(const
                    CTCBRegrowthMCMovetype&) = delete;
            CTCBRegrowthMCMovetype& operator=(const
                    CTCBRegrowthMCMovetype&) = delete;

            virtual void reset_internal() override;

        protected:

            /** Grow out domains from first in given array */
            void grow_chain(vector<Domain*> domains) final override;

            /** Calculate the CB trial weights and rosenbluth-type weight */
            vector<double> calc_bias(
                    const vector<double> bfactors, // Boltzman weights
                    const configsT& configs, // Configs
                    Domain* domain) override; // Domain being regrown

            /** Grow given staple and update fixed-end biases */
            void grow_staple_and_update_endpoints(
                    Domain* growth_domain_old);
    };

    /** CTCB regrowth of a contiguous scaffold segment and bound staples */
    class CTCBScaffoldRegrowthMCMovetype:
        public CTCBRegrowthMCMovetype {

        public:
            CTCBScaffoldRegrowthMCMovetype(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks,
                    vector<OrigamiOutputFile*> config_files,
                    string label,
                    SystemOrderParams& ops,
                    SystemBiases& biases,
                    InputParameters& params,
                    int num_excluded_staples,
                    int max_regrowth);
            CTCBScaffoldRegrowthMCMovetype(const
                    CTCBScaffoldRegrowthMCMovetype&) = delete;
            CTCBScaffoldRegrowthMCMovetype& operator=(const
                    CTCBScaffoldRegrowthMCMovetype&) = delete;

            void write_log_summary(ostream* log_entry) override;

        private:
            bool internal_attempt_move() override;
            virtual void add_tracker(bool accepted) override;

            CTCBScaffoldRegrowthTracking m_tracker {};
            unordered_map<CTCBScaffoldRegrowthTracking, MovetypeTracking>
                    m_tracking {};
    };

    /** CTCB regrowth of noncontiguous scaffold segments and bound staples */
    class CTCBJumpScaffoldRegrowthMCMovetype:
        public CTCBRegrowthMCMovetype {

        public:
            CTCBJumpScaffoldRegrowthMCMovetype(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks,
                    vector<OrigamiOutputFile*> config_files,
                    string label,
                    SystemOrderParams& ops,
                    SystemBiases& biases,
                    InputParameters& params,
                    int num_excluded_staples,
                    int max_regrowth,
                    int max_seg_regrowth);
            CTCBJumpScaffoldRegrowthMCMovetype(const
                    CTCBScaffoldRegrowthMCMovetype&) = delete;
            CTCBJumpScaffoldRegrowthMCMovetype& operator=(const
                    CTCBScaffoldRegrowthMCMovetype&) = delete;

            void write_log_summary(ostream* log_entry) override;

        private:
            bool internal_attempt_move() override;
            virtual void add_tracker(bool accepted) override;

            void set_first_seg_domain(Domain* seg);

            CTCBScaffoldRegrowthTracking m_tracker {};
            unordered_map<CTCBScaffoldRegrowthTracking, MovetypeTracking>
                    m_tracking {};
    };
}

#endif // CB_MOVETYPES_H
