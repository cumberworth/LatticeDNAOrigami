// transform_movetypes.h

#ifndef TRANSFORM_MOVETYPES_H
#define TRANSFORM_MOVETYPES_H

#include "cb_movetypes.h"

namespace movetypes {

    /**
      * Transformation of a segment and regrowth of its linkers
      */
    class CTCBLinkerRegrowthMCMovetype:
        public CTCBRegrowthMCMovetype {

        public:
            CTCBLinkerRegrowthMCMovetype(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks,
                    vector<OrigamiOutputFile*> config_files,
                    string label,
                    SystemOrderParams& ops,
                    SystemBiases& biases,
                    InputParameters& params,
                    int num_excluded_staples,
                    int max_disp,
                    int max_turns);
            CTCBLinkerRegrowthMCMovetype(const
                    CTCBLinkerRegrowthMCMovetype&) = delete;
            CTCBLinkerRegrowthMCMovetype& operator=(const
                    CTCBLinkerRegrowthMCMovetype&) = delete;

            void write_log_summary(ostream* log_entry) override;
            void reset_internal() override;

        protected:
            bool internal_attempt_move() override;
            virtual void add_tracker(bool accepted) override;

            /** Unasssign domains reset arrays for resetting */
            void reset_segment(vector<Domain*> segment, size_t last_di);

            /**
              * Select segment to be transformed and both linkers
              *
              * Returns staples to be regrown. Also sets up the constraints for
              * the fixed end biases.
              */
            virtual set<int> select_and_setup_segments(
                    vector<Domain*>& linker1,
                    vector<Domain*>& linker2,
                    vector<Domain*>& central_segment);

            virtual pair<int, int> select_internal_endpoints(
                    vector<Domain*> domains);

            /**
              * Setup constraints for fixed end biases
              *
              * Returns staples to be regrown.
              */
            set<int> setup_fixed_end_biases(
                    vector<Domain*>& linker1,
                    vector<Domain*>& linker2);

            bool domains_bound_externally(vector<Domain*> domains);
            bool scan_for_external_scaffold_domain(
                    Domain* domain,
                    vector<Domain*> domains,
                    set<int>& participating_chains);

            /** Transform the central segment */
            void transform_segment(
                    vector<Domain*> linker1,
                    vector<Domain*> linker2,
                    vector<Domain*> central_segment,
                    vector<Domain*> central_domains);

            /** Translate and rotate segment */
            bool apply_transformation(
                    vector<Domain*> central_domains,
                    VectorThree disp,
                    VectorThree center,
                    VectorThree axis,
                    int turns);

            /** Check if enough steps to regrow linkers */
            bool steps_less_than_distance(
                    vector<Domain*> linker1,
                    vector<Domain*> linker2);

            void revert_transformation(vector<Domain*> central_domains);

            int m_max_disp;
            int m_max_turns;

            CTCBLinkerRegrowthTracking m_tracker {};
            unordered_map<CTCBLinkerRegrowthTracking, MovetypeTracking> m_tracking {};
            vector<Domain*> m_linker_endpoints;
    };

    /**
      * Transformation of a segment and regrowth of its linkers
      *
      * Constrasts with parent in selection of segment to be transformed. Will
      * select domains that are contiguously bound out from the seed domain.
      */
    class ClusteredCTCBLinkerRegrowth:
        public CTCBLinkerRegrowthMCMovetype {

        public:
            ClusteredCTCBLinkerRegrowth(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks,
                    vector<OrigamiOutputFile*> config_files,
                    string label,
                    SystemOrderParams& ops,
                    SystemBiases& biases,
                    InputParameters& params,
                    int num_excluded_staples,
                    int max_disp,
                    int max_turns);
            ClusteredCTCBLinkerRegrowth(const
                    ClusteredCTCBLinkerRegrowth&) = delete;
            ClusteredCTCBLinkerRegrowth& operator=(const
                    ClusteredCTCBLinkerRegrowth&) = delete;

        protected:
            pair<int, int> select_internal_endpoints(
                    vector<Domain*> domains) override;
    };

    /**
      * Transformation of a segment and regrowth of its linkers
      *
      * Constrasts with parent in selection of segment to be transformed. Will
      * select domains that are contiguously bound out from the seed domain.
      */
    class Clustered2CTCBLinkerRegrowth:
        public CTCBLinkerRegrowthMCMovetype {

        public:
            Clustered2CTCBLinkerRegrowth(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks,
                    vector<OrigamiOutputFile*> config_files,
                    string label,
                    SystemOrderParams& ops,
                    SystemBiases& biases,
                    InputParameters& params,
                    int num_excluded_staples,
                    int max_disp,
                    int max_turns,
                    int m_max_regrowth);
            Clustered2CTCBLinkerRegrowth(const
                    ClusteredCTCBLinkerRegrowth&) = delete;
            Clustered2CTCBLinkerRegrowth& operator=(const
                    ClusteredCTCBLinkerRegrowth&) = delete;

        private:
            set<int> select_and_setup_segments(
                    vector<Domain*>& linker1,
                    vector<Domain*>& linker2,
                    vector<Domain*>& central_segment) override;
            vector<Domain*> cyclic_select_internal_endpoints();
            vector<Domain*> linear_select_internal_endpoints();
    };
}

#endif // TRANSFORM_MOVETYPES_H
