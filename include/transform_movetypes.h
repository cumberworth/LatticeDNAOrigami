// transform_movetypes.h

#ifndef TRANSFORM_MOVETYPES_H
#define TRANSFORM_MOVETYPES_H

#include "cb_movetypes.h"
#include "rg_movetypes.h"

namespace movetypes {

    /**
      * Base class for transformation and linker regrowth movetypes
      *
      * Transformation of a scaffold segment and any bound staples,
      * followed by regrowth of linker regions between it and the rest of the
      * system.
      */
    class LinkerRegrowthMCMovetype:
        virtual public CTRegrowthMCMovetype {
        public:
            LinkerRegrowthMCMovetype(
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
                    unsigned int max_regrowth,
                    unsigned int max_linker_length);
            LinkerRegrowthMCMovetype(const
                    LinkerRegrowthMCMovetype&) = delete;
            LinkerRegrowthMCMovetype& operator=(const
                    LinkerRegrowthMCMovetype&) = delete;

            void write_log_summary(ostream* log_entry) override;

        protected:
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

            virtual void revert_transformation(
                    vector<Domain*> central_domains) = 0;

            int m_max_disp;
            int m_max_turns;
            unsigned int m_max_linker_length;

            CTCBLinkerRegrowthTracking m_tracker {};
            unordered_map<CTCBLinkerRegrowthTracking, MovetypeTracking> m_tracking {};
            vector<Domain*> m_linker_endpoints;
    };

    /**
      * Base class for clustered transformation and linker regrowth movetypes
      *
      * Instead of randomly selecting a segment of the scaffold for
      * transformation, a segment is selected in which all domains are in
      * bound states, and the linker regions are in unbound (or misbound?)
      * states.
      */
    class ClusteredLinkerRegrowthMCMovetype:
        virtual public LinkerRegrowthMCMovetype {

        public:
            ClusteredLinkerRegrowthMCMovetype(
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
                    unsigned int max_linker_length);
            ClusteredLinkerRegrowthMCMovetype(const
                    ClusteredLinkerRegrowthMCMovetype&) = delete;
            ClusteredLinkerRegrowthMCMovetype& operator=(const
                    ClusteredLinkerRegrowthMCMovetype&) = delete;

        private:
            set<int> select_and_setup_segments(
                    vector<Domain*>& linker1,
                    vector<Domain*>& linker2,
                    vector<Domain*>& central_segment) override;
            vector<Domain*> cyclic_select_central_segment();
            vector<Domain*> linear_select_central_segment();
    };

    /**
      * CTCB linker regrowth movetype
      */
    class CTCBLinkerRegrowthMCMovetype:
        virtual public LinkerRegrowthMCMovetype,
        virtual public CTCBRegrowthMCMovetype {

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
                    int max_turns,
                    unsigned int max_regrowth,
                    unsigned int max_linker_length);
            CTCBLinkerRegrowthMCMovetype(const
                    CTCBLinkerRegrowthMCMovetype&) = delete;
            CTCBLinkerRegrowthMCMovetype& operator=(const
                    CTCBLinkerRegrowthMCMovetype&) = delete;

            bool internal_attempt_move() override;

        protected:
            void revert_transformation(
                    vector<Domain*> central_domains) override;
            void reset_internal() override;
    };

    /**
      * CTCB clustered linker regrowth movetype
      */
    class CTCBClusteredLinkerRegrowthMCMovetype:
        public CTCBLinkerRegrowthMCMovetype,
        public ClusteredLinkerRegrowthMCMovetype {

        public:
            CTCBClusteredLinkerRegrowthMCMovetype(
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
                    unsigned int max_linker_length);
            CTCBClusteredLinkerRegrowthMCMovetype(const
                    CTCBClusteredLinkerRegrowthMCMovetype&) = delete;
            CTCBClusteredLinkerRegrowthMCMovetype& operator=(const
                    CTCBClusteredLinkerRegrowthMCMovetype&) = delete;
            void reset_internal() override;
    };

    /**
      * CTRG linker regrowth movetype
      */
    class CTRGLinkerRegrowthMCMovetype:
        virtual public LinkerRegrowthMCMovetype,
        virtual public CTRGRegrowthMCMovetype {

        public:
            CTRGLinkerRegrowthMCMovetype(
                    OrigamiSystem& origami_system,
                    RandomGens& random_gens,
                    IdealRandomWalks& ideal_random_walks,
                    vector<OrigamiOutputFile*> config_files,
                    string label,
                    SystemOrderParams& ops,
                    SystemBiases& biases,
                    InputParameters& params,
                    int num_excluded_staples,
                    int max_num_recoils,
                    int max_c_attempts,
                    int max_disp,
                    int max_turns,
                    unsigned int max_regrowth,
                    unsigned int max_linker_length);
            CTRGLinkerRegrowthMCMovetype(const
                    CTRGLinkerRegrowthMCMovetype&) = delete;
            CTRGLinkerRegrowthMCMovetype& operator=(const
                    CTRGLinkerRegrowthMCMovetype&) = delete;

            bool internal_attempt_move() override;

        protected:
            void revert_transformation(
                    vector<Domain*> central_domains) override;
            void reset_internal() override;
    };
}

#endif // TRANSFORM_MOVETYPES_H
