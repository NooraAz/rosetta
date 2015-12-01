// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/farna/RNA_FragmentMonteCarlo.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_farna_RNA_FragmentMonteCarlo_HH
#define INCLUDED_protocols_farna_RNA_FragmentMonteCarlo_HH

#include <protocols/moves/Mover.hh>
#include <protocols/farna/RNA_FragmentMonteCarlo.fwd.hh>
#include <protocols/farna/RNA_FragmentMonteCarloOptions.fwd.hh>
#include <protocols/farna/RNA_FragmentMover.fwd.hh>
#include <protocols/farna/RNA_Fragments.fwd.hh>
#include <protocols/farna/RNA_StructureParameters.fwd.hh>
#include <protocols/farna/RNA_ChunkLibrary.fwd.hh>
#include <protocols/farna/RNA_LoopCloser.fwd.hh>
#include <protocols/farna/RNA_Minimizer.fwd.hh>
#include <protocols/farna/RNA_Relaxer.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/rigid/RigidBodyMover.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/types.hh>
#include <ObjexxFCL/format.hh>

namespace protocols {
namespace farna {

	class RNA_FragmentMonteCarlo: public protocols::moves::Mover {

	public:

		//constructor
		RNA_FragmentMonteCarlo( RNA_FragmentMonteCarloOptionsCOP = 0 );

		//destructor
		~RNA_FragmentMonteCarlo();

		/// @brief Apply the loop-rebuild protocol to the input pose
		void apply( core::pose::Pose & pose );

		std::string get_name() const { return "RNA_FragmentMonteCarlo"; }

		void
		set_refine_pose( bool const setting ){ refine_pose_ = setting; }

		void
		set_chunk_coverage( core::Real const & setting ){ chunk_coverage_ = setting; }

		void
		set_all_lores_score_final( std::list< core::Real > const & setting ){ all_lores_score_final_ = setting; }

		void
		set_denovo_scorefxn( core::scoring::ScoreFunctionCOP setting ){ denovo_scorefxn_ = setting; }

		void
		set_hires_scorefxn( core::scoring::ScoreFunctionCOP setting ){ hires_scorefxn_ = setting; }

		void
		set_out_file_tag( std::string const & setting ){ out_file_tag_ = setting; }

		std::list< core::Real > const & all_lores_score_final() const { return all_lores_score_final_; }

		core::pose::PoseOP lores_pose() { return lores_pose_; }

		core::Real lores_score_early() const { return lores_score_early_; }

		core::Real lores_score_final() const { return lores_score_final_; }

		RNA_StructureParametersCOP rna_structure_parameters() const { return rna_structure_parameters_; }
		RNA_ChunkLibraryCOP rna_chunk_library() const { return rna_chunk_library_; }

		void show(std::ostream & output) const;

	private:

		void
		initialize( core::pose::Pose & pose );

		void
		initialize_libraries( core::pose::Pose & pose );

		void
		initialize_movers();

		void
		initialize_score_functions();

		void
		initialize_parameters();

		void
		do_random_moves( core::pose::Pose & pose );

		void
		randomize_rigid_body_orientations( core::pose::Pose & pose );

		void
		update_denovo_scorefxn_weights( Size const r );

		Size
		figure_out_constraint_separation_cutoff( Size const r, Size const  max_dist );

		void
		update_pose_constraints( Size const r, core::pose::Pose & pose );

		void
		update_frag_size( Size const r );

		void
		random_fragment_trial( core::pose::Pose & pose );

		bool
		random_chunk_trial( core::pose::Pose & pose );

		void
		random_jump_trial( core::pose::Pose & pose );

		void
		RNA_move_trial( core::pose::Pose & pose );

		void
		setup_rigid_body_mover( core::pose::Pose const & pose, core::Size const r );

		bool
		check_score_filter( core::Real const lores_score_, std::list< core::Real > & all_lores_score_ );

		void
		apply_chem_shift_data( core::pose::Pose & pose );

		void
		setup_monte_carlo_cycles( core::pose::Pose const & pose );

		void
		final_score( core::pose::Pose & pose );

	private:

		// The parameters in this OptionsCOP should not change:
		RNA_FragmentMonteCarloOptionsCOP options_;
		std::string out_file_tag_;

		// Movers (currently must be set up outside, but should write auto-setup code)
		RNA_StructureParametersOP rna_structure_parameters_;
		RNA_ChunkLibraryOP rna_chunk_library_;
		RNA_FragmentMoverOP rna_fragment_mover_;
		RNA_LoopCloserOP rna_loop_closer_;
		RNA_MinimizerOP rna_minimizer_;
		RNA_RelaxerOP rna_relaxer_;
		protocols::rigid::RigidBodyPerturbMoverOP rigid_body_mover_;

		core::scoring::ScoreFunctionCOP denovo_scorefxn_;
		core::scoring::ScoreFunctionCOP hires_scorefxn_;
		core::scoring::ScoreFunctionCOP chem_shift_scorefxn_;
		core::scoring::ScoreFunctionCOP final_scorefxn_;
		core::scoring::ScoreFunctionOP working_denovo_scorefxn_;

		// Parameters that change during run:
		protocols::moves::MonteCarloOP monte_carlo_;
		Size const monte_carlo_cycles_max_default_;
		Size monte_carlo_cycles_;
		Size rounds_;
		Size frag_size_;
		bool do_close_loops_;
		bool refine_pose_;
		core::Real jump_change_frequency_;
		core::Real lores_score_early_;
		core::Real lores_score_final_;
		core::Real chunk_coverage_;

		std::list< core::Real > all_lores_score_final_;
		core::scoring::constraints::ConstraintSetOP constraint_set_;
		core::pose::PoseOP lores_pose_;

	};

} //farna
} //protocols

#endif