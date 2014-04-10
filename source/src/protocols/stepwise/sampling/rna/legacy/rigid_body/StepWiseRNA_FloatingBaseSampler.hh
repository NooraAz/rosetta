// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/legacy/rigid_body/StepWiseRNA_FloatingBaseSampler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_FloatingBaseSampler_HH
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_FloatingBaseSampler_HH

#include <protocols/moves/MoverForPoseList.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_ModelerOptions.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Classes.hh>
#include <protocols/stepwise/sampling/rna/rigid_body/FloatingBaseClasses.hh>
#include <protocols/stepwise/sampling/rna/sugar/SugarModeling.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_PoseSelection.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_BaseCentroidChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_VDW_BinChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/AtrRepChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/ChainClosureChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/ChainClosableGeometryChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/o2prime/O2PrimePacker.fwd.hh>
#include <protocols/stepwise/sampling/rna/phosphate/MultiPhosphateSampler.fwd.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamer.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>


#ifdef WIN32
	#include <protocols/stepwise/sampling/rna/checker/ChainClosureChecker.hh>
	#include <protocols/stepwise/sampling/rna/checker/ChainClosableGeometryChecker.hh>
	#include <protocols/stepwise/sampling/rna/checker/AtrRepChecker.hh>
#endif

// using namespace core;
// using namespace protocols::stepwise::sampling::rna;
// using namespace protocols::stepwise::sampling::rna::rigid_body;
//
// Commented out because “using namespace X” in header files outside of class declaration is explicitly forbidden
// by our coding convention due to problems it create on modern compilers and because of the name clashing.
// For more information please see: https://wiki.rosettacommons.org/index.php/Coding_conventions#Using

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace legacy {
namespace rigid_body {

	class StepWiseRNA_FloatingBaseSampler: public protocols::moves::MoverForPoseList {

	public:

		//constructor
		StepWiseRNA_FloatingBaseSampler( StepWiseRNA_JobParametersCOP & job_parameters_ );

		//destructor
		~StepWiseRNA_FloatingBaseSampler();

		virtual void apply( core::pose::Pose & pose_to_visualize );

		virtual std::string get_name() const;

		using MoverForPoseList::apply;

	public:

		void set_silent_file( std::string const & setting ){ silent_file_ = setting; }

		void set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

		utility::vector1< core::pose::PoseOP > &		get_pose_list();

		void
		set_base_centroid_checker( checker::RNA_BaseCentroidCheckerOP & checker );

		void
		set_user_input_VDW_bin_checker( checker::RNA_VDW_BinCheckerOP const & user_input_VDW_bin_checker );

		void
		set_anchor_sugar_modeling( sugar::SugarModeling const & anchor_sugar_modeling );

		void
		set_options( StepWiseRNA_ModelerOptionsCOP options );

	private:

		void
		initialize_poses_and_stubs_and_checkers( core::pose::Pose & pose  );

		void
		reinstantiate_backbone_and_add_constraint_at_moving_res(core::pose::Pose & pose, core::Size const & five_prime_chain_break_res );

		void
		initialize_euler_angle_grid_parameters();

		void
		initialize_xyz_grid_parameters();

		void
		initialize_rigid_body_sampler(core::pose::Pose const & pose );

		bool
		break_early_for_integration_tests();

		core::Size
		screen_anchor_sugar_conformation(core::pose::Pose & pose, std::string & tag );

		core::kinematics::Stub
		get_reference_stub(core::pose::Pose const & pose ) const;

		void
		initialize_other_residues_base_list(core::pose::Pose const & pose );

		void
		update_base_bin_map(protocols::stepwise::sampling::rna::rigid_body::BaseBin const & base_bin );

		void
		update_base_bin_map( utility::vector1< core::Real > const & rigid_body_values );

		void
		output_count_data();

		void
		instantiate_moving_sugar_and_o2prime(core::pose::Pose & pose );

		void
		virtualize_moving_sugar_and_o2prime(core::pose::Pose & pose );

		bool
		check_moving_sugar(core::pose::Pose & pose );

	private:

		StepWiseRNA_JobParametersCOP job_parameters_; //need to use the full_to_sub map...should convert to const style.. Parin Feb 28, 2010
		core::Size const moving_res_; // Might not corresponds to user input.
		core::Size const moving_suite_; // dofs betweeen this value and value+1 actually move.
		bool const is_prepend_;
		bool const is_internal_; // no cutpoints before or after moving_res.
		core::Size const gap_size_; /* If this is zero or one, need to screen or closable chain break */
		core::Size const gap_size_to_anchor_;
		core::Size const five_prime_chain_break_res_;
		core::Size const chain_break_reference_res_;
		core::Size const reference_res_; //the last static_residues that this attach to the moving residues
		core::Size const floating_base_five_prime_chain_break_;
		core::Size const floating_base_three_prime_chain_break_;
		bool const is_dinucleotide_;
		bool const close_chain_to_distal_;
		bool const close_chain_to_anchor_;

		utility::vector1< core::pose::PoseOP > pose_list_;

		core::scoring::ScoreFunctionOP scorefxn_;

		protocols::rotamer_sampler::rigid_body::RigidBodyRotamerOP sampler_;

		StepWiseRNA_CountStruct count_data_;
		std::string silent_file_;
		bool native_rmsd_screen_;
		core::Size num_pose_kept_to_use_;

		StepWiseRNA_ModelerOptionsCOP options_;

		sugar::SugarModeling anchor_sugar_modeling_;

		checker::AtrRepCheckerOP atr_rep_checker_;
		checker::AtrRepCheckerOP atr_rep_checker_with_instantiated_sugar_;
		utility::vector1< checker::AtrRepCheckerOP > atr_rep_checkers_for_anchor_sugar_models_;
		checker::RNA_VDW_BinCheckerOP VDW_bin_checker_;
		checker::RNA_VDW_BinCheckerOP user_input_VDW_bin_checker_;
		checker::ChainClosableGeometryCheckerOP chain_closable_geometry_to_distal_checker_, chain_closable_geometry_to_anchor_checker_;
		checker::ChainClosureCheckerOP chain_break_to_distal_checker_, chain_break_to_anchor_checker_;
		checker::RNA_BaseCentroidCheckerOP base_centroid_checker_;

		core::pose::PoseOP screening_pose_, sugar_screening_pose_;
		utility::vector1 < core::conformation::ResidueOP > moving_rsd_at_origin_list_, screening_moving_rsd_at_origin_list_, sugar_screening_moving_rsd_at_origin_list_;
		StepWiseRNA_PoseSelectionOP pose_selection_;

		o2prime::O2PrimePackerOP o2prime_packer_;
		phosphate::MultiPhosphateSamplerOP phosphate_sampler_;

		core::kinematics::Stub reference_stub_;
		utility::vector1 < core::kinematics::Stub > other_residues_base_list_;
		std::map< protocols::stepwise::sampling::rna::rigid_body::BaseBin, int, protocols::stepwise::sampling::rna::rigid_body::compare_base_bin > base_bin_map_;

		core::Size max_ntries_; // for choose_random

		int euler_angle_bin_min_;
		int euler_angle_bin_max_;
		core::Real euler_angle_bin_size_;
		int euler_z_bin_min_;
		int euler_z_bin_max_;
		core::Real euler_z_bin_size_;
		int centroid_bin_min_;
		int centroid_bin_max_;
		core::Real centroid_bin_size_;
		core::Distance max_distance_;
		core::Real max_distance_squared_;

		bool try_sugar_instantiation_;
		core::Distance o2prime_instantiation_distance_cutoff_;
	};

} //rigid_body
} //legacy
} //rna
} //sampling
} //stepwise
} //protocols

#endif
