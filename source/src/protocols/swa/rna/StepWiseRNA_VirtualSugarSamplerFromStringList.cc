// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/rna/StepWiseRNA_VirtualSugarSamplerFromStringList.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu

#include <protocols/swa/rna/StepWiseRNA_VirtualSugarSamplerFromStringList.hh>
#include <protocols/swa/rna/StepWiseRNA_VirtualSugarSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_VirtualSugarUtil.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/rna/StepWiseRNA_OutputData.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/SugarModeling.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/conformation/Conformation.hh>

#include <fstream>
#include <ObjexxFCL/string.functions.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.swa.rna.StepWiseRNA_VirtualSugarSamplerFromStringList" );

namespace protocols {
namespace swa {
namespace rna {

	//Constructor
	StepWiseRNA_VirtualSugarSamplerFromStringList::StepWiseRNA_VirtualSugarSamplerFromStringList( StepWiseRNA_JobParametersCOP & job_parameters,
																																																utility::vector1< std::string > const sample_virtual_sugar_string_list):
		job_parameters_( job_parameters ),
		sample_virtual_sugar_string_list_( sample_virtual_sugar_string_list ),
		silent_file_out_( "default.out" ),
		use_phenix_geo_( false ),
		legacy_mode_( false ),
		integration_test_mode_( false ),
		tag_( "" )
	{
	}

	//Destructor
	StepWiseRNA_VirtualSugarSamplerFromStringList::~StepWiseRNA_VirtualSugarSamplerFromStringList()
	{}

	/////////////////////
	std::string
	StepWiseRNA_VirtualSugarSamplerFromStringList::get_name() const {
		return "StepWiseRNA_VirtualSugarSamplerFromStringList";
}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_VirtualSugarSamplerFromStringList::apply( core::pose::Pose & pose ){

		using namespace ObjexxFCL;
		using namespace core::io::silent;
		using namespace core::id;
		using namespace core::scoring;
		using namespace core::conformation;

		output_title_text( "Enter StepWiseRNA_VirtualSugarSampler::sample_virtual_sugar", TR );

		// why all this rigamarole? -- rhiju
		ConformationOP copy_conformation = pose.conformation().clone();
		pose::Pose new_pose;
		new_pose.set_new_conformation( copy_conformation );
		pose = new_pose;
		tag_into_pose( pose, "" );

		///////////////////////////////////////////////////////////////
		pose::Pose const pose_save = pose;
		pose = pose_save; //this recopy is useful for triggering graphics.

		utility::vector1< SugarModeling > sugar_modeling_list = setup_sugar_modeling_list( pose );
		if ( empty_sugar_modeling_list( sugar_modeling_list ) ) return;

		for ( Size n = 1; n <= sugar_modeling_list.size(); n++ ){

			TR << TR.Blue << "Sampling sugar " << n <<  " out of " << sugar_modeling_list.size() << TR.Reset << std::endl;
			SugarModeling & curr_modeling = sugar_modeling_list[n];
			curr_modeling.set_base_and_pucker_state( pose, job_parameters_ );

			StepWiseRNA_VirtualSugarSampler virtual_sugar_sampler( job_parameters_, curr_modeling );
			virtual_sugar_sampler.set_tag( "VIRT_RIBOSE_NUM_" + string_of( n ) );
			virtual_sugar_sampler.set_scorefxn( scorefxn_ );
			virtual_sugar_sampler.set_integration_test_mode( integration_test_mode_ );
			virtual_sugar_sampler.set_use_phenix_geo( use_phenix_geo_ );
			virtual_sugar_sampler.set_legacy_mode( legacy_mode_ );
			virtual_sugar_sampler.set_virtual_sugar_is_from_prior_step( false );

			virtual_sugar_sampler.apply( pose );

			std::sort( curr_modeling.pose_list.begin(), curr_modeling.pose_list.end(), sort_pose_by_score );
			if ( empty_pose_data_list( curr_modeling.pose_list, n, "curr_modeling" ) ) return;
		}

		///////////////////////////////////////////////////////////////////////////////////
		utility::vector1< PoseOP > pose_data_list;
		enumerate_starting_pose_data_list( pose_data_list, sugar_modeling_list, pose );
		minimize_all_sampled_floating_bases( pose, sugar_modeling_list, pose_data_list, scorefxn_, job_parameters_, false /*virtualize_other_partition*/ );

		output_pose_data_list( pose_data_list );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< SugarModeling >
	StepWiseRNA_VirtualSugarSamplerFromStringList::setup_sugar_modeling_list( pose::Pose const & pose ) const {

		using namespace ObjexxFCL;

		utility::vector1< SugarModeling > sugar_modeling_list;

		std::string const & working_sequence = job_parameters_->working_sequence();

		for ( Size n = 1; n <= sample_virtual_sugar_string_list_.size(); n++ ){

			utility::vector1< std::string > const tokenize_list = tokenize( sample_virtual_sugar_string_list_[n], "-" );
			if ( tokenize_list.size() != 2 ) utility_exit_with_message( "tokenize_list != 2" );

			if ( tokenize_list[2] != "A" && tokenize_list[2] != "P" ){
				utility_exit_with_message( "tokenize_list[2] != \"A\" && tokenize_list[2] != \"P\" ( " + tokenize_list[2] + " )" );
			}

			bool const is_prepend = ( tokenize_list[2] == "P" );
			Size const full_sugar_res = string_to_int( tokenize_list[1] );
			Size const full_bulge_res = ( is_prepend ) ? full_sugar_res + 1 : full_sugar_res - 1;
			Size const full_ref_res   = ( is_prepend ) ? full_sugar_res + 2 : full_sugar_res - 2;

			TR.Debug << "Case: " << sample_virtual_sugar_string_list_[n];
			TR.Debug << " full_sugar_res = " << full_sugar_res << " full_bulge_res = " << full_bulge_res << " full_ref_res = " << full_ref_res;
			output_boolean( " is_prepend = ", is_prepend, TR.Debug );

			runtime_assert ( pose.total_residue() == working_sequence.size() );

			if ( check_is_working_res( full_sugar_res, job_parameters_ ) ){

				Size const working_sugar_res = check_validity_and_get_working_res( full_sugar_res, job_parameters_ );
				bool const sugar_is_virtual = pose.residue( working_sugar_res ).has_variant_type( "VIRTUAL_RIBOSE" );

				TR.Debug << " | working_sugar_res = " << working_sugar_res;
				output_boolean( " sugar_is_virtual = ", sugar_is_virtual, TR.Debug );

				if ( sugar_is_virtual ){

					Size const working_bulge_res = check_validity_and_get_working_res( full_bulge_res, job_parameters_ );
					Size const working_ref_res   = check_validity_and_get_working_res( full_ref_res, job_parameters_ );

					TR.Debug << " | working_bulge_res = " << working_bulge_res << " working_ref_res = " << working_ref_res;

					if ( is_prepend ){
						runtime_assert( working_sugar_res == ( working_bulge_res - 1 ) );
						runtime_assert( working_sugar_res == ( working_ref_res - 2 ) );
					} else{
						runtime_assert( working_sugar_res == ( working_bulge_res + 1 ) );
						runtime_assert( working_sugar_res == ( working_ref_res + 2 ) );
					}
					runtime_assert( pose.residue( working_bulge_res ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) );

					SugarModeling curr_modeling = SugarModeling();
					curr_modeling = SugarModeling( working_sugar_res, working_ref_res );
					sugar_modeling_list.push_back( curr_modeling );
				}

			} else{
				TR.Debug << " | full_sugar_res is not a working res! ";
			}
			TR.Debug << std::endl;
		}

		return sugar_modeling_list;
	}

	/////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_VirtualSugarSamplerFromStringList::empty_sugar_modeling_list( utility::vector1< SugarModeling > const &  sugar_modeling_list  ){

		TR.Debug << "num_virtual_sugar = " << sugar_modeling_list.size() << std::endl;
		if ( sugar_modeling_list.size() == 0 ){
			TR.Debug << "no_virtual_sugar ( sugar_modeling_list.size() == 0 ). EARLY RETURN/NO OUTPUT SILENT_FILE!" << std::endl;

			std::ofstream outfile;
			outfile.open( silent_file_out_.c_str() ); //Opening the file with this command removes all prior content..
			outfile << "no_virtual_sugar ( sugar_modeling_list.size() == 0 ).\n";
			outfile.flush();
			outfile.close();

			return true;
		}

		return false;
	}

	/////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_VirtualSugarSamplerFromStringList::empty_pose_data_list( utility::vector1< pose::PoseOP > const & pose_list, Size const n, std::string tag ){

	if (  pose_list.size() == 0 ){
		TR.Debug << "Case n = " << n << " is_sugar_virt == True but " << tag << ".pose_list.size() == 0. EARLY RETURN!" << std::endl;

		std::ofstream outfile;
		outfile.open( silent_file_out_.c_str() ); //Opening the file with this command removes all prior content..
		outfile << "num_virtual_sugar != 0 but for one of the sampled virtual_sugar," << tag << ".pose_list.size() == 0.\n";
		outfile.flush();
		outfile.close();
		return true;
	}

	return false;
}

	/////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_VirtualSugarSamplerFromStringList::output_pose_data_list( utility::vector1< pose::PoseOP > & pose_data_list ){

		for ( Size n = 1; n <= pose_data_list.size(); n++ ){
			Pose & pose = ( *pose_data_list[n] ); //set viewer_pose;
			std::string pose_tag = tag_ + "_sample_sugar" + tag_from_pose( *pose_data_list[n]);
			if ( job_parameters_->gap_size() == 0 ) utility_exit_with_message( "job_parameters_->gap_size() == 0" );
			( *scorefxn_ )( pose );

			output_data( silent_file_out_, pose_tag, false, pose, job_parameters_->working_native_pose(), job_parameters_ );
		}
	}

	//////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_VirtualSugarSamplerFromStringList::set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
		scorefxn_ = scorefxn;
	}

} //rna
} //swa
} //protocols
