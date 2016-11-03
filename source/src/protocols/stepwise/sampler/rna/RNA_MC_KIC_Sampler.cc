// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/rna/RNA_MC_KIC_Sampler.cc
/// @brief Sample and torsions and close an RNA loop.
/// @author Fang-Chieh Chou


// Unit headers
#include <protocols/stepwise/sampler/rna/RNA_MC_KIC_Sampler.hh>

// Package headers
#include <protocols/stepwise/sampler/StepWiseSamplerSizedComb.hh>
#include <protocols/stepwise/sampler/MC_Comb.hh>
#include <protocols/stepwise/sampler/rna/RNA_KinematicCloser_DB.hh>
#include <protocols/stepwise/sampler/rna/RNA_ChiStepWiseSampler.hh>
#include <protocols/stepwise/sampler/rna/RNA_SugarStepWiseSampler.hh>
#include <protocols/stepwise/sampler/screener/RNA_TorsionScreener.hh>
#include <protocols/stepwise/sampler/MC_OneTorsion.hh>
#include <protocols/stepwise/sampler/rna/RNA_MC_Sugar.hh>
#include <protocols/stepwise/monte_carlo/mover/TransientCutpointHandler.hh>
#include <numeric/random/WeightedSampler.hh>
#include <numeric/random/random.hh>
#include <core/id/NamedAtomID.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/chemical/rna/RNA_SamplerUtil.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/rna/RNA_SuiteName.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

using namespace core;
using namespace core::chemical::rna;
using namespace core::pose::rna;
using namespace protocols::stepwise::sampler::screener;
using namespace protocols::stepwise::monte_carlo::mover;
static THREAD_LOCAL basic::Tracer TR( "protocols.sampler.rna.RNA_MC_KIC_Sampler" );

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rna {

RNA_MC_KIC_Sampler::RNA_MC_KIC_Sampler(
	core::pose::PoseOP const & ref_pose,
	core::Size const moving_suite,
	core::Size const chainbreak_suite,
	bool const change_ft
):
	MC_StepWiseSampler(),
	moving_suite_( moving_suite ),
	chainbreak_suite_( chainbreak_suite ),
	init_pucker_( NORTH),
	pucker_flip_rate_( 0.1 ),
	base_state_( ANY_CHI ), // ANY_CHI, ANTI, SYN, NO_CHI
	sample_nucleoside_( moving_suite + 1 ), // default, may be replaced.
	max_tries_( 100 ),
	verbose_( false ),
	extra_epsilon_( false ),
	extra_chi_( false ),
	skip_same_pucker_( false ),
	idealize_coord_( false ),
	torsion_screen_( true ),
	random_chain_closed_( true ),
	did_close( false ),
	screener_( screener::RNA_TorsionScreenerOP( new RNA_TorsionScreener ) ),
	cutpoint_handler_(TransientCutpointHandlerOP ( new TransientCutpointHandler( moving_suite, chainbreak_suite, change_ft )))
{
	ref_pose_.reset(new core::pose::Pose(*ref_pose));
}

//////////////////////////////////////////////////////////////////////////
void RNA_MC_KIC_Sampler::init() {
	using namespace core::id;

	/////// Backbone MC Sampler ////////

	TorsionID epsilon_ID (TorsionID( moving_suite_, BB, EPSILON));
	TorsionIDs.push_back(epsilon_ID);
	MC_OneTorsionOP epsilon_torsion( new MC_OneTorsion( epsilon_ID, ref_pose_->torsion( epsilon_ID) ));
	bb_samplers_.push_back(epsilon_torsion);

	TorsionID zeta_ID (TorsionID( moving_suite_, BB, ZETA));
	TorsionIDs.push_back(zeta_ID);
	MC_OneTorsionOP zeta_torsion( new MC_OneTorsion( zeta_ID, ref_pose_->torsion( zeta_ID) ));
	bb_samplers_.push_back(zeta_torsion);

	TorsionID alpha1_ID (TorsionID( moving_suite_ + 1, BB, ALPHA));
	TorsionIDs.push_back(alpha1_ID);
	MC_OneTorsionOP alpha1_torsion( new MC_OneTorsion( alpha1_ID, ref_pose_->torsion( alpha1_ID) ));
	bb_samplers_.push_back(alpha1_torsion);

	TorsionID alpha2_ID (TorsionID( chainbreak_suite_ + 1, BB, ALPHA));
	TorsionIDs.push_back(alpha2_ID);
	MC_OneTorsionOP alpha2_torsion( new MC_OneTorsion( alpha2_ID, ref_pose_->torsion( alpha2_ID) ));
	bb_samplers_.push_back(alpha2_torsion);

	//Add all the pivot torsions to the list of the torsion IDs so that we have a
	//full list of the torsions that could be modified by this sampler

	TorsionID pivot1 (TorsionID( moving_suite_ + 1, BB, BETA ));
	TorsionIDs.push_back(pivot1);
	TorsionID pivot2 (TorsionID( moving_suite_ + 1, BB, GAMMA ));
	TorsionIDs.push_back(pivot2);
	TorsionID pivot3 (TorsionID( moving_suite_ + 1, BB, EPSILON ));
	TorsionIDs.push_back(pivot3);
	TorsionID pivot4 (TorsionID( moving_suite_ + 1, BB, ZETA ));
	TorsionIDs.push_back(pivot4);
	TorsionID pivot5 (TorsionID( chainbreak_suite_ + 1, BB, BETA ));
	TorsionIDs.push_back(pivot5);
	TorsionID pivot6 (TorsionID( chainbreak_suite_ + 1, BB, GAMMA ));
	TorsionIDs.push_back(pivot6);

	// Now find the initial torsions
	core::Real torsion;
	for ( Size i = 1; i <= TorsionIDs.size(); ++i ) {
		torsion = ref_pose_->torsion(TorsionIDs[i]);
		initial_torsions_.push_back( torsion );
		angle_min_.push_back(-180);
		angle_max_.push_back(180);
	}


	//Let's not worry about the sugar rotamers for now (1/15/15)

	for ( Size i = 1; i <= bb_samplers_.size(); ++i ) {
		bb_samplers_[i]->init();
	}

	////////// Make a stored loop closer /////////
	stored_loop_closer_ = RNA_KinematicCloser_DBOP( new RNA_KinematicCloser_DB(
		ref_pose_, moving_suite_, chainbreak_suite_ ) );

	////////// Loop Closer //////////
	loop_closer_ = RNA_KinematicCloser_DBOP( new RNA_KinematicCloser_DB(
		ref_pose_, moving_suite_, chainbreak_suite_ ) );
	loop_closer_->set_verbose( verbose_ );

	set_init( true );
}
//////////////////////////////////////////////////////////////////////////
void RNA_MC_KIC_Sampler::reset( pose::Pose & /*pose*/ ) {
	runtime_assert( is_init() );
	// //This resets the driver torsions
	// for ( Size i = 1; i <= bb_samplers_.size(); ++i ) {
	//  bb_samplers_[i]->reset();
	//  bb_samplers_[i]->apply( pose );
	// }
}
////////////////////////////////////////////////////////////////////////
void RNA_MC_KIC_Sampler::operator++() {
	return;
}
////////////////////////////////////////////////////////////////////////
void RNA_MC_KIC_Sampler::next( const Pose & pose ) {
	runtime_assert( not_end() );
	Pose pose_copy = pose;
	//Need to update the stored angles in the samplers in case they were changed in a different sampler
	for ( Size i = 1; i <= bb_samplers_.size(); ++i ) {
		bb_samplers_[i]->set_angle( pose_copy.torsion( TorsionIDs[i] ) ) ;
	}

	//First need to make a cutpoint
	//Then find the solutions for the closed pose
	cutpoint_handler_->put_in_cutpoints( pose_copy );
	stored_loop_closer_->init(pose_copy, pose);
	stored_jacobians_ = stored_loop_closer_->get_all_jacobians();
	for ( Size i = 1; i <= max_tries_; ++i ) {
		for ( Size i = 1; i <= bb_samplers_.size(); ++i ) {
			++( *bb_samplers_[i]);
			bb_samplers_[i]->apply( pose_copy );
		}

		loop_closer_->init(pose_copy, pose);
		//if there are 0 solutions, try again
		if ( !loop_closer_->not_end() ) {
			continue;
		}
		cutpoint_handler_->take_out_cutpoints( pose_copy );
		did_close =  true ;
		return;
	}

	cutpoint_handler_->take_out_cutpoints( pose_copy );
	did_close = false ;

	TR.Debug << "Chain not closable after " << max_tries_ << " tries!" << std::endl;
}
///////////////////////////////////////////////////////////////////////////
bool RNA_MC_KIC_Sampler::check_closed() const {
	return did_close;
}
///////////////////////////////////////////////////////////////////////////
bool RNA_MC_KIC_Sampler::check_moved() const {
	return did_move;
}
///////////////////////////////////////////////////////////////////////////
bool RNA_MC_KIC_Sampler::not_end() const {
	//since it's a MC sampler, not enumerative, never ends
	runtime_assert( is_init() );
	return true;
}
///////////////////////////////////////////////////////////////////////////
void RNA_MC_KIC_Sampler::apply() {
	apply( *ref_pose_ );
}
///////////////////////////////////////////////////////////////////////////
void RNA_MC_KIC_Sampler::apply( pose::Pose & pose_in ) {
	runtime_assert( is_init() );
	Pose pose = pose_in;
	cutpoint_handler_->put_in_cutpoints( pose );
	used_current_solution_ = false;
	did_move = false;

	if ( did_close ) {
		current_jacobians_ = loop_closer_->get_all_jacobians();
		all_jacobians_ = current_jacobians_;
		copy_stored_jacobians_ = stored_jacobians_;
		all_jacobians_.insert(all_jacobians_.end(),
			copy_stored_jacobians_.begin(), copy_stored_jacobians_.end());
		jacobian_sampler_ = numeric::random::WeightedSampler( all_jacobians_ );
		solution_ = jacobian_sampler_.random_sample( numeric::random::rg().uniform() );
		if ( solution_ <= current_jacobians_.size() ) {
			for ( Size i = 1; i <= bb_samplers_.size(); ++i ) {
				bb_samplers_[i]->apply( pose );
			}
			loop_closer_->apply( pose, solution_ );
			used_current_solution_ = true;
			did_move = true;
		} else {
			Real calculated_jacobian = get_jacobian( pose );
			stored_loop_closer_->apply( pose, (solution_ - current_jacobians_.size()) );
			//Check whether this is actually a new pose
			//This isn't really a good way to check if the pose is the same...change this!!
			Real picked_jacobian = stored_jacobians_[ solution_ - current_jacobians_.size()];
			if ( std::abs(calculated_jacobian - picked_jacobian) > 0.0000000001 ) { did_move = true; }
			//used_current solution remains false, so don't update stored angles, don't update stored loop closer
		}
	} else {
		all_jacobians_ = stored_jacobians_;
		jacobian_sampler_ = numeric::random::WeightedSampler( all_jacobians_ );
		solution_ = jacobian_sampler_.random_sample( numeric::random::rg().uniform() );
		Real calculated_jacobian = get_jacobian( pose );
		stored_loop_closer_->apply( pose, solution_ );
		//Check whether this is actually a new pose
		//This isn't really a good way to check if the pose is the same...change this!!
		Real picked_jacobian = stored_jacobians_[ solution_ - current_jacobians_.size()];
		if ( std::abs(calculated_jacobian - picked_jacobian) > 0.0000000001 ) { did_move = true; }
		//used_current solution remains false, so don't update stored angles, don't update stored loop closer
	}
	cutpoint_handler_->take_out_cutpoints( pose );
	if ( did_move ) {
		if ( !( check_angles_in_range( pose ) ) ) {
			did_move = false;
			return;
		}
	}
	pose_in = pose;
}
///////////////////////////////////////////////////////////////////////////
Real RNA_MC_KIC_Sampler::get_jacobian( pose::Pose & pose ) {
	jacobian_ = loop_closer_->get_jacobian(pose);
	return jacobian_;
}
///////////////////////////////////////////////////////////////////////////
Real RNA_MC_KIC_Sampler::vector_sum( utility::vector1< core::Real > const & vector ) {
	sum_ = 0;
	for ( Size i = 1; i<=vector.size(); ++i ) {
		sum_ += vector[i];
	}
	return sum_;
}
///////////////////////////////////////////////////////////////////////////
void RNA_MC_KIC_Sampler::update() {
	return;
}
///////////////////////////////////////////////////////////////////////////
void RNA_MC_KIC_Sampler::update( const Pose & /*pose*/ ) {
	runtime_assert( is_init() );
	if ( used_current_solution_ ) {
		for ( Size i = 1; i <= bb_samplers_.size(); ++i ) {
			bb_samplers_[i]->update();
		}
	}
}
/////////////////////////////////////////////////////////////////////////
void RNA_MC_KIC_Sampler::set_gaussian_stdev( Real const setting ) {
	gaussian_stdev_ = setting;
	if ( is_init() ) {
		for ( Size i = 1; i <= bb_samplers_.size(); ++i ) {
			bb_samplers_[i]->set_gaussian_stdev( gaussian_stdev_ );
		}
	}
}
/////////////////////////////////////////////////////////////////////////
void RNA_MC_KIC_Sampler::set_angle_range_from_init_torsions( core::Real const range ) {
	runtime_assert( is_init() );
	for ( Size i = 1; i <= bb_samplers_.size(); ++i ) {
		bb_samplers_[i]->set_angle_range( initial_torsions_[i] - range, initial_torsions_[i] + range );
	}
	for ( Size i = bb_samplers_.size() + 1; i <= initial_torsions_.size(); ++i ) {
		angle_min_[i] = initial_torsions_[i] - range;
		angle_max_[i] = initial_torsions_[i] + range;
	}
}
/////////////////////////////////////////////////////////////////////////
bool RNA_MC_KIC_Sampler::check_angles_in_range( const Pose & pose ) {
	// Just need to check the torsions for the pivot angles, the ones for the driver angles
	// are already guaranteed to be in the correct range by the set_angle_range fn
	// (see set_angle_range_from_init_torsions)
	for ( Size i = bb_samplers_.size() + 1; i <= TorsionIDs.size(); i++ ) {
		if ( angle_max_[i] - angle_min_[i] >= 360 ) continue;
		angle_ = pose.torsion( TorsionIDs[i] );
		while ( angle_ < angle_min_[i] ) angle_ += 360;
		while ( angle_ >= angle_max_[i] ) angle_ -= 360;
		if ( angle_ >= angle_min_[i] ) continue;
		return false;
	}
	return true;
}
/////////////////////////////////////////////////////////////////////////
/// @brief Name of the class
std::string
RNA_MC_KIC_Sampler::get_name() const {
	return "RNA_MC_KIC_Sampler moving_suite:" + utility::to_string(moving_suite_) + " chainbreak_suite:" +
		utility::to_string(chainbreak_suite_);
}

} //rna
} //sampler
} //stepwise
} //protocols
