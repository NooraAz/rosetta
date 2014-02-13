// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/FragmentCM.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/movers/FragmentCM.hh>

// Package headers
#include <core/environment/DofPassport.hh>
#include <core/environment/LocalPosition.hh>
#include <protocols/environment/DofUnlock.hh>

#include <protocols/environment/claims/TorsionClaim.hh>
#include <protocols/environment/claims/CutBiasClaim.hh>

// Project headers
#include <core/kinematics/MoveMap.hh>

#include <core/fragment/FragSet.hh>

#include <core/pose/Pose.hh>

#include <protocols/simple_moves/FragmentMover.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

static basic::Tracer tr("protocols.environment.movers.FragmentCM", basic::t_info);

namespace protocols {
namespace environment {

using namespace core::environment;

FragmentCM::FragmentCM():
  ClaimingMover(),
  mover_( NULL ),
  label_( "BASE" )
{}

FragmentCM::FragmentCM( simple_moves::FragmentMoverOP mover,
                        std::string const& label ):
  ClaimingMover(),
  label_( label )
{
  set_mover( mover );
}

FragmentCM::~FragmentCM() {}

void FragmentCM::set_label( std::string const& label ){
  label_ = label;
}

void FragmentCM::set_mover( simple_moves::FragmentMoverOP mover ){
  mover_ = mover;
  type( get_name() );
}


claims::EnvClaims FragmentCM::yield_claims( core::pose::Pose& ){
  using namespace claims;
  claims::EnvClaims claim_list;

  TorsionClaimOP new_claim = new TorsionClaim( this, label(),
                                               std::make_pair( mover()->fragments()->min_pos(),
                                                               mover()->fragments()->max_pos() ) );

  new_claim->init_strength( DOES_NOT_INITIALIZE );
  new_claim->ctrl_strength( CAN_CONTROL );
  claim_list.push_back( new_claim );


  core::fragment::SecondaryStructureOP ss = new core::fragment::SecondaryStructure( *mover()->fragments() );
  CutBiasClaimOP cb_claim = new CutBiasClaim( this, label(), *ss );
  claim_list.push_back( cb_claim );

  return claim_list;
}

void FragmentCM::initialize( Pose& pose ){
  DofUnlock activation( pose.conformation(), passport() );
  mover()->apply_at_all_positions( pose );
}

void FragmentCM::apply( Pose& pose ) {
  if( pose.conformation().is_protected() ) {
    DofUnlock activation( pose.conformation(), passport() );
    mover()->apply( pose );
  } else {
    mover()->apply( pose );
  }
}

void FragmentCM::passport_updated() {
  //TODO: account for possible offset shifts.
  if( has_passport() ){
    mover()->set_movemap( passport()->render_movemap() );
  } else {
    mover()->set_movemap( new core::kinematics::MoveMap() );
  }
}

std::string FragmentCM::get_name() const {
  return "FragmentCM("+mover()->get_name()+")";
}


} // environment
} // protocols
