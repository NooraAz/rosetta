// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/simple_moves/RepulsiveOnlyMover.cc
/// @brief calculate repulsive energy only for certain residues
/// @author Will Sheffler
/// @author Ray Wang (wangyr@uw.edu)

// Unit Headers
#include <protocols/simple_moves/RepulsiveOnlyMover.hh>

// Package headers 

// Project headers 
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// tracer
#include <basic/Tracer.hh>

// Utility Headers

// C++ Headers
#include <iostream>

#include <utility/vector1.hh>


// ObjexxFCL Headers

namespace protocols {
namespace simple_moves {

static basic::Tracer TR( "protocol.mover.RepulsiveOnlyMover" );


/// RepulsiveOnlyMover; based on the protocols::moves::Mover basis class
RepulsiveOnlyMover::RepulsiveOnlyMover() : protocols::moves::Mover() {}
RepulsiveOnlyMover::~RepulsiveOnlyMover() {}


// @brief apply function here
void
RepulsiveOnlyMover::apply( core::pose::Pose & pose ) {
	if (	basic::options::option[basic::options::OptionKeys::in::replonly_loops]() ){
		for ( core::Size i=1; i<=pose.n_residue(); i++ ){   
			if ( pose.secstruct(i)=='L' ){
				if ( ! pose.residue(i).has_variant_type( core::chemical::REPLONLY ) ){
					core::pose::remove_lower_terminus_type_from_pose_residue( pose, i );
					core::pose::remove_upper_terminus_type_from_pose_residue( pose, i );
					core::pose::add_variant_type_to_pose_residue( pose, "REPLONLY", i );
				}
			}else {
				if ( pose.residue(i).has_variant_type( core::chemical::REPLONLY )){
					core::pose::remove_variant_type_from_pose_residue( pose, "REPLONLY", i );
				}
			}
		}
	}
	if ( basic::options::option[ basic::options::OptionKeys::in::replonly_residues ].user() ){
		utility::vector1<Size> replonly_rsd = basic::options::option[ basic::options::OptionKeys::in::replonly_residues ]();
		//TR << "RepulsiveOnly protocols::moves::Mover has been called" << std::endl;
		for ( core::Size i=1; i<=replonly_rsd.size(); i++ ){   
			core::chemical::ResidueType const & gly( pose.residue(i).residue_type_set().name_map("GLY") );
			core::pose::replace_pose_residue_copying_existing_coordinates( pose, replonly_rsd[i], gly );
		//TR << replonly_rsd[i] << " has been changed as GLY" << std::endl;
			if ( ! pose.residue( replonly_rsd[i] ).has_variant_type( core::chemical::REPLONLY ) ){
				core::pose::remove_lower_terminus_type_from_pose_residue( pose, replonly_rsd[i] );
				core::pose::remove_upper_terminus_type_from_pose_residue( pose, replonly_rsd[i] );
				core::pose::add_variant_type_to_pose_residue( pose, "REPLONLY", replonly_rsd[i] );
			} 
		}   
	}   
}   

std::string
RepulsiveOnlyMover::get_name() const {
	return "RepulsiveOnlyMover";
}


} // moves
} // protocols

