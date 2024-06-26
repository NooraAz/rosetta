// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file MinMover.cc
/// @brief
/// @author Ingemar Andre

// Unit headers
#include <protocols/symmetric_docking/SymFoldandDockMoveRbJumpMover.hh>
#include <protocols/symmetry/SetupForSymmetryMover.hh>

// Package headers
#include <core/pose/symmetry/util.hh>

#include <protocols/symmetric_docking/SymFoldandDockCreators.hh>

#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

// Utility Headers
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


namespace protocols {
namespace symmetric_docking {

static basic::Tracer TR( "protocols.moves.symmetry.SymFoldandDockMoveRbJumpMover" );


SymFoldandDockMoveRbJumpMover::SymFoldandDockMoveRbJumpMover()
: Mover("SymFoldandDockMoveRbJumpMover")
{}

void
SymFoldandDockMoveRbJumpMover::apply( core::pose::Pose & pose )
{
	using namespace core::conformation::symmetry;
	protocols::symmetry::SetupForSymmetryMover setup;
	setup.apply( pose );
	core::pose::symmetry::find_new_symmetric_jump_residues( pose );
}


void
SymFoldandDockMoveRbJumpMover::parse_my_tag(
	utility::tag::TagCOP const,
	basic::datacache::DataMap &
)
{
}





std::string SymFoldandDockMoveRbJumpMover::get_name() const {
	return mover_name();
}

std::string SymFoldandDockMoveRbJumpMover::mover_name() {
	return "SymFoldandDockMoveRbJumpMover";
}

void SymFoldandDockMoveRbJumpMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "This seems to help set up Symmetric Fold and Dock. It symmetrizes the asymmetric input pose and then finds the jump residues between rigid bodies.", attlist );
}

std::string SymFoldandDockMoveRbJumpMoverCreator::keyname() const {
	return SymFoldandDockMoveRbJumpMover::mover_name();
}

protocols::moves::MoverOP
SymFoldandDockMoveRbJumpMoverCreator::create_mover() const {
	return utility::pointer::make_shared< SymFoldandDockMoveRbJumpMover >();
}

void SymFoldandDockMoveRbJumpMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SymFoldandDockMoveRbJumpMover::provide_xml_schema( xsd );
}


} // symmetric_docking
} // protocols
