// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Monica Berrondo
/// @author Modified by Sergey Lyskov

// Unit headers
#include <protocols/moves/RepackSidechainsMover.hh>
#include <protocols/moves/RepackSidechainsMoverCreator.hh>
// AUTO-REMOVED #include <protocols/moves/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/pack_missing_sidechains.hh>

#include <core/id/AtomID_Mask.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
// AUTO-REMOVED #include <utility/string_util.hh> // string_split

// option key includes
// AUTO-REMOVED #include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/options/BooleanVectorOption.hh>
#include <basic/Tracer.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

namespace protocols {
namespace moves {

using basic::Warning;
using basic::t_warning;
static basic::Tracer TR("protocols.moves.RepackSidechainsMover");

/// RepackSidechainsMover

std::string
RepackSidechainsMoverCreator::keyname() const
{
	return RepackSidechainsMoverCreator::mover_name();
}

protocols::moves::MoverOP
RepackSidechainsMoverCreator::create_mover() const {
	return new RepackSidechainsMover;
}

std::string
RepackSidechainsMoverCreator::mover_name()
{
	return "RepackSidechainsMover";
}

RepackSidechainsMover::RepackSidechainsMover() :
 	Mover("RepackSidechainsMover"),
	scorefxn_(0)
{}

// RepackSidechainsMover::RepackSidechainsMover( std::string const & type_name ) :
// 	Mover( type_name ),
// 	scorefxn_(0)
// {}

	// constructors with arguments
RepackSidechainsMover::RepackSidechainsMover(
	ScoreFunctionCOP scorefxn
) :
	Mover("RepackSidechainsMover"),
	scorefxn_( scorefxn )
{}

RepackSidechainsMover::RepackSidechainsMover( RepackSidechainsMover const & other ) :
	//utility::pointer::ReferenceCount(),
	Mover( other )
{
	scorefxn_ = other.scorefxn();
}

void
RepackSidechainsMover::apply( Pose & pose )
{
// repack missing sidechains
	core::id::AtomID_Mask all_atoms( true );
	core::pose::initialize_atomid_map( all_atoms, pose );

	//build a PackerTask to control rotamer_trials
	core::pack::task::PackerTaskOP task = core::pack::task::TaskFactory::create_packer_task( pose );
	task->initialize_from_command_line();
	task->restrict_to_repacking();

	utility::vector1_bool repackable;
	bool something_to_pack = core::pack::figure_out_repackable_residues( pose, all_atoms, repackable );
	if (!something_to_pack) return;

	//task is set up
	task->restrict_to_residues(repackable);
	core::pack::pack_rotamers( pose, *scorefxn_, task );

}

std::string
RepackSidechainsMover::get_name() const {
	return RepackSidechainsMoverCreator::mover_name();
}


///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
RepackSidechainsMover::parse_my_tag(
	TagPtr const tag,
	DataMap & datamap,
	Filters_map const & filters,
	Movers_map const & movers,
	Pose const & pose
)
{
	parse_score_function( tag, datamap, filters, movers, pose );
	//	parse_task_operations( tag, datamap, filters, movers, pose );
}

///@brief parse "scorefxn" XML option (can be employed virtually by derived Packing movers)
void
RepackSidechainsMover::parse_score_function(
	TagPtr const tag,
	DataMap const & datamap,
	Filters_map const &,
	Movers_map const &,
	Pose const &
)
{
	core::scoring::ScoreFunctionOP new_score_function( protocols::rosetta_scripts::parse_score_function( tag, datamap ) );
	if ( new_score_function == 0 ) return;
	set_scorefxn( new_score_function );
}

///@brief required in the context of the parser/scripting scheme
MoverOP
RepackSidechainsMover::fresh_instance() const
{
	return new RepackSidechainsMover;
}

///@brief required in the context of the parser/scripting scheme
MoverOP
RepackSidechainsMover::clone() const
{
	return new RepackSidechainsMover( *this );
}

// setters
void RepackSidechainsMover::set_scorefxn( ScoreFunctionCOP sf )
{
	scorefxn_ = sf;
}

} // moves
} // protocols

