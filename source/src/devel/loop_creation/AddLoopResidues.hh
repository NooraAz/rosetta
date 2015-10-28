// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file AddLoopResidues.hh
///
/// @brief A Mover that uses loophash to find fragments that can
/// bridge a gap with minimal modifications to the original pose.
/// @author Tim Jacobs


#ifndef INCLUDED_devel_loop_creation_AddLoopResidues_HH
#define INCLUDED_devel_loop_creation_AddLoopResidues_HH

//Unit
#include <devel/loop_creation/AddLoopResidues.fwd.hh>
#include <protocols/moves/Mover.hh>

//Core
#include <core/pose/Pose.hh>

#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>

//Utility
#include <utility/vector1.hh>

//C++
#include <set>

namespace devel {
namespace loop_creation {

class AddLoopResidues : public protocols::moves::Mover {
	
public:

	AddLoopResidues();
	
	protocols::moves::MoverOP
	clone() const;

	protocols::moves::MoverOP
	fresh_instance() const;
	
	std::string
	get_name() const;

	///@details we accumulate state, so we definitely need to reinitialize
	virtual bool
	reinitialize_for_new_input() const;

	virtual void
	apply( Pose & pose );
	
	void
	update_anchors(
		utility::vector1<core::Size> & loop_anchors,
		protocols::loops::Loop const & new_loop,
		core::Size index_of_new_loop
	);

	void
	dump_loops_file(
		std::string filename,
		protocols::loops::Loops loops
	);

	protocols::loops::Loops
	update_loops(
		protocols::loops::Loop const & new_loop,
		protocols::loops::Loops const & all_loops
	);

	void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const &filters,
		protocols::moves::Movers_map const &movers,
		Pose const & pose
	);
	
private:

	core::Size asym_size_;
	utility::vector1<core::Size> loop_sizes_;
	utility::vector1<core::Size> loop_anchors_;

};
	
} // loop_creation
} // devel

#endif