// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/StarTreeBuilder.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_NONLOCAL_STARTREEBUILDER_HH_
#define PROTOCOLS_NONLOCAL_STARTREEBUILDER_HH_

// Unit headers
#include <protocols/nonlocal/StarTreeBuilder.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/loops/Loop.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>

// Package headers
#include <protocols/nonlocal/TreeBuilder.hh>

namespace protocols {
namespace nonlocal {

class StarTreeBuilder : public TreeBuilder {
 public:
  StarTreeBuilder();

  /// @brief Constructs a star fold tree by placing a virtual residue at
  /// <chunks> center of mass and adding jumps from it to a stochastically
  /// chosen residue in each chunk.
  void set_up(const protocols::loops::Loops& chunks, core::pose::Pose* pose);

  /// @brief Removes the virtual residue added to <pose> in calls to build()
  void tear_down(core::pose::Pose* pose);

 private:
  /// @brief Stochastically selects an anchor position on [fragment.start(), fragment.stop()]
  /// according to per-residue structural conservation
  core::Size choose_conserved_position(const protocols::loops::Loop& chunk,
                                       const core::pose::Pose& pose) const;

  /// @brief Index of the virtual residue we added to the pose in build()
  int virtual_res_;
};

}  // namespace nonlocal
}  // namespace protocols

#endif  // PROTOCOLS_NONLOCAL_STARTREEBUILDER_HH_
