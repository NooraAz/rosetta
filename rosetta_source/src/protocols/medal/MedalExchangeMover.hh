// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/medal/MedalExchangeMover.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef PROTOCOLS_MEDAL_MEDAL_EXCHANGE_MOVER_HH_
#define PROTOCOLS_MEDAL_MEDAL_EXCHANGE_MOVER_HH_

// Unit header
#include <protocols/medal/MedalExchangeMover.fwd.hh>

// C/C++ headers
#include <string>

// Utility headers
#include <utility/vector1.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace medal {

class MedalExchangeMover : public protocols::moves::Mover {
 public:
  MedalExchangeMover();
  void apply(core::pose::Pose& pose);

  /// @brief Uses the copy constructor to create a new instance
  protocols::moves::MoverOP clone() const;

  /// @brief Uses the no-argument constructor to create a new instance
  protocols::moves::MoverOP fresh_instance() const;

  /// @brief Returns the name of this mover
  std::string get_name() const;

 private:
  /// @brief Creates a coordinate constraint between a fixed virtual atom and the
  /// CA atom of each aligned residue. The constraint is evaluated using a BoundFunc
  /// functional form, with lower bound set to the initial distance - delta and
  /// upper bound set to the initial distance + delta. The value of delta, as well
  /// as the standard deviation, are controllable via the command line.
  void setup_constraints(const core::pose::Pose& pose,
                         protocols::loops::LoopsCOP aligned,
                         core::scoring::constraints::ConstraintSetOP constraints) const;

  /// @brief Computes the probability of selecting each residue as a candidate
  /// for fragment insertion. P(unaligned) = 0, P(aligned) > 0.
  void setup_sampling_probs(core::Size num_residues,
                            const core::kinematics::FoldTree& tree,
                            protocols::loops::LoopsCOP aligned,
                            utility::vector1<double>* probs) const;

  core::fragment::FragSetOP fragments_;
  core::fragment::SecondaryStructureOP pred_ss_;
};

}  // namespace medal
}  // namespace protocols

#endif  // PROTOCOLS_MEDAL_MEDAL_EXCHANGE_MOVER_HH_
