// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/StarTreeBuilder.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <protocols/nonlocal/StarTreeBuilder.hh>

// C/C++ headers
#include <utility>

// Objexx headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

// Utility headers
#include <numeric/xyzVector.hh>
#include <numeric/random/WeightedReservoirSampler.hh>
#include <utility/vector1.hh>

// Project headers
#include <core/types.hh>
#include <core/fragment/SecondaryStructure.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/datacache/StructuralConservationStore.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Package headers
#include <protocols/nonlocal/CutFinder.hh>

namespace protocols {
namespace nonlocal {

StarTreeBuilder::StarTreeBuilder() : virtual_res_(-1) {}

/// Note: assumes <chunks> are sorted in increasing order of start position
void StarTreeBuilder::set_up(const protocols::loops::Loops& chunks, core::pose::Pose* pose) {
  using core::Size;
  using core::Real;
  using core::fragment::SecondaryStructure;
  using core::fragment::SecondaryStructureOP;
  using protocols::loops::Loop;
  using protocols::loops::Loops;
  using utility::vector1;

  assert(pose);
  assert(chunks.num_loop());

  // Number of residues before addition of virtual residue
  Size num_residues = pose->total_residue();

  // Add a virtual residue at the center of mass (updates the fold tree)
  numeric::xyzVector<Real> center;
  chunks.center_of_mass(*pose, &center);
  core::pose::addVirtualResAsRoot(center, *pose);

  // Initialize member variable <virtual_res_> with the index of the newly added
  // virtual residue. Subsequently, <virtual_res_> can serve as a proxy for
  // <num_residues>
  virtual_res_ = pose->total_residue();
  bool has_conservation = pose->data().has(core::pose::datacache::CacheableDataType::STRUCTURAL_CONSERVATION);

  vector1<std::pair<int, int> > jumps;
  for (Loops::const_iterator i = chunks.begin(); i != chunks.end(); ++i) {
    const Loop& chunk = *i;

    // When available, use structural conservation to inform jump selection.
    // Otherwise, choose the midpoint of the region.
    Size anchor_position = has_conservation ?
        choose_conserved_position(chunk, *pose) :
        (chunk.start() + chunk.stop()) / 2;

    // virtual residue => anchor position
    jumps.push_back(std::make_pair(virtual_res_, anchor_position));
  }

  // TODO(cmiles)
  //  - Calculate RMSD over jump residues (jump residue + 1 residue on each side)
  //  - Comment in silent file
  //  - Tracer output

  // Insert cutpoints between adjacent jumps
  SecondaryStructureOP ss = new SecondaryStructure(*pose);
  vector1<int> cuts;
  for (Size i = 2; i <= jumps.size(); ++i) {
    Size cutpoint = CutFinder::choose_cutpoint(jumps[i-1].second + 1, jumps[i].second - 1, ss);
    cuts.push_back(cutpoint);
  }

  // Remember to include the original cutpoint at the end of the chain
  // (before the virtual residue)
  cuts.push_back(num_residues);

  // Construct the star fold tree from the set of jumps and cuts above.
  // Reorder the resulting fold tree so that <virtual_res> is the root.
  ObjexxFCL::FArray2D_int ft_jumps(2, jumps.size());
  for (Size i = 1; i <= jumps.size(); ++i) {
    ft_jumps(1, i) = std::min(jumps[i].first, jumps[i].second);
    ft_jumps(2, i) = std::max(jumps[i].first, jumps[i].second);
  }

  ObjexxFCL::FArray1D_int ft_cuts(cuts.size());
  for (Size i = 1; i <= cuts.size(); ++i) {
    ft_cuts(i) = cuts[i];
  }

  // Try to build the star fold tree from jumps and cuts
  core::kinematics::FoldTree tree(pose->fold_tree());
  bool status = tree.tree_from_jumps_and_cuts(virtual_res_,   // nres_in
                                              jumps.size(),   // num_jump_in
                                              ft_jumps,       // jump_point
                                              ft_cuts,        // cuts
                                              virtual_res_);  // root
  if (!status)
    utility_exit_with_message("StarTreeBuilder: failed to build fold tree from cuts and jumps");

  // Update the pose's fold tree
  pose->fold_tree(tree);
}

core::Size StarTreeBuilder::choose_conserved_position(const protocols::loops::Loop& chunk,
                                                      const core::pose::Pose& pose) const {
  utility::vector1<core::Size> positions;
  utility::vector1<core::Real> weights;
  for (core::Size j = chunk.start(); j <= chunk.stop(); ++j) {
    positions.push_back(j);
    weights.push_back(core::pose::structural_conservation(pose, j));
  }

  // Perform weighted selection to choose (a single) jump position
  utility::vector1<core::Size> samples;
  numeric::random::WeightedReservoirSampler<Size> sampler(1);
  for (core::Size j = 1; j <= positions.size(); ++j) {
    sampler.consider_sample(positions[j], weights[j]);
  }
  sampler.samples(&samples);
  return samples[1];
}

void StarTreeBuilder::tear_down(core::pose::Pose* pose) {
  assert(pose);
  if (virtual_res_ == -1) {
    return;
  } else {
    pose->conformation().delete_residue_slow(virtual_res_);
    virtual_res_ = -1;
  }
}

}  // namespace nonlocal
}  // namespace protocols
