// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/rna/RNA_IdealCoord.hh
/// @brief Apply ideal RNA geometry to a residue or a pose
/// @author Fang-Chieh Chou

#ifndef INCLUDED_core_pose_rna_RNA_IdealCoord_HH
#define INCLUDED_core_pose_rna_RNA_IdealCoord_HH

// Unit headers
#include <core/pose/rna/RNA_IdealCoord.fwd.hh>

// Package headers
#include <core/pose/Pose.fwd.hh>

// Project headers
#include <core/conformation/Residue.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <string>

namespace core {
namespace pose {
namespace rna {

class RNA_IdealCoord : public utility::pointer::ReferenceCount {
public:

	RNA_IdealCoord();
	~RNA_IdealCoord();

	//Apply ideal coords to one residue. Keep the backbone torsion values by default
	void apply( Pose & pose, Size const seqpos, Size pucker, bool const keep_backbone_torsion = true ) const;

	//Apply ideal coords to whole pose.
	//pucker_conformations: 0 for keeping pucker, 1 for North, 2 for South
	void apply( Pose & pose, utility::vector1 < Size > const & puckers, bool const keep_backbone_torsion = true ) const;

	//Apply ideal coords to whole pose. Keep all pucker state
	void apply( Pose & pose, bool const keep_backbone_torsion = true ) const;

private:
	void init();
	bool is_torsion_exists(Pose const & pose, id::TorsionID const & torsion_id) const;
	utility::vector1 < PoseOP > ref_pose_list_;
	std::string const path_;
	Real delta_cutoff_;
};

}
}
}

#endif
