// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/FragmentAllAtomCrmsd.hh
/// @brief  Object that scores a fragment by its crmsd to the native
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_FragmentAllAtomCrmsd_hh
#define INCLUDED_protocols_frag_picker_scores_FragmentAllAtomCrmsd_hh

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

// mini

#include <core/types.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <core/pose/Pose.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoringMethod.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace ObjexxFCL;

typedef utility::vector1<utility::vector1<Real> > Matrix;

/// @brief  scores a fragment by its crmsd to the given reference structure
class FragmentAllAtomCrmsd: public FragmentScoringMethod {
public:

	/// @brief  creates a crmsd-based scoring function.
	/// @details fragments will be compared to a given pose, which should have the same number of residues a the query sequence
	FragmentAllAtomCrmsd(Size, Real, bool, std::string, core::pose::PoseOP);

	/// @brief  creates a crmsd-based scoring function.
	/// @details fragments will be compared to given coordinates, which should have the same number of residues a the query sequence
	FragmentAllAtomCrmsd(Size, Real, bool, std::string, utility::vector1< utility::vector1<Real> >);

	~FragmentAllAtomCrmsd();

	// Undefined, commenting out to fix PyRosetta build  void do_caching(VallChunkOP);
	// Undefined, commenting out to fix PyRosetta build  void clean_up();
	bool score(FragmentCandidateOP, FragmentScoreMapOP);
	// Undefined, commenting out to fix PyRosetta build  bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);

private:
	std::string cached_scores_id_;
	Size total_ref_resids_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// Size frag_resids_;
	std::string query_sequence_;
	pose::PoseOP reference_pose_;
	pose::PoseOP fragment_pose_;
	FArray2D_double reference_coordinates_;
	FArray2D_double chunk_coordinates_;
	FArray2D_double fragment_coordinates_;
	FArray2D_double frag_pos_coordinates_;
	FArray1D_double weights_;

	void fill_coords(core::pose::Pose const &, FArray2_double &, Size, std::string);
};

/// @brief  Maker class that produces a new FragmentAllAtomCrmsd object
class MakeFragmentAllAtomCrmsd: public MakeFragmentScoringMethod {
public:

	MakeFragmentAllAtomCrmsd() :
		MakeFragmentScoringMethod("FragmentAllAtomCrmsd") {
	}

	FragmentScoringMethodOP make(Size, Real, bool, FragmentPickerOP, std::string);
};

} // scores
} // frag_picker
} // protocols

#endif // INCLUDED_protocols_frag_picker_scores_FragmentAllAtomCrmsd_HH
