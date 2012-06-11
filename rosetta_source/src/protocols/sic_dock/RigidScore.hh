// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_protocols_sic_dock_RigidScore_hh
#define INCLUDED_protocols_sic_dock_RigidScore_hh

#include <protocols/sic_dock/RigidScore.fwd.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <core/id/AtomID_Map.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/sic_dock/xyzStripeHashPose.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <protocols/loophash/LoopHashLibrary.hh>

namespace protocols {
namespace sic_dock {

	struct Vec3 { numeric::xyzVector<core::Real> a,b,c; };
	typedef utility::vector1<Vec3> Vec3s;

class RigidScore : public utility::pointer::ReferenceCount {
protected:
	typedef core::Real Real;
	typedef core::Size Size;
	typedef core::pose::Pose Pose;
	typedef core::kinematics::Stub Stub;
	typedef numeric::xyzVector<Real> Vec;
	typedef numeric::xyzMatrix<Real> Mat;
	typedef utility::vector1<Vec> Vecs;
	typedef utility::vector1<Real> Reals;
	typedef utility::vector1<Size> Sizes;
	typedef utility::vector1<Stub> Stubs;
	typedef utility::vector1<RigidScoreCOP> Scores;
public:

	virtual ~RigidScore() {}

	virtual
	core::Real
	score(
		Stub const & x1,
		Stub const & x2
	) const = 0;

};


class CBScore : public RigidScore {
public:
	CBScore(Real clash_dis, Real contact_dis);
	virtual ~CBScore(){}
	CBScore(
		Pose const & pose1,
		Pose const & pose2,
		Real clash_dis,
		Real contact_dis
	);
	core::Real score( Stub const & x1, Stub const & x2 ) const;
//private:
	bool const hash_pose1_;
	core::Real const clash_dis_, contact_dis_;
	Reals const weights_;
	Vecs const points_;
	xyzStripeHashPose const xyzhash_;
	// Pose const & pose1_,pose2_;
};


class LinkerScore : public RigidScore {
public:
	LinkerScore(
		Pose const & pose1,
		Pose const & pose2,
		Size max_loop_len
	);
	virtual ~LinkerScore(){}
	core::Real  score( Stub const & x1, Stub const & x2 ) const;
	void dump_linkers( Stub const & x1, Stub const & x2 ) const;
private:
	protocols::loophash::LoopHashLibraryOP loop_hash_library_;
	Sizes const loopsizes_;
	Pose const & pose1_,pose2_;
	Vec3s lowers1_,uppers1_,lowers2_,uppers2_;
	Real max_dis2_;
};



class EdgeStandScore : public RigidScore {
public:
	EdgeStandScore();
	virtual ~EdgeStandScore(){}
	core::Real score( Stub const & x1, Stub const & x2 ) const;
private:
	Vecs donors,acceptors;
};

class HelixScore : public RigidScore {
public:
	HelixScore();
	virtual ~HelixScore(){}
	core::Real score( Stub const & x1, Stub const & x2 ) const;
private:
};

class BuriedPolarScore : public RigidScore {
public:
	BuriedPolarScore();
	virtual ~BuriedPolarScore(){}
	core::Real score( Stub const & x1, Stub const & x2 ) const;
private:
	Vecs polars;
};

////// composite scores

class JointScore : public RigidScore {
public:
	JointScore(){}
	JointScore(
		Scores scores,
		Reals weights
	);
	void add_score(RigidScoreCOP score, Real weight);
	virtual ~JointScore(){}
	core::Real score( Stub const & x1, Stub const & x2 ) const;
private:
	Scores scores_;
	Reals weights_;

};



// class CachedScore : public RigidScore {
// public:
// 	CachedScore(RigidScoreCOP score);
// 	virtual ~CachedScore(){}
// 	core::Real score( Stub const & x1, Stub const & x2 ) const;
// private:
// 	RigidScoreCOP score_;
// 	// some kind of 6 dof hash
// };






} // namespace sic_dock
} // namespace protocols

#endif
