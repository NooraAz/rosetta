// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_Util.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_SWA_Util_HH
#define INCLUDED_protocols_swa_SWA_Util_HH


#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <string>
#include <map>
// AUTO-REMOVED #include <core/chemical/AA.hh>
#include <core/io/silent/RNA_SilentStruct.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
// AUTO-REMOVED #include <numeric/angle.functions.hh> // Need this to prevent the compiling error: 'principal_angle_degrees' is not a member of 'numeric' Oct 14, 2009
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
#include <core/id/AtomID_Map.fwd.hh>
// AUTO-REMOVED #include <set>

//Auto Headers
#include <core/id/AtomID.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>



using namespace core;
typedef  numeric::xyzMatrix< Real > Matrix;

namespace protocols {
namespace swa {

	typedef std::map< std::string, pose::PoseOP > PoseList;

	Size
	make_cut_at_moving_suite( pose::Pose & pose, Size const & moving_suite );

	bool
	Is_close_chain_break(pose::Pose const & pose);

	bool
	Contain_seq_num(Size const & seq_num, utility::vector1< Size > const & residue_list);

	void
	Output_boolean(std::string const & tag, bool boolean);

	void
	Output_boolean(bool boolean);

	void
	Output_movemap(kinematics::MoveMap const & mm, Size const total_residue);

	void
	Figure_out_moving_residues( kinematics::MoveMap & mm, pose::Pose const & pose,
															utility::vector1< Size > const & fixed_res,
															bool const move_takeoff_torsions = true	);

	void
	pdbslice( core::pose::Pose & new_pose,
						core::pose::Pose const & pose,
						utility::vector1< core::Size > const & slice_res );

	void
	pdbslice( pose::Pose & pose,
						utility::vector1< Size > const & slice_res );

	id::AtomID_Map< id::AtomID >
 	create_alignment_id_map(	pose::Pose & mod_pose, pose::Pose const & ref_pose,
														utility::vector1< Size > const & superimpose_res );


	utility::vector1< Size > const
	convert_to_working_res( utility::vector1< Size > const & res_vector,
													utility::vector1< Size > const & working_res ) ;


	scoring::constraints::ConstraintSetOP
	constraint_set_slice( scoring::constraints::ConstraintSetOP & cst_set,
												utility::vector1< Size > const & slice_res,
												pose::Pose const & pose,
												pose::Pose const & full_pose );

	utility::vector1< std::string > load_s_and_l();


	std::string
	get_file_name( std::string const & silent_file, std::string const & tag );

	void
	check_scorefxn_has_constraint_terms_if_pose_has_constraints( pose::Pose const & pose,
																															 scoring::ScoreFunctionOP & scorefxn );

	utility::vector1< Size >
	merge_vectors( utility::vector1< Size > const & vec1,
								 utility::vector1< Size > const & vec2 );

	void
	get_euler_angles( Real & alpha, Real & beta, Real & gamma, Matrix M1, Matrix M2, bool const verbose = true );

	void
	create_euler_rotation(
												Matrix & M,
												Real const & alpha,
												Real const & beta,
												Real const & gamma,
												Vector const & /* axis1 not actually used*/,
												Vector const & axis2,
												Vector const & axis3
												);

	void
	create_euler_rotation(
												Matrix & M,
												Real const & alpha,
												Real const & beta,
												Real const & gamma );

	void
	translate( pose::Pose & pose, Vector const shift,
						 pose::Pose const & ref_pose,
						 utility::vector1< Size > const & moving_res );

	void
	rotate( pose::Pose & pose, Matrix const M,
					pose::Pose const & ref_pose,
					utility::vector1< Size > const & moving_res,
					Vector const & centroid );

	void
	rotate( pose::Pose & pose, Matrix const M,
					pose::Pose const & ref_pose,
					utility::vector1< Size > const & moving_res );

	void
	get_base_centroid_and_rotation_matrix( pose::Pose const & pose, Size const i, Vector & centroid, Matrix & M );

	void
	translate_and_rotate_residue_to_origin( pose::Pose & pose, Size const i, utility::vector1< Size > const moving_res, bool const do_not_rotate = false );

	void
	translate_and_rotate_residue_to_origin( pose::Pose & pose, Size const i );

}
}

#endif
