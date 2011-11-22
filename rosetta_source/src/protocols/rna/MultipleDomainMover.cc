// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_protocols
/// @brief protocols that are specific to MultipleDomainMover
/// @detailed
/// @author Rhiju Das

#include <protocols/rna/MultipleDomainMover.hh>

#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

#include <protocols/coarse_rna/CoarseRNA_LoopCloser.hh>
#include <protocols/moves/RigidBodyMover.hh>
#include <protocols/rna/AllowInsert.hh>

// Utility headers
#include <utility/vector1.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>

//C++ headers
#include <vector>
#include <string>
#include <sstream>

using namespace core;
using basic::T;

typedef  numeric::xyzMatrix< Real > Matrix;
static numeric::random::RandomGenerator RG(1277320);  // <- Magic number, do not change it!

namespace protocols {
namespace rna {

/////////////////////////////////////
	MultipleDomainMover::MultipleDomainMover( pose::Pose const & pose, protocols::coarse_rna::CoarseRNA_LoopCloserOP rna_loop_closer  ):
	verbose_( false ),
	rot_mag_( 10.0 ),
	trans_mag_( 0.5 ),
	num_domains_( 0 ),
	rna_loop_closer_( rna_loop_closer )
{
	Mover::type("MultipleDomainMover");
	initialize( pose, rna_loop_closer_->allow_insert() );
}


////////////////////////////////////////////////////////////////
void
MultipleDomainMover::apply( core::pose::Pose & pose )
{
	apply_and_return_jump( pose );
}

///////////////////////////////////////////////////////////////////////////////
	////////////////////////
std::string
MultipleDomainMover::get_name() const {
	return "MultipleDomainMover";
}


////////////////////////////////////////////////////////////////
Size
MultipleDomainMover::apply_and_return_jump( core::pose::Pose & pose )
{
  Size const n = static_cast<Size> ( RG.uniform() * num_domains_ ) + 1;
	return apply_at_domain( pose, n );
}


////////////////////////////////////////////////////////////////
Size
MultipleDomainMover::apply_at_domain( core::pose::Pose & pose, Size const & n )
{
	rb_movers_[ n ]->apply( pose );
	slide_back_to_origin( pose );
	return jump_numbers_[ n ];
}


////////////////////////////////////////////////////////////////
void
MultipleDomainMover::initialize( pose::Pose const & pose, AllowInsertOP allow_insert ){
	//	std::cout << "HELLO! " << pose.residue_type( pose.total_residue() ).name3() << std::endl;
	if ( pose.residue_type( pose.total_residue() ).name3() != "XXX" /*virtual residue*/ ) return;
	setup_jump_numbers_and_partner( pose );
	setup_ok_for_centroid_calculation( allow_insert );
	initialize_rigid_body_movers();
}

/////////////////////////////////////
void
MultipleDomainMover::randomize_pose_rigid_bodies( pose::Pose & pose ){
	randomize_orientations( pose );
	try_to_slide_into_contact( pose );
	close_all_loops( pose );
	slide_back_to_origin( pose );
}

/////////////////////////////////////
Vector
MultipleDomainMover::get_centroid( pose::Pose const & pose ){

	Vector cen( 0.0, 0.0, 0.0 );
	Size nres( 0 );
	// Look at all residues except anchor (last residue).
	for ( Size i = 1; i < pose.total_residue(); i++ ){
		if ( ok_for_centroid_calculation_[ i ] ){
			cen += pose.xyz( core::id::AtomID( 1, i ) );
			nres += 1;
		}
	}
	cen /= nres;
	return cen;
}

/////////////////////////////////////
void
MultipleDomainMover::slide_back_to_origin( pose::Pose & pose ){

	using namespace core::kinematics;

	if ( num_domains_ < 2 ) return;

	Vector cen = get_centroid( pose );
	//	std::cout << "CENTROID1: " << cen(1) << ' ' << cen(2) << ' ' << cen(3) << std::endl;

	for ( Size n = 1; n <= num_domains_; n++ ){
		Size jumpno( jump_numbers_[ n ] );
		Jump j( pose.jump( jumpno ) );
		Stub stub = pose.conformation().upstream_jump_stub( jumpno );
		Vector new_translation = j.get_translation() - stub.M.transposed() * cen;
		j.set_translation( new_translation);
		pose.set_jump( jumpno, j );
	}

	//	cen = get_centroid( pose );
	//	std::cout << "CENTROID2: " << cen(1) << ' ' << cen(2) << ' ' << cen(3) << std::endl;

}



/////////////////////////////////////
void MultipleDomainMover::setup_jump_numbers_and_partner( pose::Pose const & pose ){

	using namespace protocols::moves;

	// create rigid body mover for segment 1
	Size const virt_res = pose.total_residue();

	jump_numbers_.clear();
	partner_.clear();

	for ( Size n = 1; n <= pose.fold_tree().num_jump(); n++ ) {
		Size const i = pose.fold_tree().upstream_jump_residue( n );
		Size const j = pose.fold_tree().downstream_jump_residue( n );
		if ( i == virt_res || j == virt_res ){

			jump_numbers_.push_back( n );

			if ( i == virt_res ) {
				//partner.push_back( partner_downstream );
				partner_.push_back( partner_upstream );
			} else {
				partner_.push_back( partner_downstream );
				//partner.push_back( partner_upstream );
			}
		}
	}

	num_domains_ = jump_numbers_.size();
}

/////////////////////////////////////
void
MultipleDomainMover::setup_ok_for_centroid_calculation( AllowInsertOP & allow_insert ){
	// Need to find jump number [Can try alternate constructor later]
	ok_for_centroid_calculation_.clear();

	for ( Size i = 1; i <= allow_insert->nres(); i++) {
		ok_for_centroid_calculation_.push_back(  !allow_insert->get( core::id::AtomID( 2, i) ) );
		if (verbose_) std::cout <<  "OK " << i << " " << ok_for_centroid_calculation_[ i ] << std::endl;
	}

}

/////////////////////////////////////
void
MultipleDomainMover::randomize_orientations( pose::Pose & pose ) {

	using namespace protocols::moves;

	for ( Size n = 1; n<= num_domains_; n++ ){
		RigidBodyRandomizeMover rb( pose, jump_numbers_[ n ], partner_[ n ] );
		rb.apply( pose );
	}

	if (verbose_) pose.dump_pdb( "random.pdb" );

}


/////////////////////////////////////////////////////////////////////////////////////
// slide into contact - actually just try to get overlap at cutpoint_closed
////////////////////////////////////////////////////////////////////////////////////
void MultipleDomainMover::try_to_slide_into_contact( pose::Pose & pose ) {

	using namespace core::kinematics;

	utility::vector1< bool > slid_into_contact( pose.total_residue(), false );

	for ( Size n = 2; n <= num_domains_; n++ ){ // no need to move domain 1.

		Size jumpno( jump_numbers_[ n ] );

		// find cutpoint_closed
		rna_loop_closer_->apply_after_jump_change( pose, jumpno );
		utility::vector1< Size > const & cutpos_list = rna_loop_closer_->cutpos_list();

		Size cutpoint_closed( 0 );
		for (Size i = 1; i <= cutpos_list.size() ; i++ ){
			if ( !slid_into_contact[ cutpos_list[ i ] ]  ){
				cutpoint_closed = cutpos_list[ i ]; break;
			}
		}
		if (verbose_) std::cout << "CUTPOINT_CLOSED: " << cutpoint_closed << std::endl;
		if ( cutpoint_closed == 0 ) continue;

		// find locations of overlap atoms. 	// figure out translation
		Vector diff = pose.residue( cutpoint_closed ).xyz( "OVL1" ) - pose.residue( cutpoint_closed+1 ).xyz( " P  " );
		if ( verbose_) std::cout << "OLD VEC " <<  diff.length() << std::endl;
		// apply translation
		// This could be generalized. Note that this relies on one partner being
		// the VRT residue at the origin!!!
		Jump j( pose.jump( jumpno ) );
		if (verbose_){
			std::cout << "OLD JUMP  " << j << std::endl;
			std::cout << "UPSTREAM DOWNSTREAM " << pose.fold_tree().upstream_jump_residue( jumpno ) <<' ' <<  pose.fold_tree().downstream_jump_residue( jumpno ) << std::endl;
		}

		ObjexxFCL::FArray1D< bool > const & partition_definition = rna_loop_closer_->partition_definition();
		int sign_jump = ( partition_definition( cutpoint_closed ) == partition_definition( pose.fold_tree().downstream_jump_residue( jumpno ) )  ? -1: 1);

		Stub stub = pose.conformation().upstream_jump_stub( jumpno );
		if (verbose_) std::cout << ' ' << stub.M.col_x()[0]<< ' ' << stub.M.col_x()[1]<< ' ' << stub.M.col_x()[2] << std::endl;
		Vector new_translation = j.get_translation() + sign_jump * stub.M.transposed() * diff;
		j.set_translation( new_translation );//+ Vector( 50.0, 0.0, 0.0 ) );
		pose.set_jump_now( jumpno, j );
		if (verbose_) std::cout << "NEW JUMP  " << pose.jump( jumpno ) << std::endl;
		Vector diff2 = pose.residue( cutpoint_closed ).xyz( "OVL1" ) - pose.residue( cutpoint_closed+1 ).xyz( " P  " );
		if (verbose_) std::cout << "NEW VEC " <<  diff2.length() << std::endl << std::endl;

		slid_into_contact[ cutpoint_closed ] = true;
	}

	if (verbose_) pose.dump_pdb( "slide.pdb" );

}

//////////////////////////////////////////////////////////////////////////////////
void
MultipleDomainMover::close_all_loops( pose::Pose & pose ) {
	for ( Size n = 2; n <= num_domains_; n++ ) {
		rna_loop_closer_->apply_after_jump_change( pose, n );
	}
	if ( verbose_ ) pose.dump_pdb( "closed.pdb" );
}



/////////////////////////////////////
void
MultipleDomainMover::initialize_rigid_body_movers(){

	using namespace protocols::moves;

	rb_movers_.clear();

	std::cout << "NUM_DOMAINS: " << num_domains_ << std::endl;
	for ( Size n = 1; n <= num_domains_; n++ ){
		rb_movers_.push_back( new RigidBodyPerturbMover( jump_numbers_[ n ],
																										 rot_mag_, trans_mag_,
																										 partner_[ n ], ok_for_centroid_calculation_ ) );
	}

}

/////////////////////////////////////
void
MultipleDomainMover::update_rot_trans_mag( Real const & rot_mag, Real const & trans_mag ){
	rot_mag_ = rot_mag;
	trans_mag_ = trans_mag;
	for ( Size n = 1; n <= num_domains_; n++ ){
		rb_movers_[ n ]->rot_magnitude( rot_mag_ );
		rb_movers_[ n ]->trans_magnitude( rot_mag_ );
	}
}



} // namespace rna
} // namespace protocols
