// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// awatkins: based heavily on kdrew/oop_creator.cc

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/ncbb/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>


#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>

#include <core/kinematics/MoveMap.hh>


// Mover headers
#include <protocols/ncbb/util.hh>
#include <protocols/minimization_packing/MinMover.hh>
//#include <protocols/ncbb/hbs/HbsRandomSmallMover.hh>
#include <protocols/ncbb/hbs/HbsPatcher.hh>
#include <protocols/simple_moves/chiral/ChiralMover.hh>


//Basic headers

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
//#include <basic/options/keys/OptionKeys.hh>
//#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <sstream>

#include <basic/options/keys/OptionKeys.hh> // AUTO IWYU For OptionKeys
#include <utility/file/FileName.fwd.hh> // AUTO IWYU For FileName

using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace pose;
using namespace protocols;
using namespace protocols::ncbb;
using namespace protocols::moves;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::hbs;
using namespace protocols::simple_moves::chiral;
using namespace core::pack::task;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::id;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

static basic::Tracer TR("PeptoidHBS");

int
main( int argc, char* argv[] )
{
	try {

		devel::init(argc, argv);

		// create score function
		scoring::ScoreFunctionOP score_fxn( scoring::get_score_function() );
		scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn(*score_fxn);

		if ( score_fxn->has_zero_weight( atom_pair_constraint ) ) {
			score_fxn->set_weight( atom_pair_constraint, 0.1 );
		}

		if ( score_fxn->has_zero_weight( dihedral_constraint ) ) {
			score_fxn->set_weight( dihedral_constraint, 0.1 );
		}

		if ( score_fxn->has_zero_weight( angle_constraint ) ) {
			score_fxn->set_weight( angle_constraint, 0.1 );
		}

		chemical::ResidueTypeSetCOP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

		Pose pose;
		make_pose_from_sequence( pose, "X[201:NtermPeptoidFull]A[DALA]GAAAAAA[ALA:MethylatedCtermProteinFull]", *restype_set, false );
		//make_pose_from_sequence( pose, "A[ALA:NtermProteinFull]AGAAAAAA[ALA:MethylatedCtermProteinFull]", *restype_set, false );
		scoring::constraints::add_fa_constraints_from_cmdline_to_pose(pose);

		using namespace core::scoring;
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;

		for ( Size i = 1; i <= pose.size(); ++i ) {

			if ( i  < pose.size() - 3 ) {
				pose.add_constraint( utility::pointer::make_shared< AtomPairConstraint >( AtomID( pose.residue( i ).atom_index( "O" ), i ),
					AtomID( pose.residue( i+4 ).atom_index( "H" ), i+4 ),
					utility::pointer::make_shared< HarmonicFunc >( 1.8, 0.1 ) ) );
			}
			// constrain non pre-peptoid omega
			if ( i < pose.size() && !pose.residue(i+1).type().is_peptoid() ) {
				pose.add_constraint( utility::pointer::make_shared< DihedralConstraint >( AtomID( pose.residue( i ).atom_index( "CA" ), i ),
					AtomID( pose.residue( i ).atom_index( "C"  ), i ),
					AtomID( pose.residue( i+1 ).atom_index( "N"  ), i+1 ),
					AtomID( pose.residue( i+1 ).atom_index( "CA"  ), i+1 ),
					utility::pointer::make_shared< CircularHarmonicFunc >( 3.14159, 0.04 ) ) );
			}


			// super loose dihedral constraints non pre-peptoid omega
			if ( i > 1 && i < pose.size() ) {
				pose.add_constraint( utility::pointer::make_shared< DihedralConstraint >( AtomID( pose.residue( i ).atom_index( "N" ), i ),
					AtomID( pose.residue( i ).atom_index( "CA"  ), i ),
					AtomID( pose.residue( i ).atom_index( "C"  ), i ),
					AtomID( pose.residue( i+1 ).atom_index( "N"  ), i+1 ),
					utility::pointer::make_shared< CircularHarmonicFunc >( -45.0*3.14159/180, 0.2 ) ) );
			}

			if ( i > 1 ) {
				pose.add_constraint( utility::pointer::make_shared< DihedralConstraint >( AtomID( pose.residue( i-1 ).atom_index( "C" ), i-1 ),
					AtomID( pose.residue( i ).atom_index( "N"  ), i ),
					AtomID( pose.residue( i ).atom_index( "CA"  ), i ),
					AtomID( pose.residue( i ).atom_index( "C"  ), i ),
					utility::pointer::make_shared< CircularHarmonicFunc >( -65.0*3.14159/180, 0.2 ) ) );
			}

			pose.set_phi( i, -65 );
			pose.set_psi( i, -45 );
			pose.set_omega( i, 180 );
		}
		pose.add_constraint( utility::pointer::make_shared< AtomPairConstraint >( AtomID( pose.residue( pose.size() - 3 ).atom_index( "O" ), pose.size() - 3 ),
			AtomID( pose.residue( pose.size() ).atom_index( "HM" ), pose.size() ),
			utility::pointer::make_shared< HarmonicFunc >( 1.8, 0.1 ) ) );

		pose.add_constraint( utility::pointer::make_shared< DihedralConstraint >( AtomID( pose.residue( pose.size() ).atom_index( "N" ), pose.size() ),
			AtomID( pose.residue( pose.size() ).atom_index( "CA"  ), pose.size() ),
			AtomID( pose.residue( pose.size() ).atom_index( "C"  ), pose.size() ),
			AtomID( pose.residue( pose.size() ).atom_index( "NM"  ), pose.size() ),
			utility::pointer::make_shared< CircularHarmonicFunc >( -45.0*3.14159/180, 0.2 ) ) );

		pose.dump_pdb( "first.pdb");

		// create move map for minimization
		kinematics::MoveMapOP mm( new kinematics::MoveMap() );
		mm->set_chi( true );
		for ( Size i = 1; i <= pose.size(); ++i ) {
			mm->set_bb( i, true);
		}

		TR << *mm << std::endl;
		//pose.conformation().declare_chemical_bond( 1, "CYH", 3, "CZH" );

		// create minimization mover
		minimization_packing::MinMoverOP minM( new protocols::minimization_packing::MinMover( mm, score_fxn, "lbfgs_armijo_nonmonotone", 0.01, true ) );
		minM->cartesian( true );
		for ( Real wt = 0.1; wt <= 1.0; wt += 0.1 ) {
			TR << "Weight " << wt << std::endl;
			score_fxn->set_weight( atom_pair_constraint, wt );
			score_fxn->set_weight( dihedral_constraint, wt/4 );
			score_fxn->set_weight( angle_constraint, wt );
			//score_fxn->set_weight( atom_pair_constraint, 0 );
			//score_fxn->set_weight( dihedral_constraint, 0 );
			score_fxn->set_weight( angle_constraint, 0 );

			minM->apply( pose );
		}

		pose.dump_pdb( "second.pdb");

		HbsPatcher hp( 1 );
		hp.apply( pose );
		pose::ncbb::initialize_ncbbs(pose);

		pose.dump_pdb( "third.pdb");
		pose.conformation().declare_chemical_bond( 1, "CYH", 3, "CZH" );

		// create minimization mover
		for ( Real wt = 0.1; wt <= 1; wt += 0.1 ) {
			score_fxn->set_weight( atom_pair_constraint, wt );
			score_fxn->set_weight( dihedral_constraint, wt/4 );
			score_fxn->set_weight( angle_constraint, wt );
			//score_fxn->set_weight( atom_pair_constraint, 0 );
			//score_fxn->set_weight( dihedral_constraint, 0 );
			score_fxn->set_weight( angle_constraint, 0 );

			minM->apply( pose );
		}
		pose.dump_pdb( "fourth.pdb");

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;

}//main
