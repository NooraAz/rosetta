// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/pack/dunbrack/DunbrackEnergyTests_oligourea.cxxtest.hh
/// @brief  Test the Dunbrack energy with a residue type in which only some of the mainchain dihedrals affect rotamers (oligoureas).
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/dunbrack/DunbrackEnergy.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>

// Protocols Headers
#include <protocols/cyclic_peptide/PeptideStubMover.hh> //To make it easier to build poses

// Utility, etc Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("DunbrackEnergyTests_oligourea");


class DunbrackEnergyTests_oligourea : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-out:levels core.pack.rotamer_set.ContinuousRotamerSet:500 -extra_res_fa core/pack/dunbrack/dummy_type.params" );
	}

	void tearDown(){

	}

	void test_beta_aa_dunbrack_scoring(){
		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "GLY:NtermProteinFull", 1, true, "", 0, 1, "");
		builder.add_residue("Append", "B3C", 2, false, "N", 0, 1, "C");
		builder.add_residue("Append", "GLY:CtermProteinFull", 3, false, "N", 0, 2, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		pose.set_phi(2, 45);
		pose.set_theta(2, 33);
		pose.set_psi(2, -63);
		pose.set_omega(2, 180);
		pose.set_chi(1, 2, 33.7);
		pose.set_chi(2, 2, 62.7);

		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::fa_dun, 1.0 );

		core::Real const score1( (*sfxn)(pose) );
		core::Real score2;
		{
			core::pose::Pose pose2(pose);
			pose2.set_phi(2, 65);
			score2 = (*sfxn)(pose2);
		}
		core::Real score3;
		{
			core::pose::Pose pose3(pose);
			pose3.set_theta(2, -12);
			score3 = (*sfxn)(pose3);
		}
		core::Real score4;
		{
			core::pose::Pose pose4(pose);
			pose4.set_psi(2, 89);
			score4 = (*sfxn)(pose4);
		}
		core::Real score5;
		{
			core::pose::Pose pose5(pose);
			pose5.set_omega(2, -27);
			score5 = (*sfxn)(pose5);
		}

		TR.precision( 10 );
		TR << "Scores are: " << score1 << " " << score2 << " " << score3 << " " << score4 << " " << score5 << std::endl;

		TS_ASSERT_LESS_THAN( 0.1, std::abs(score2-score1) );
		TS_ASSERT_LESS_THAN( 0.1, std::abs(score3-score1) );
		TS_ASSERT_LESS_THAN( 0.1, std::abs(score4-score1) );
		TS_ASSERT_DELTA( score5, score1, 0.00001 );
	}

	void test_oligourea_dunbrack_scoring(){
		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "GLY:NtermProteinFull", 1, true, "", 0, 1, "");
		builder.add_residue("Append", "OU3_VAL", 2, false, "N", 0, 1, "C");
		builder.add_residue("Append", "GLY:CtermProteinFull", 3, false, "N", 0, 2, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		pose.set_phi(2, 45);
		pose.set_theta(2, 33);
		pose.set_psi(2, -63);
		pose.set_mu(2, 180);
		pose.set_omega(2, 180);
		pose.set_chi(1, 2, 33.7);

		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::fa_dun, 1.0 );

		core::Real const score1( (*sfxn)(pose) );
		core::Real score2;
		{
			core::pose::Pose pose2(pose);
			pose2.set_phi(2, 65);
			score2 = (*sfxn)(pose2);
		}
		core::Real score3;
		{
			core::pose::Pose pose3(pose);
			pose3.set_theta(2, -12);
			score3 = (*sfxn)(pose3);
		}
		core::Real score4;
		{
			core::pose::Pose pose4(pose);
			pose4.set_psi(2, 89);
			score4 = (*sfxn)(pose4);
		}
		core::Real score5;
		{
			core::pose::Pose pose5(pose);
			pose5.set_mu(2, -77);
			score5 = (*sfxn)(pose5);
		}
		core::Real score6;
		{
			core::pose::Pose pose6(pose);
			pose6.set_omega(2, -27);
			score6 = (*sfxn)(pose6);
		}

		TR.precision( 10 );
		TR << "Scores are: " << score1 << " " << score2 << " " << score3 << " " << score4 << " " << score5 << " " << score6 << std::endl;

		TS_ASSERT_LESS_THAN( 0.1, std::abs(score2-score1) );
		TS_ASSERT_LESS_THAN( 0.1, std::abs(score3-score1) );
		TS_ASSERT_LESS_THAN( 0.1, std::abs(score4-score1) );
		TS_ASSERT_DELTA( score5, score1, 0.00001 );
		TS_ASSERT_DELTA( score6, score1, 0.00001 );
	}

	void test_oligourea_dunbrack_scoring2(){
		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "GLY:NtermProteinFull", 1, true, "", 0, 1, "");
		builder.add_residue("Append", "OU3_DUMMYTYPE", 2, false, "N", 0, 1, "C");
		builder.add_residue("Append", "GLY:CtermProteinFull", 3, false, "N", 0, 2, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		pose.set_phi(2, 45);
		pose.set_theta(2, 33);
		pose.set_psi(2, -63);
		pose.set_mu(2, 180);
		pose.set_omega(2, 180);
		pose.set_chi(1, 2, 33.7);

		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::fa_dun, 1.0 );

		core::Real const score1( (*sfxn)(pose) );
		core::Real score2;
		{
			core::pose::Pose pose2(pose);
			pose2.set_phi(2, 65);
			score2 = (*sfxn)(pose2);
		}
		core::Real score3;
		{
			core::pose::Pose pose3(pose);
			pose3.set_theta(2, -12);
			score3 = (*sfxn)(pose3);
		}
		core::Real score4;
		{
			core::pose::Pose pose4(pose);
			pose4.set_psi(2, 89);
			score4 = (*sfxn)(pose4);
		}
		core::Real score5;
		{
			core::pose::Pose pose5(pose);
			pose5.set_mu(2, -77);
			score5 = (*sfxn)(pose5);
		}
		core::Real score6;
		{
			core::pose::Pose pose6(pose);
			pose6.set_omega(2, -27);
			score6 = (*sfxn)(pose6);
		}

		TR.precision( 10 );
		TR << "Scores are: " << score1 << " " << score2 << " " << score3 << " " << score4 << " " << score5 << " " << score6 << std::endl;

		TS_ASSERT_LESS_THAN( 0.1, std::abs(score2-score1) );
		TS_ASSERT_LESS_THAN( 0.1, std::abs(score3-score1) );
		TS_ASSERT_DELTA( score4, score1, 0.00001 );
		TS_ASSERT_LESS_THAN( 0.1, std::abs(score5-score1) );
		TS_ASSERT_DELTA( score6, score1, 0.00001 );
	}

	void test_oligourea_dunbrack_derivs(){
		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "GLY:NtermProteinFull", 1, true, "", 0, 1, "");
		builder.add_residue("Append", "OU3_VAL", 2, false, "N", 0, 1, "C");
		builder.add_residue("Append", "GLY:CtermProteinFull", 3, false, "N", 0, 2, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		pose.set_phi(2, 45);
		pose.set_theta(2, 33);
		pose.set_psi(2, -63);
		pose.set_mu(2, -175);
		pose.set_omega(2, 175);
		pose.set_chi(1, 2, 33.7);

		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::fa_dun, 1.0 );

		(*sfxn)(pose); //Ensure that scoring objects are set up.

		core::pack::dunbrack::DunbrackEnergy dun;

		core::Real const score1a (dun.eval_dof_derivative( core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("N"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 1 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
		core::Real const score1b (dun.eval_dof_derivative( core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("CA"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 2 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
		core::Real const score1c (dun.eval_dof_derivative( core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("CM"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 3 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
		core::Real const score1d (dun.eval_dof_derivative( core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("NU"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 4 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
		core::Real const score1e (dun.eval_dof_derivative( core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("C"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 5 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );

		TR.precision( 10 );
		TR << "Derivatives are: " << score1a << " " << score1b << " " << score1c << " " << score1d << " " << score1e << std::endl;

		TS_ASSERT_LESS_THAN( 0.1, std::abs(score1a) );
		TS_ASSERT_LESS_THAN( 0.1, std::abs(score1b) );
		TS_ASSERT_LESS_THAN( 0.1, std::abs(score1c) );
		TS_ASSERT_DELTA( score1d, 0.0, 0.00001 );
		TS_ASSERT_DELTA( score1e, 0.0, 0.00001 );
	}

	void test_oligourea_dunbrack_derivs2(){
		protocols::cyclic_peptide::PeptideStubMover builder;
		builder.add_residue("Append", "GLY:NtermProteinFull", 1, true, "", 0, 1, "");
		builder.add_residue("Append", "OU3_DUMMYTYPE", 2, false, "N", 0, 1, "C");
		builder.add_residue("Append", "GLY:CtermProteinFull", 3, false, "N", 0, 2, "C");
		core::pose::Pose pose;
		builder.apply(pose); //Build the peptide.

		pose.set_phi(2, 45);
		pose.set_theta(2, 33);
		pose.set_psi(2, -63);
		pose.set_mu(2, -175);
		pose.set_omega(2, 175);
		pose.set_chi(1, 2, 33.7);

		core::scoring::ScoreFunctionOP sfxn( new core::scoring::ScoreFunction );
		sfxn->set_weight( core::scoring::fa_dun, 1.0 );

		(*sfxn)(pose); //Ensure that scoring objects are set up.

		core::pack::dunbrack::DunbrackEnergy dun;

		core::Real const score1a (dun.eval_dof_derivative( core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("N"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 1 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
		core::Real const score1b (dun.eval_dof_derivative( core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("CA"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 2 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
		core::Real const score1c (dun.eval_dof_derivative( core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("CM"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 3 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
		core::Real const score1d (dun.eval_dof_derivative( core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("NU"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 4 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );
		core::Real const score1e (dun.eval_dof_derivative( core::id::DOF_ID( core::id::AtomID( pose.residue_type(2).atom_index("C"), 2 ), core::id::PHI ), core::id::TorsionID( 2, core::id::BB, 5 ), pose, *sfxn, pose.energies().onebody_energies(2) ) );

		TR.precision( 10 );
		TR << "Derivatives are: " << score1a << " " << score1b << " " << score1c << " " << score1d << " " << score1e << std::endl;

		TS_ASSERT_LESS_THAN( 0.1, std::abs(score1a) );
		TS_ASSERT_LESS_THAN( 0.1, std::abs(score1b) );
		TS_ASSERT_DELTA( score1c, 0.0, 0.00001 );
		TS_ASSERT_LESS_THAN( 0.1, std::abs(score1d) );
		TS_ASSERT_DELTA( score1e, 0.0, 0.00001 );
	}

};
