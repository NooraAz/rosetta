// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/match/ProteinSCSampler.cxxtest.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Utility headers

/// Project headers
#include <core/types.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/kinematics/FoldTree.hh>
// C++ headers

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.hh>




// --------------- Test Class --------------- //

class FoldTreeFromSSTests : public CxxTest::TestSuite {

public:


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.


	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.


	// --------------- Test Cases --------------- //
	void test_hello_world(){
		TS_ASSERT ( true );
  }

  void test_known_strings(){
    TS_ASSERT_EQUALS(identify_secondary_structure_spans(s1_), s1_pairs);
    TS_ASSERT_EQUALS(identify_secondary_structure_spans(s2_), s2_pairs);
    TS_ASSERT_EQUALS(identify_secondary_structure_spans(s3_), s3_pairs);
  }

  // secondary structure string -> fold tree 
  core::kinematics::FoldTree
	fold_tree_from_dssp_string ( std::string const & dssp ) {

    // Adding peptide edges
    //// From mid to both ways, for non-loops
    for (auto e = edges.begin(); e!= edges.end(); e++ ){
      int midpoint = static_cast<int>( (e->first + e->second ) / 2 );
      foldtree.add_edge(midpoint, e->first, core::kinematics::Edge::PEPTIDE );
      foldtree.add_edge(midpoint, e->second, core::kinematics::Edge::PEPTIDE );
    }
    //// From mid to both ways, for loops
		for (auto g = gaps.begin(); g!= gaps.end(); g++){
			int midloop = static_cast<int>( (g->first + g->second ) / 2 );
      foldtree.add_edge(midloop, g->first, core::kinematics::Edge::PEPTIDE );
      foldtree.add_edge(midloop, g->second, core::kinematics::Edge::PEPTIDE );
    }


    return foldtree;
  } // end of fold_tree_from_dssp_string

  private:
		std::string s1_ = "EEEEEHHHHHHHHEEEEEIGNOREEEEEEHHHHHHHHHHHEEEEEHHHH";
		std::string s2_ = "HHHHHHH HHHHHHHHHHHHHHHHHHHHHHHHEEEEEEEEEEHHHHHHHEEEEHHH";
		std::string s3_ = "EEEEEEEEE EEEEEEEE EEEEEEEEE H EEEEE H H H EEEEEEEE";
		std::string s4_ = "   EEEEEEE    EEEEEEE         EEEEEEEEE    EEEEEEEEEE   HHHHHH         EEEEEEEEE         EEEEE     ";
    utility::vector1< std::pair< core::Size, core::Size > > s1_pairs_ = 
      {{1,5},{6,13},{14,18},{24,29},{30,40},{41,45},{46,49}};
    utility::vector1< std::pair< core::Size, core::Size > > s2_pairs_ = 
      {{1,7},{9,32},{33,42},{43,49},{50,53},{54,56}};
    utility::vector1< std::pair< core::Size, core::Size > > s3_pairs_ = 
			{{1,9},{11,18},{20,28},{30,30},{32,36},{38,38},{40,40},{42,42},{44,51}};
};
