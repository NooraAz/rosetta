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

  utility::vector1< std::pair< core::Size, core::Size > >
  identify_secondary_structure_spans( std::string const & ss_string ) {
    utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries;
    core::Size strand_start = -1;
    for ( core::Size ii = 0; ii < ss_string.size(); ++ii ) {
      if ( ss_string[ ii ] == 'E' || ss_string[ ii ] == 'H'  ) {
        if ( int( strand_start ) == -1 ) {
          strand_start = ii;
        } else if ( ss_string[ii] != ss_string[strand_start] ) {
          ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
          strand_start = ii;
        }
       } else {
         if ( int( strand_start ) != -1 ) {
           ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
           strand_start = -1;
         }
       }
    }

    if ( int( strand_start ) != -1 ) {
    // last residue was part of a ss-eleemnt
      ss_boundaries.push_back( std::make_pair( strand_start+1, ss_string.size() ));
    }
    for ( core::Size ii = 1; ii <= ss_boundaries.size(); ++ii ) {
       std::cout << "SS Element " << ii << " from residue "
       << ss_boundaries[ ii ].first << " to "
       << ss_boundaries[ ii ].second << std::endl;
      }
    return ss_boundaries;
    }
    
  void test_known_strings(){
    std::string s1 = "EEEEEHHHHHHHHEEEEEIGNOREEEEEEHHHHHHHHHHHEEEEEHHHH";
    utility::vector1< std::pair< core::Size, core::Size > > s1_pairs = {{1,5},{6,13},{14,18},{24,29},{30,40},{41,45},{46,49}};
    std::string s2 = "HHHHHHH HHHHHHHHHHHHHHHHHHHHHHHHEEEEEEEEEEHHHHHHHEEEEHHH";
    utility::vector1< std::pair< core::Size, core::Size > > s2_pairs = {{1,7},{9,32},{33,42},{43,49},{50,53},{54,56}};
    std::string s3 = "EEEEEEEEE EEEEEEEE EEEEEEEEE H EEEEE H H H EEEEEEEE";
    utility::vector1< std::pair< core::Size, core::Size > > s3_pairs = {{1,9},{11,18},{20,28},{30,30},{32,36},{38,38},{40,40},{42,42},{44,51}};

    TS_ASSERT_EQUALS(identify_secondary_structure_spans(s1), s1_pairs);
    TS_ASSERT_EQUALS(identify_secondary_structure_spans(s2), s2_pairs);
    TS_ASSERT_EQUALS(identify_secondary_structure_spans(s3), s3_pairs);
  }
  core::kinematics::FoldTree fold_tree_from_dssp_string(std::string const & dssp){
    // convert ss string to foldtree
    core::kinematics::FoldTree foldtree;
    utility::vector1< std::pair< core::Size, core::Size > > edges = identify_secondary_structure_spans(dssp);
    for (auto e = edges.begin(); e!= edges.end(); e++ ){
      foldtree.add_edge(e->first, e->second, core::kinematics::Edge::PEPTIDE );
    }
    int midppint_of_first_block = static_cast<int>( (edges[1].first + edges[1].second ) / 2 );
    int n = 1;
    for (auto e = edges.begin() + 1; e!= edges.end(); e++){
      int midpoint = static_cast<int>( (e->first + e->second ) / 2 );
      foldtree.add_edge(midppint_of_first_block, midpoint, n);
      n = n + 1;
    }
    return foldtree;
  }

  core::kinematics::FoldTree fold_tree_from_ss(core::pose::Pose const & pose){
    core::scoring::dssp::Dssp dssp = core::scoring::dssp::Dssp( pose );
    std::string secondary_structure = dssp.get_dssp_secstruct();
    return fold_tree_from_dssp_string( secondary_structure );
  }
 
    
};
