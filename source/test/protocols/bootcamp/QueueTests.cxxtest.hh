// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  p/r/o/t/o/c/o/l/s///b/o/o/t/c/a/m/p/QueueTests.cxxtest.hh
/// @brief  lab3
/// @author Noora Azadvari (eazadvar@uoregon.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>
#include <protocols/bootcamp/Queue.hh>

// Project Headers


// Core Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("QueueTests");


class QueueTests : public CxxTest::TestSuite {
	//Define Variables
	private:
    protocols::bootcamp::Queue Q_; 

public:

	void setUp() {
		core_init();
    protocols::bootcamp::Queue myq;
    bool isQ = myq.is_empty();
    TS_TRACE(isQ);
    core::Size q_size = Q_.size();
    TS_TRACE(" size of Q in the beginning is ");
    TS_TRACE( q_size );
    Q_.enqueue( std::string("first instance") );
    Q_.enqueue( std::string("second instance") );
    Q_.enqueue( std::string("third instance") );
    Q_.enqueue( std::string("fourth instance") );
    core::Size q_size2 = Q_.size();
    TS_TRACE(" size of after dding 4 elements is ");
    TS_TRACE( q_size2 );
   
	}

	void tearDown() {

	}



	void test_first() {
    TS_TRACE( "Running my first unit test!" );
    TS_ASSERT( true );
	}


};
