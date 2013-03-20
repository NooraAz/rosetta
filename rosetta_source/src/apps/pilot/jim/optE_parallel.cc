// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

// libRosetta headers
//#include <basic/options/option.hh>

/// MPI
#ifdef USEMPI
#include <mpi.h>
#endif

/// Core headers
#include <devel/init.hh>


// Protocol headers
#include <protocols/optimize_weights/IterativeOptEDriver.hh>
#include <utility/excn/Exceptions.hh>
// #include <utility/file/file_sys_util.hh>
// #include <utility/exit.hh>


int
main( int argc, char * argv [] )
{
    try {
    #ifdef USEMPI
    	MPI_Init(&argc, &argv);
    #endif


    	devel::init( argc, argv );

    	protocols::optimize_weights::IterativeOptEDriver driver;
    	driver.go();

    #ifdef USEMPI
    	MPI_Finalize();
    #endif
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
    return 0;
}
