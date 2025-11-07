// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
//
// // vi: set ts=2 noet:
//
// //
//
// // (c) Copyright Rosetta Commons Member Institutions.
//
// // (c) This file is part of the Rosetta software suite and is made available under license.
//
// // (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
//
// // (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <iostream>
#include <vector>
#include <numeric> 
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <devel/init.hh>
#include <utility/pointer/owning_ptr.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/bootcamp/BootCampMover.hh>

 

int main( int argc, char ** argv ) {
  devel::init( argc, argv );
  utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value();
  if ( filenames.size() > 0 ) {
    std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
    core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );
		protocols::bootcamp::BootCampMover mover;
		mover.apply( *mypose );
  } else {
    std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
    return 1;
  }
  return 0;
}
