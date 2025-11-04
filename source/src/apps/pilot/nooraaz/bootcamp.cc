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
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <numeric/random/random.hh>
#include <protocols/monte_carlo/GenericMonteCarloMover.hh>
#include <protocols/moves/MonteCarlo.hh>

int main( int argc, char ** argv ) {
  std::cout << "Hello World!" << std::endl;
  devel::init( argc, argv );
  utility::vector1< std::string > filenames = basic::options::option[ basic::options::OptionKeys::in::file::s ].value();
  if ( filenames.size() > 0 ) {
    std::cout << "You entered: " << filenames[ 1 ] << " as the PDB file to be read" << std::endl;
    // 1. Initialize a pose, a score function object, and then a MonteCarloObject
    core::pose::PoseOP mypose = core::import_pose::pose_from_file( filenames[1] );
    double N = static_cast<double>( mypose->size() );
    std::cout << "Pose length is : " << N << std::endl;
    core::scoring::ScoreFunctionCOP sfxn = core::scoring::get_score_function();
    protocols::moves::MonteCarloOP mover = utility::pointer::make_shared<protocols::moves::MonteCarlo>( *mypose, *sfxn, 0.5 );
    core::Real score = sfxn->score( *mypose );
    std::cout << "The initial score is: " << score << std::endl;
    
    for (int j=0;j<5;j++){
        std::vector<double> scores;
        std::vector<bool>   accept;
				for (int i=0;i<100;i++){
							// 3. perturb the phi/psi values for your Pose
							core::Real random_uni = numeric::random::uniform();
							core::Size random_res_num = static_cast< core::Size > ( random_uni * N + 1 );
							core::Real random_pert1 = numeric::random::gaussian();
							core::Real random_pert2 = numeric::random::gaussian();
							core::Real orig_phi = mypose->phi( random_res_num );
							core::Real orig_psi = mypose->psi( random_res_num );
							mypose->set_phi( random_res_num, orig_phi + random_pert1 );
							mypose->set_psi( random_res_num, orig_psi + random_pert2 );

							// boltzmann minimizer
							bool if_accept = mover->boltzmann( *mypose );

							// print out the score
							core::Real score = sfxn->score( *mypose );
							scores.push_back(score);
							accept.push_back(if_accept);
				}
			  double average_score = std::accumulate(scores.begin(), scores.end(), 0.0) / scores.size() ;
			  double success_rate  = std::accumulate(accept.begin(), accept.end(), 0.0) / accept.size() ;
			  std::cout << "The average score after 100 steps is : " << average_score << std::endl;
			  std::cout << "The acceptane rate is : " << success_rate << std::endl;
        
    }
  } else {
    std::cout << "You didn’t provide a PDB file with the -in::file::s option" << std::endl;
    return 1;
  }
  return 0;
}
