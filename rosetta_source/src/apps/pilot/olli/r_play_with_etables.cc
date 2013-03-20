// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/olli/r_broker.cc
/// @author Oliver Lange

// keep these headers first for compilation with Visual Studio C++


// Etable Headers
#include <core/scoring/etable/EtableOptions.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>

#include <basic/database/open.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>


//Rosetta Scripts Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/viewer/viewers.hh>
#include <core/types.hh>

#include <devel/init.hh>
#include <basic/options/option.hh>

// Utility Headers

// Unit Headers
#include <protocols/moves/Mover.hh>

#include <basic/options/keys/parser.OptionKeys.gen.hh>
#include <basic/options/keys/jd2.OptionKeys.gen.hh>
// AUTO-REMOVED #include <utility/excn/Exceptions.hh>

#include <utility/vector1.hh>
#include <utility/excn/EXCN_Base.hh>
#include <basic/Tracer.hh>


// Utility headers
#include <devel/init.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <utility/exit.hh>

#include <basic/Tracer.hh>

// C++ headers
#include <iostream>

#include <utility/vector1.hh>


static basic::Tracer tr("main");


using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;


void*
my_main( void *)
{
  protocols::moves::MoverOP mover;//note that this is not instantiated and will crash if the job distributor actually tries to use it. That means that this can only be used with parser=true
  protocols::jd2::JobDistributor::get_instance()->go(mover);
  return 0 ;
}


int main( int argc, char * argv [] ) {
	try{
  devel::init( argc, argv ); //d

  chemical::AtomTypeSetCAP cptr_atom_type_set( chemical::ChemicalManager::get_instance()->atom_type_set( chemical::FA_STANDARD ) );
  chemical::AtomTypeSetAP local_atom_type_set( const_cast< chemical::AtomTypeSet* >( cptr_atom_type_set() ));
  local_atom_type_set->add_parameters_from_file(basic::database::full_name("chemical/atom_type_sets/fa_standard/extras/extra_soft_rep_params.txt"));
  scoring::etable::EtableOptions options;


  using namespace scoring;
  using namespace etable;
  chemical::AtomTypeSetCAP ptr_atom_type_set( local_atom_type_set.get() );

  options.lj_switch_dis2sigma = 0.6*0.9+0.1*0.91;
  scoring::etable::EtableOP etable0( new Etable( ptr_atom_type_set, options, std::string("SOFT10") ) );
  scoring::ScoringManager::get_instance()->add_etable( "FA_STANDARD_SOFT10", etable0 );

  options.lj_switch_dis2sigma = 0.6*0.8+0.2*0.91;
  scoring::etable::EtableOP etable1( new Etable( ptr_atom_type_set, options, std::string("SOFT20") ) );
  scoring::ScoringManager::get_instance()->add_etable( "FA_STANDARD_SOFT20", etable1 );

  options.lj_switch_dis2sigma = 0.6*0.5+0.5*0.91;
  scoring::etable::EtableOP etable2( new Etable( ptr_atom_type_set, options, "SOFT50" ) );
  scoring::ScoringManager::get_instance()->add_etable( "FA_STANDARD_SOFT50", etable2 );

  options.lj_switch_dis2sigma = 0.6*0.25+0.75*0.91;
  scoring::etable::EtableOP etable3( new Etable( ptr_atom_type_set, options, "SOFT75" ) );
  scoring::ScoringManager::get_instance()->add_etable( "FA_STANDARD_SOFT75", etable3 );

  scoring::ScoreFunctionOP scorefxn( scoring::getScoreFunction() );

  // pose::Pose pose;
  //  import_pose::pose_from_pdb( pose, option[ in::file::s ]()[1] );

  //  tr.Info << "full score: " << (*scorefxn)( pose ) << std::endl;

  scorefxn->set_etable( "FA_STANDARD_SOFT10" );
  //  tr.Info << "zak1 score: " << (*scorefxn)( pose ) << std::endl;

  scorefxn->set_etable( "FA_STANDARD_SOFT20" );
  //  tr.Info << "zak1 score: " << (*scorefxn)( pose ) << std::endl;

  scorefxn->set_etable( "FA_STANDARD_SOFT50" );
  //  tr.Info << "zak2 score: " << (*scorefxn)( pose ) << std::endl;

  scorefxn->set_etable( "FA_STANDARD_SOFT75" );
  //  tr.Info << "zak3 score: " << (*scorefxn)( pose ) << std::endl;


  //copied from rosetta_scripts
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  bool const view( option[ parser::view ] );
  protocols::moves::MoverOP mover;//note that this is not instantiated and will crash if the job distributor actually tries to use it.
  //That means that this can only be used with parser=true
  option[ jd2::dd_parser ].value( true ); // So here we fix that. jd2_parser app makes no sense without this option=true
  if( !option[ jd2::ntrials ].user() ) {
    // when using rosteta_scripts we want ntrials to be set to 1 if the user forgot to specify. o/w infinite loops might
    // occur. We don't want ntrials to be set as default to 1, b/c other protocols might want it to behave differently
    option[ jd2::ntrials ].value( 1 );
  }

  try {
    if( view )
      protocols::viewer::viewer_main( my_main );
    else{
      protocols::jd2::JobDistributor::get_instance()->go( mover );
    }
  } catch( utility::excn::EXCN_Base& excn ) {
    basic::Error()
      << "ERROR: Exception caught by rosetta_scripts application:"
      << excn << std::endl;
    assert(false); // core dump in debug mode
    std::exit( 1 );
  }
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl; 
	} 
}
