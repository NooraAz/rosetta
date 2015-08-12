// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    interface_statistics.cc
/// @brief   Computing interface statistics for membrane protein interfaces
/// @details last Modified: 5/14/15
/// @author  JKLeman (julia.koehler1982@gmail.com)

// App headers
#include <devel/init.hh>

// Project headers
#include <core/conformation/membrane/MembraneInfo.hh>
#include <protocols/scoring/Interface.fwd.hh>
#include <protocols/scoring/Interface.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <protocols/moves/Mover.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <protocols/toolbox/pose_metric_calculators/PiPiCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

// Package Headers
#include <core/scoring/sasa/SasaCalc.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <core/pose/util.hh>
#include <utility/string_util.hh>
#include <core/scoring/sasa/util.hh>


static basic::Tracer TR( "apps.pilot.jkleman.interface_statistics" );

using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;
using namespace protocols::moves;
using namespace protocols::scoring;
using namespace protocols::toolbox::pose_metric_calculators;

class MPInterfaceStatistics : public Mover {
	
public:
	
	////////////////////
	/// Constructors ///
	////////////////////

	/// @brief Construct a Default Membrane Position Mover
	MPInterfaceStatistics( Size jump );
	
	/// @brief Copy Constructor
	/// @details Make a deep copy of this mover object
	MPInterfaceStatistics( MPInterfaceStatistics const & src );
	
	/// @brief Assignment Operator
	/// @details Make a deep copy of this mover object, overriding the assignment operator
	MPInterfaceStatistics &
	operator=( MPInterfaceStatistics const & src );
	
	/// @brief Destructor
	~MPInterfaceStatistics();

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this mover
	virtual std::string get_name() const;

	/// @brief Get Membrane protein interface statistics
	virtual void apply( Pose & pose );

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

//	/// @brief Create a Clone of this mover
//	virtual protocols::moves::MoverOP clone() const;
//
//	/// @brief Create a Fresh Instance of this Mover
//	virtual protocols::moves::MoverOP fresh_instance() const;
//
//	/// @brief Pase Rosetta Scripts Options for this Mover
//	void parse_my_tag(
//					  utility::tag::TagCOP tag,
//					  basic::datacache::DataMap &,
//					  protocols::filters::Filters_map const &,
//					  protocols::moves::Movers_map const &,
//					  core::pose::Pose const &
//					  );
	
private: // methods

	/// @brief Register Options with JD2
	void register_options();

	/// @brief Initialize Mover options from the comandline
	void init_from_cmd();

	/// @brief Get size of the interface
	Real get_size( Pose & pose );

	/// @brief Get charges for all states
	utility::vector1< SSize > get_charge( Pose & pose );

	/// @brief Get charges for all states
	utility::vector1< Real > get_avg_charge( Pose & pose );

		/// @brief Get number of charged residues for all states
	utility::vector1< Size > get_number_charges( Pose & pose );

	/// @brief Get number of charged residues for all states
	utility::vector1< Real > get_avg_number_charges( Pose & pose );

	/// @brief Get hydrophobicity for all states (UHS, Koehler & Meiler, Proteins, 2008)
	utility::vector1< Real > get_hydrophobicity( Pose & pose );

	/// @brief Get average hydrophobicity for all states (UHS, Koehler & Meiler, Proteins, 2008)
	utility::vector1< Real > get_avg_hydrophobicity( Pose & pose );

	/// @brief Get number of contacts
	utility::vector1< Size > get_number_of_residues( Pose & pose );

	/// @brief Get surface complementarity
	Real get_surface_complementarity( Pose & pose );

	/// @brief Get size of residues for all states
	utility::vector1< Real > get_residue_size( Pose & pose );

	/// @brief Get number of Hbonds across interface
	Size get_number_hbonds( Pose & pose );

	/// @brief Get number of pi-stacking across interface
	Size get_number_pi_stacking( Pose & pose );

	/// @brief Get in/out orientation
	Real get_inout_orientation( Pose & pose );

	/// @brief Add output to scorefile
	void add_to_scorefile( Pose & pose );

	/// @brief fill vectors with data
	void fill_vectors_with_data( Pose & pose );

private: // data
	
	/// @brief jump
	Size jump_;
	
	/// @brief Interface
	Interface interface_;
	
	/// @brief Partners;
	std::string partners_;
	
	/// @brief Are residues in the membrane?
	utility::vector1< bool > in_mem_;
	
//	/// @brief What is the charge of the residues?
//	utility::vector1< bool > charge_;
	
	/// @brief Are the residues in the interface?
	utility::vector1< bool > intf_;
	
	/// @brief Are the residues exposed?
	utility::vector1< bool > exposed_;
};

//////////////////////////////////////////////////////////////////////

////////////////////
/// Constructors ///
////////////////////

/// @brief Construct a Default Membrane Position Mover
MPInterfaceStatistics::MPInterfaceStatistics( Size jump ) :
	Mover(),
	jump_( jump ),
	interface_(),
	partners_()
{}

/// @brief Copy Constructor
/// @details Make a deep copy of this mover object
MPInterfaceStatistics::MPInterfaceStatistics( MPInterfaceStatistics const & src ) :
	Mover( src ),
	jump_( src.jump_ ),
	interface_( src.interface_ ),
	partners_( src.partners_ )
{}

/// @brief Assignment Operator
/// @details Make a deep copy of this mover object, overriding the assignment operator
MPInterfaceStatistics &
MPInterfaceStatistics::operator=( MPInterfaceStatistics const & src )
{
	// Abort self-assignment.
	if (this == &src) {
	return *this;
	}

	// Otherwise, create a new object
	return *( new MPInterfaceStatistics( *this ) );
}

/// @brief Destructor
MPInterfaceStatistics::~MPInterfaceStatistics() {}

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this mover
std::string
MPInterfaceStatistics::get_name() const {
	return "MPInterfaceStatistics";
}

/// @brief Get Membrane protein interface statistics
void MPInterfaceStatistics::apply( Pose & pose ) {
	
	using namespace protocols::membrane;
	
	TR << "Calling MPInterfaceStatistics" << std::endl;

	register_options();
	init_from_cmd();

	// add membrane
	AddMembraneMoverOP addmem( new AddMembraneMover() );
	addmem->apply( pose );

	// Check the pose is a membrane protein
	if (! pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Cannot apply membrane move to a non-membrane pose!" );
	}

	// create scorefunction and score the pose
	core::scoring::ScoreFunctionOP highres_sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012.wts" );
	( *highres_sfxn )( pose );

	// calculate interface
	interface_ = Interface( jump_ );
	interface_.calculate( pose );
	
	// fill vectors with data
	fill_vectors_with_data( pose );

	// add measures to scorefile
	add_to_scorefile( pose );

	// other ideas:
	// Compute a ruggedness factor of the interface???
	// Compute an average deviation from a plane along the interface
	// Get contact preferences for AAs: Yan & Jernigan, 2008, ProteinJ
	
} // apply

/// @brief Register Options with JD2
void MPInterfaceStatistics::register_options() {
	
	using namespace basic::options;
	option.add_relevant( OptionKeys::docking::partners );

} // register options

/// @brief Initialize Mover options from the comandline
void MPInterfaceStatistics::init_from_cmd() {
	
	using namespace basic::options;
	
	// docking partners
	if( option[ OptionKeys::docking::partners ].user() ) {
		partners_ = option[ OptionKeys::docking::partners ]();
	}
} // init from commandline

/////////////////////////////////
///// Rosetta Scripts Methods ///
/////////////////////////////////
//
///// @brief Create a Clone of this mover
//protocols::moves::MoverOP
//MPInterfaceStatistics::clone() const {
//	return ( protocols::moves::MoverOP( new MPInterfaceStatistics( *this ) ) );
//}
//
///// @brief Create a Fresh Instance of this Mover
//protocols::moves::MoverOP
//MPInterfaceStatistics::fresh_instance() const {
//	return protocols::moves::MoverOP( new MPInterfaceStatistics() );
//}
//
///// @brief Pase Rosetta Scripts Options for this Mover
//void
//MPInterfaceStatistics::parse_my_tag(
//									   utility::tag::TagCOP tag,
//									   basic::datacache::DataMap &,
//									   protocols::filters::Filters_map const &,
//									   protocols::moves::Movers_map const &,
//									   core::pose::Pose const &
//									   ) {
//		// TODO
//
//}
//
///// @brief Create a new copy of this mover
//protocols::moves::MoverOP
//MPInterfaceStatisticsCreator::create_mover() const {
//	return protocols::moves::MoverOP( new MPInterfaceStatistics );
//}
//
///// @brief Return the Name of this mover (as seen by Rscripts)
//std::string
//MPInterfaceStatisticsCreator::keyname() const {
//	return MPInterfaceStatisticsCreator::mover_name();
//}
//
///// @brief Mover name for Rosetta Scripts
//std::string
//MPInterfaceStatisticsCreator::mover_name() {
//	return "MPInterfaceStatistics";
//}

//////////////////////////////////////////////////////////////////////

typedef utility::pointer::shared_ptr< MPInterfaceStatistics > MPInterfaceStatisticsOP;

//////////////////////////////////////////////////////////////////////

/// @brief Get size of the interface in square Angstrom
Real MPInterfaceStatistics::get_size( Pose & pose ) {

	using namespace numeric;
	using namespace core::pose;
	using namespace core::scoring::sasa;
	using namespace protocols::scoring;
	
	// get per-residue SASA
	SasaCalc calc = SasaCalc();
	core::Real total_sasa( calc.calculate( pose ) );
	
	// partition pose by jump
	Pose partner1, partner2;
	partition_pose_by_jump( pose, jump_, partner1, partner2 );

	// calculate SASAs of individual partners
	core::Real p1_sasa( calc.calculate( partner1 ) );
	core::Real p2_sasa( calc.calculate( partner2 ) );
	
	// calculate interface
	core::Real intf_sasa = - ( total_sasa - p1_sasa - p2_sasa) / 2;

	return intf_sasa;

}// get_size

	//////////////////////////////////////////////////////////////////////

/// @brief Get charge of the interface
utility::vector1< SSize > MPInterfaceStatistics::get_charge( Pose & pose ) {

	SSize tm_int( 0 );
	SSize tm_nonint( 0 );
	SSize tm_core( 0 );
	SSize sol_int( 0 );
	SSize sol_nonint( 0 );
	SSize sol_core( 0 );

	// go through residues
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {
	
		std::string res = pose.residue( i ).name3();
		
		// is residue positive?
		bool pos( false );
		if ( res == "ARG" || res == "LYS" ) {
			pos = true;
		}

		// is residue negative?
		bool neg( false );
		if ( res == "ASP" || res == "GLU" ) {
			neg = true;
		}

		// interface, membrane, positive / negative
		if ( exposed_[i] && intf_[i] && in_mem_[i] && pos ) {
			++tm_int;
		}
		else if ( exposed_[i] && intf_[i] && in_mem_[i] && neg ) {
			--tm_int;
		}
		
		// interface, sol, positive / negative
		else if ( exposed_[i] && intf_[i] && !in_mem_[i] && pos ) {
			++sol_int;
		}
		else if ( exposed_[i] && intf_[i] && !in_mem_[i] && neg ) {
			--sol_int;
		}

		// not interface, membrane, positive / negative
		else if ( exposed_[i] && !intf_[i] && in_mem_[i] && pos ) {
			++tm_nonint;
		}
		else if ( exposed_[i] && !intf_[i] && in_mem_[i] && neg ) {
			--tm_nonint;
		}

		// not interface, sol, positive / negative
		else if ( exposed_[i] && !intf_[i] && !in_mem_[i] && pos ) {
			++sol_nonint;
		}
		else if ( exposed_[i] && !intf_[i] && !in_mem_[i] && neg ) {
			--sol_nonint;
		}
		
		// if core, membrane, positive / negative
		else if ( !exposed_[i] && in_mem_[i] && pos ) {
			++tm_core;
		}
		else if ( !exposed_[i] && in_mem_[i] && neg ) {
			--tm_core;
		}

		// if core, solution, positive / negative
		else if ( !exposed_[i] && !in_mem_[i] && pos ) {
			++sol_core;
		}
		else if ( !exposed_[i] && !in_mem_[i] && neg ) {
			--sol_core;
		}
		
	}

	// write counts into vector
	utility::vector1< SSize > charge;
	charge.push_back( tm_int );
	charge.push_back( tm_nonint );
	charge.push_back( tm_core );
	charge.push_back( sol_int );
	charge.push_back( sol_nonint );
	charge.push_back( sol_core );

	return charge;

} // get charge

	//////////////////////////////////////////////////////////////////////

/// @brief Get average charge of the interface, meaning averaged over the
///			number of residues
utility::vector1< Real > MPInterfaceStatistics::get_avg_charge( Pose & pose ) {

	utility::vector1< SSize > charge( get_charge( pose ) );
	utility::vector1< Size > nres( get_number_of_residues( pose ) );
	utility::vector1< Real > avg_charge;
	
	// get average
	for ( Size i = 1; i <= nres.size(); ++i ) {
		avg_charge.push_back( static_cast< Real >( charge[ i ] ) / static_cast< Real >( nres[ i ] ) );
	}
	
	return avg_charge;

} // get avg charge

	//////////////////////////////////////////////////////////////////////

/// @brief Get number of charged residues in the interface
utility::vector1< Size > MPInterfaceStatistics::get_number_charges( Pose & pose ) {

	Size tm_int( 0 );
	Size tm_nonint( 0 );
	Size tm_core( 0 );
	Size sol_int( 0 );
	Size sol_nonint( 0 );
	Size sol_core( 0 );
	
	// go through residues
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {
		
		std::string res = pose.residue( i ).name3();

		// is residue charged?
		bool charged( false );
		if ( res == "ARG" || res == "LYS" || res == "ASP" || res == "GLU" ) {
			charged = true;
		}
		
		// interface, membrane, charged
		if ( exposed_[i] && intf_[i] && in_mem_[i] && charged ) {
			++tm_int;
		}
		
		// interface, sol, charged
		else if ( exposed_[i] && intf_[i] && !in_mem_[i] && charged ) {
			++sol_int;
		}
		
		// not interface, membrane, charged
		else if ( exposed_[i] && !intf_[i] && in_mem_[i] && charged ) {
			++tm_nonint;
		}
	
		// not interface, sol, charged
		else if ( exposed_[i] && !intf_[i] && !in_mem_[i] && charged ) {
			++sol_nonint;
		}
		
		// core, membrane
		else if ( !exposed_[i] && in_mem_[i] && charged ) {
			++tm_core;
		}
		
		// core, membrane
		else if ( !exposed_[i] && !in_mem_[i] && charged ) {
			++sol_core;
		}
		
	}
	
	// write counts into vector
	utility::vector1< Size > ncharges;
	ncharges.push_back( tm_int );
	ncharges.push_back( tm_nonint );
	ncharges.push_back( tm_core );
	ncharges.push_back( sol_int );
	ncharges.push_back( sol_nonint );
	ncharges.push_back( sol_core );
	
	return ncharges;

} // get number of charges

	//////////////////////////////////////////////////////////////////////

/// @brief Get average number of charges, meaning averaged over the
///			number of residues
utility::vector1< Real > MPInterfaceStatistics::get_avg_number_charges( Pose & pose ) {

	utility::vector1< Size > ncharges( get_number_charges( pose ) );
	utility::vector1< Size > nres( get_number_of_residues( pose ) );
	utility::vector1< Real > avg_ncharges;

	// get average
	for ( Size i = 1; i <= nres.size(); ++i ) {
		avg_ncharges.push_back( static_cast< Real >( ncharges[ i ] ) / static_cast< Real >( nres[ i ] ) );
	}

	return avg_ncharges;

} // get avg charge

	//////////////////////////////////////////////////////////////////////

/// @brief Get hydrophobicity of the interface (UHS, Koehler & Meiler, Proteins, 2008)
utility::vector1< Real > MPInterfaceStatistics::get_hydrophobicity( Pose & pose ) {

	// hydrophobicity scale
	std::map< char, Real > hydro;
	
	// polar
	hydro[ 'C' ] = 0.01;
	hydro[ 'N' ] = 0.50;
	hydro[ 'Q' ] = 0.46;
	hydro[ 'S' ] = 0.06;
	hydro[ 'T' ] = -0.01;
	
	// charged
	hydro[ 'D' ] = 0.73;
	hydro[ 'E' ] = 0.70;
	hydro[ 'K' ] = 0.90;
	hydro[ 'R' ] = 0.55;
	
	// hydrophobic
	hydro[ 'A' ] = -0.16;
	hydro[ 'G' ] = -0.20;
	hydro[ 'I' ] = -0.39;
	hydro[ 'L' ] = -0.30;
	hydro[ 'M' ] = -0.20;
	hydro[ 'P' ] = 0.50;
	hydro[ 'V' ] = -0.25;
	
	// aromatic
	hydro[ 'F' ] = -0.46;
	hydro[ 'H' ] = 0.38;
	hydro[ 'W' ] = -0.03;
	hydro[ 'Y' ] = -0.12;

	Real tm_int( 0 );
	Real tm_nonint( 0 );
	Real tm_core( 0 );
	Real sol_int( 0 );
	Real sol_nonint( 0 );
	Real sol_core( 0 );
	
	// go through residues
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {

		// interface, mem
		if ( exposed_[i] && intf_[i] && in_mem_[i] ) {
			tm_int += hydro[ pose.residue( i ).name1() ];
		}

		// interface, sol
		else if ( exposed_[i] && intf_[i] && !in_mem_[i] ) {
			sol_int += hydro[ pose.residue( i ).name1() ];
		}

		// not interface, mem
		else if ( exposed_[i] && !intf_[i] && in_mem_[i] ) {
			tm_nonint += hydro[ pose.residue( i ).name1() ];
		}

		// not interface, sol
		else if ( exposed_[i] && !intf_[i] && !in_mem_[i] ) {
			sol_nonint += hydro[ pose.residue( i ).name1() ];
		}
		
		// core, mem
		else if ( !exposed_[i] && in_mem_[i] ) {
			tm_core += hydro[ pose.residue( i ).name1() ];
		}
		
		// core, sol
		else if ( !exposed_[i] && !in_mem_[i] ) {
			sol_core += hydro[ pose.residue( i ).name1() ];
		}
	}

	// write counts into vector
	utility::vector1< Real > total_hydro;
	total_hydro.push_back( tm_int );
	total_hydro.push_back( tm_nonint );
	total_hydro.push_back( tm_core );
	total_hydro.push_back( sol_int );
	total_hydro.push_back( sol_nonint );
	total_hydro.push_back( sol_core );

	return total_hydro;

} // get hydrophobicity

	//////////////////////////////////////////////////////////////////////

/// @brief Get avg hydrophobicity for all states (UHS, Koehler & Meiler, Proteins, 2008)
utility::vector1< Real > MPInterfaceStatistics::get_avg_hydrophobicity( Pose & pose ) {

	utility::vector1< Real > hydro( get_hydrophobicity( pose ) );
	utility::vector1< Size > nres( get_number_of_residues( pose ) );
	utility::vector1< Real > avg_hydro;
	
	// get average
	for ( Size i = 1; i <= nres.size(); ++i ) {
		avg_hydro.push_back( hydro[ i ] / nres[ i ] );
	}
	
	return avg_hydro;

} // average hydrophobicity

	//////////////////////////////////////////////////////////////////////
	

/// @brief Get number of residues in each bin
utility::vector1< Size > MPInterfaceStatistics::get_number_of_residues( Pose & pose ) {

	Size tm_int( 0 );
	Size tm_nonint( 0 );
	Size tm_core( 0 );
	Size sol_int( 0 );
	Size sol_nonint( 0 );
	Size sol_core( 0 );
	
	// go through residues
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {
		
		// interface, mem
		if ( exposed_[i] && intf_[i] && in_mem_[i] ) {
			++tm_int;
		}

		// interface, sol
		else if ( exposed_[i] && intf_[i] && !in_mem_[i] ) {
			++sol_int;
		}

		// non interface, mem
		else if ( exposed_[i] && !intf_[i] && in_mem_[i] ) {
			++tm_nonint;
		}
		
		// not interface, sol
		else if ( exposed_[i] && !intf_[i] && !in_mem_[i] ) {
			++sol_nonint;
		}

		// core, mem
		else if ( !exposed_[i] && in_mem_[i] ) {
			++tm_core;
		}
		
		// core, sol
		else if ( !exposed_[i] && !in_mem_[i] ) {
			++sol_core;
		}
	}
	
	// write counts into vector
	utility::vector1< Size > nresidues;
	nresidues.push_back( tm_int );
	nresidues.push_back( tm_nonint );
	nresidues.push_back( tm_core );
	nresidues.push_back( sol_int );
	nresidues.push_back( sol_nonint );
	nresidues.push_back( sol_core );

	return nresidues;
	
} // get number of contacts in interface

	//////////////////////////////////////////////////////////////////////

/// @brief Get surface complementarity
Real MPInterfaceStatistics::get_surface_complementarity( Pose & pose ) {

	using namespace core::scoring::sc;
	ShapeComplementarityCalculator shape;
	shape.Init();
	
	return shape.CalcSc( pose, jump_, 1 );

} // get surface complementarity of interface

	//////////////////////////////////////////////////////////////////////

/// @brief Get avg residue in the interface
utility::vector1 < Real > MPInterfaceStatistics::get_residue_size( Pose & pose ) {

	Real tm_int( 0 );
	Real tm_nonint( 0 );
	Real tm_core( 0 );
	Real sol_int( 0 );
	Real sol_nonint( 0 );
	Real sol_core( 0 );

	// go through residues
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {

		// interface, mem
		if ( exposed_[i] && intf_[i] && in_mem_[i] ) {
			tm_int += static_cast< Real >( pose.residue( i ).natoms() );
		}
		
		// interface, sol
		else if ( exposed_[i] && intf_[i] && !in_mem_[i] ) {
			sol_int += static_cast< Real >( pose.residue( i ).natoms() );
		}

		// not interface, mem
		else if ( exposed_[i] && !intf_[i] && in_mem_[i] ) {
			tm_nonint += static_cast< Real >( pose.residue( i ).natoms() );
		}
		
		// not interface, sol
		else if ( exposed_[i] && !intf_[i] && !in_mem_[i] ) {
			sol_nonint += static_cast< Real >( pose.residue( i ).natoms() );
		}

		// core, mem
		else if ( !exposed_[i] && in_mem_[i] ) {
			tm_core += static_cast< Real >( pose.residue( i ).natoms() );
		}
		
		// core, sol
		else if ( !exposed_[i] && !in_mem_[i] ) {
			sol_core += static_cast< Real >( pose.residue( i ).natoms() );
		}
	}
	
	utility::vector1< Size > nres( get_number_of_residues( pose ) );

	TR << "number of residues:" << std::endl;
	TR << "tm int " << nres[1] << std::endl;
	TR << "tm nonint " << nres[2] << std::endl;
	TR << "tm core " << nres[3] << std::endl;
	TR << "sol int " << nres[4] << std::endl;
	TR << "sol nonint " << nres[5] << std::endl;
	TR << "sol core " << nres[6] << std::endl;
	
	TR << "residue size:" << std::endl;
	TR << "tm int " << tm_int << std::endl;
	TR << "tm nonint " << tm_nonint << std::endl;
	TR << "tm core " << tm_core << std::endl;
	TR << "sol int " << sol_int << std::endl;
	TR << "sol nonint " << sol_nonint << std::endl;
	TR << "sol core " << sol_core << std::endl;

	// get average residue size
	tm_int /= static_cast< Real > ( nres[1] );
	tm_nonint /= static_cast< Real > ( nres[2] );
	tm_core /= static_cast< Real > ( nres[3] );
	sol_int /= static_cast< Real > ( nres[4] );
	sol_nonint /= static_cast< Real > ( nres[5] );
	sol_core /= static_cast< Real > ( nres[6] );

	// write counts into vector
	utility::vector1< Real > avg_size;
	avg_size.push_back( tm_int );
	avg_size.push_back( tm_nonint );
	avg_size.push_back( tm_core );
	avg_size.push_back( sol_int );
	avg_size.push_back( sol_nonint );
	avg_size.push_back( sol_core );
	
	return avg_size;
	
} // get residue size

	//////////////////////////////////////////////////////////////////////

/// @brief Get number of Hbonds across interface
Size MPInterfaceStatistics::get_number_hbonds( Pose & pose ) {

	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::scoring::hbonds;
	
	// get the partners
	utility::vector1< std::string > partners( utility::string_split( partners_, '_' ) );
	utility::vector1< std::string > partner1( utility::split( partners[1] ) );
	utility::vector1< std::string > partner2( utility::split( partners[2] ) );

	// get the Hbonds
	hbonds::HBondDatabaseCOP hb_database(hbonds::HBondDatabase::get_database());
	hbonds::HBondSet hb_set;
	pose.update_residue_neighbors();
	hbonds::fill_hbond_set( pose, false, hb_set, false );
//	hb_set.show(pose, true, std::cout);
	
	Size nhbonds( 0 );
	
	// go through Hbonds
	for ( Size i = 1; i <= hb_set.nhbonds(); ++i ) {

		HBond hbond = hb_set.hbond( i );
		Size don = hbond.don_res();
		Size acc = hbond.acc_res();

		// check whether donor/acc is in partner1
		// doesn't check whether the residues are defined as interface (they should anyway because they are between the two partners)
		
		bool don_in_1( false );
		bool acc_in_1( false );
		
		// donor in partner1?
		for ( Size j = 1; j <= partner1.size(); ++j ) {
			don_in_1 = res_in_chain( pose, don, partner1[ j ] );
		}

		// acceptor in partner1?
		for ( Size j = 1; j <= partner1.size(); ++j ) {
			acc_in_1 = res_in_chain( pose, acc, partner1[ j ] );
		}
		
//		TR << "hbond: don " << don << ", acc " << acc;
//		TR << ", don: " << don_in_1 << ", acc: " << acc_in_1 << std::endl;

		// if donor in 1 and acceptor not, then we have Hbond across interface
		if ( don_in_1 == true && acc_in_1 == false ) {
			++nhbonds;
		}
		// similarly, if donor not in 1 but acc is, we also have Hbond across interface
		else if ( don_in_1 == false && acc_in_1 == true ) {
			++nhbonds;
		}
	} // go through Hbonds
	

	return nhbonds;
} // number hbonds across interface

	//////////////////////////////////////////////////////////////////////

/// @brief Get number of pi-stacking across interface
Size MPInterfaceStatistics::get_number_pi_stacking( Pose & pose ) {

	core::pose::metrics::PoseMetricCalculatorOP pi_pi_calculator ( new protocols::toolbox::pose_metric_calculators::PiPiCalculator() );
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "pi_pi", pi_pi_calculator );

	TR << "WARNING: PLEASE CHECK THE NUMBER OF YOUR Pi-stacking INTERACTIONS VISUALLY TO CONFIRM!" << std::endl;

	return utility::string2Size( pose.print_metric( "pi_pi", "pi_pi" ) );
}

	//////////////////////////////////////////////////////////////////////

/// @brief Get in/out orientation
Real MPInterfaceStatistics::get_inout_orientation( Pose & pose ) {
	
	if ( pose.residue( nres_protein( pose ) ).atom("CA").xyz().z() > 0.0 ){
		return 0.0;
	}
	return 1.0;

}// get in/out orientation

	//////////////////////////////////////////////////////////////////////

/// @brief Add output to scorefile
void MPInterfaceStatistics::add_to_scorefile( Pose & pose ) {
	
	// get job
	protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );
	
	// get size of the interface
	Real size = get_size( pose );
	job->add_string_real_pair( "Isize", size );
	TR << "intf size " << size << std::endl;

	// Get surface complementarity
	Real shape = get_surface_complementarity( pose );
	job->add_string_real_pair( "shape", shape );
	TR << "shape_compl " << shape << std::endl;

	// Get number of Hbonds across interface
	Size nhbonds = get_number_hbonds( pose );
	job->add_string_real_pair( "nhbonds", nhbonds );
	TR << "nhbonds " << nhbonds << std::endl;
	
	// Get number of pi-stacking across interface
	Size npi = get_number_pi_stacking( pose );
	job->add_string_real_pair( "npi", npi );
	TR << "pi stacking " << npi << std::endl;
	
	// Get in/out orientation
	Real inout = get_inout_orientation( pose );
	job->add_string_real_pair( "inout", inout );
	TR << "inout " << inout << std::endl;
	
	// Get charge in the interface
	utility::vector1< SSize > charge = get_charge( pose );
	job->add_string_real_pair( "chrg_tm_int", charge[1] );
	job->add_string_real_pair( "chrg_tm_nint", charge[2] );
	job->add_string_real_pair( "chrg_tm_core", charge[3] );
	job->add_string_real_pair( "chrg_so_int", charge[4] );
	job->add_string_real_pair( "chrg_so_nint", charge[5] );
	job->add_string_real_pair( "chrg_so_core", charge[6] );
	TR << "charge tm int " << charge[1] << std::endl;
	TR << "charge tm nonint " << charge[2] << std::endl;
	TR << "charge tm core " << charge[3] << std::endl;
	TR << "charge sol int " << charge[4] << std::endl;
	TR << "charge sol nonint " << charge[5] << std::endl;
	TR << "charge sol core " << charge[6] << std::endl;

	// Get avg charge in the interface
	utility::vector1< Real > avg_charge = get_avg_charge( pose );
	job->add_string_real_pair( "avg_chrg_tm_int", avg_charge[1] );
	job->add_string_real_pair( "avg_chrg_tm_nint", avg_charge[2] );
	job->add_string_real_pair( "avg_chrg_tm_core", avg_charge[3] );
	job->add_string_real_pair( "avg_chrg_so_int", avg_charge[4] );
	job->add_string_real_pair( "avg_chrg_so_nint", avg_charge[5] );
	job->add_string_real_pair( "avg_chrg_so_core", avg_charge[6] );
	TR << "avg charge tm int " << avg_charge[1] << std::endl;
	TR << "avg charge tm nonint " << avg_charge[2] << std::endl;
	TR << "avg charge tm core " << avg_charge[3] << std::endl;
	TR << "avg charge sol int " << avg_charge[4] << std::endl;
	TR << "avg charge sol nonint " << avg_charge[5] << std::endl;
	TR << "avg charge sol core " << avg_charge[6] << std::endl;
	
	// Get number of charged residues in the interface
	utility::vector1< Size > ncharges = get_number_charges( pose );
	job->add_string_real_pair( "nchrg_tm_int", ncharges[1] );
	job->add_string_real_pair( "nchrg_tm_nint", ncharges[2] );
	job->add_string_real_pair( "nchrg_tm_core", ncharges[3] );
	job->add_string_real_pair( "nchrg_so_int", ncharges[4] );
	job->add_string_real_pair( "nchrg_so_nint", ncharges[5] );
	job->add_string_real_pair( "nchrg_so_core", ncharges[6] );
	TR << "ncharges tm int " << ncharges[1] << std::endl;
	TR << "ncharges tm nonint " << ncharges[2] << std::endl;
	TR << "ncharges tm core " << ncharges[3] << std::endl;
	TR << "ncharges sol int " << ncharges[4] << std::endl;
	TR << "ncharges sol nonint " << ncharges[5] << std::endl;
	TR << "ncharges sol core " << ncharges[6] << std::endl;

	// Get avg number of charged residues in the interface
	utility::vector1< Real > avg_ncharges = get_avg_number_charges( pose );
	job->add_string_real_pair( "avg_nchrg_tm_int", avg_ncharges[1] );
	job->add_string_real_pair( "avg_nchrg_tm_nint", avg_ncharges[2] );
	job->add_string_real_pair( "avg_nchrg_tm_core", avg_ncharges[3] );
	job->add_string_real_pair( "avg_nchrg_so_int", avg_ncharges[4] );
	job->add_string_real_pair( "avg_nchrg_so_nint", avg_ncharges[5] );
	job->add_string_real_pair( "avg_nchrg_so_core", avg_ncharges[6] );
	TR << "avg ncharges tm int " << avg_ncharges[1] << std::endl;
	TR << "avg ncharges tm nonint " << avg_ncharges[2] << std::endl;
	TR << "avg ncharges tm core " << avg_ncharges[3] << std::endl;
	TR << "avg ncharges sol int " << avg_ncharges[4] << std::endl;
	TR << "avg ncharges sol nonint " << avg_ncharges[5] << std::endl;
	TR << "avg ncharges sol core " << avg_ncharges[6] << std::endl;
	
	// Get hydrophobicity of the interface
	utility::vector1< Real > hydro = get_hydrophobicity( pose );
	job->add_string_real_pair( "hydro_tm_int", hydro[1] );
	job->add_string_real_pair( "hydro_tm_nint", hydro[2] );
	job->add_string_real_pair( "hydro_tm_core", hydro[3] );
	job->add_string_real_pair( "hydro_so_int", hydro[4] );
	job->add_string_real_pair( "hydro_so_nint", hydro[5] );
	job->add_string_real_pair( "hydro_so_core", hydro[6] );
	TR << "hydro tm int " << hydro[1] << std::endl;
	TR << "hydro tm nonint " << hydro[2] << std::endl;
	TR << "hydro tm core " << hydro[3] << std::endl;
	TR << "hydro sol int " << hydro[4] << std::endl;
	TR << "hydro sol nonint " << hydro[5] << std::endl;
	TR << "hydro sol core " << hydro[6] << std::endl;

	// Get average hydrophobicity of the interface
	utility::vector1< Real > avg_hydro = get_avg_hydrophobicity( pose );
	job->add_string_real_pair( "avg_hydro_tm_int", avg_hydro[1] );
	job->add_string_real_pair( "avg_hydro_tm_nint", avg_hydro[2] );
	job->add_string_real_pair( "avg_hydro_tm_core", avg_hydro[3] );
	job->add_string_real_pair( "avg_hydro_so_int", avg_hydro[4] );
	job->add_string_real_pair( "avg_hydro_so_nint", avg_hydro[5] );
	job->add_string_real_pair( "avg_hydro_so_core", avg_hydro[6] );
	TR << "avg hydro tm int " << avg_hydro[1] << std::endl;
	TR << "avg hydro tm nonint " << avg_hydro[2] << std::endl;
	TR << "avg hydro tm core " << avg_hydro[3] << std::endl;
	TR << "avg hydro sol int " << avg_hydro[4] << std::endl;
	TR << "avg hydro sol nonint " << hydro[5] << std::endl;
	TR << "avg hydro sol core " << avg_hydro[6] << std::endl;
	
	// Get number of contacts
	utility::vector1< Size > nresidues = get_number_of_residues( pose );
	job->add_string_real_pair( "nres_tm_int", nresidues[1] );
	job->add_string_real_pair( "nres_tm_nint", nresidues[2] );
	job->add_string_real_pair( "nres_tm_core", nresidues[3] );
	job->add_string_real_pair( "nres_so_int", nresidues[4] );
	job->add_string_real_pair( "nres_so_nint", nresidues[5] );
	job->add_string_real_pair( "nres_so_core", nresidues[6] );
	TR << "nresidues tm int " << nresidues[1] << std::endl;
	TR << "nresidues tm nonint " << nresidues[2] << std::endl;
	TR << "nresidues tm core " << nresidues[3] << std::endl;
	TR << "nresidues sol int " << nresidues[4] << std::endl;
	TR << "nresidues sol nonint " << nresidues[5] << std::endl;
	TR << "nresidues sol core " << nresidues[6] << std::endl;
	
	// Get size of residues in the interface
	utility::vector1< Real > res_size = get_residue_size( pose );
	job->add_string_real_pair( "ressize_tm_int", res_size[1] );
	job->add_string_real_pair( "ressize_tm_nint", res_size[2] );
	job->add_string_real_pair( "ressize_tm_core", res_size[3] );
	job->add_string_real_pair( "ressize_so_int", res_size[4] );
	job->add_string_real_pair( "ressize_so_nint", res_size[5] );
	job->add_string_real_pair( "ressize_so_core", res_size[6] );
	TR << "avg res size tm int " << res_size[1] << std::endl;
	TR << "avg res size tm nonint " << res_size[2] << std::endl;
	TR << "avg res size tm core " << res_size[3] << std::endl;
	TR << "avg res size sol int " << res_size[4] << std::endl;
	TR << "avg res size sol nonint " << res_size[5] << std::endl;
	TR << "avg res size sol core " << res_size[6] << std::endl;
	
} // add_to_scorefile

/// @brief Fill vectors with data
void MPInterfaceStatistics::fill_vectors_with_data( Pose & pose ) {

	using namespace core::scoring::sasa;

	// get relative per residue SASA
	utility::vector1< Real > sasa( rel_per_res_sc_sasa( pose ) );

	// go through residues
	for ( Size i = 1; i <= nres_protein( pose ); ++i ) {
		
		// is residue in the membrane?
		bool in_mem = pose.conformation().membrane_info()->in_membrane( i );
		in_mem_.push_back( in_mem );

		// is the residue in the interface?
		bool intf = interface_.is_interface( i );
		intf_.push_back( intf );

		// is the residue exposed?
		bool exposed( false );
		if ( sasa[i] >= 0.5 ) {
			exposed = true;
		}
		exposed_.push_back( exposed );
	}
		
} // fill vectors with data

////////////////////////////// MAIN ////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		using namespace protocols::jd2;
		
		// initialize options, RNG, and factory-registrators
		devel::init(argc, argv);
		
		MPInterfaceStatisticsOP mpis( new MPInterfaceStatistics( 1 ) );
		JobDistributor::get_instance()->go(mpis);
	}
	catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
