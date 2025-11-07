// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/bootcamp/BootCampMover.cc
/// @brief Lab6 Mover
/// @author Noora Azadvari (eazadvar@uoregon.edu)

// Unit headers
#include <protocols/bootcamp/BootCampMover.hh>
#include <protocols/bootcamp/BootCampMoverCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <numeric/random/random.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/dssp/Dssp.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/memory.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// Citation Manager
#include <utility/vector1.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

static basic::Tracer TR( "protocols.bootcamp.BootCampMover" );

namespace protocols {
namespace bootcamp {

	/////////////////////
	/// Constructors  ///
	/////////////////////

/// @brief Default constructor
BootCampMover::BootCampMover():
	protocols::moves::Mover( BootCampMover::mover_name() )
{
	std::cout << "BootCampMover constructed!\n"; 
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
BootCampMover::~BootCampMover(){}

////////////////////////////////////////////////////////////////////////////////
	/// Mover Methods ///
	/////////////////////
// Secondary structure -> vector of edges for non-loop spans
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
		ss_boundaries.push_back( std::make_pair( strand_start+1, ss_string.size() ));
	}

	for ( core::Size ii = 1; ii <= ss_boundaries.size(); ++ii ) {
		std::cout << "SS Element " << ii << " from residue "
		<< ss_boundaries[ ii ].first << " to "
		<< ss_boundaries[ ii ].second << std::endl;
	}

	return ss_boundaries;
} // END of function : identify_secondary_structure_spans

// Edges -> foldtree
core::kinematics::FoldTree
fold_tree_from_edges( utility::vector1< std::pair< core::Size, core::Size > > const & edges ) {
	core::kinematics::FoldTree foldtree;
  // Midpoint of the first non-loop span
  int mid_first = static_cast<int>( (edges[1].first + edges[1].second) / 2  );
  // Creating vector of loop edges
	utility::vector1< std::pair< core::Size, core::Size > > gaps;
  for (auto e = edges.begin(); e!= edges.end() -1 ; e++){
			auto next_e = e+1;
			bool there_is_a_gap = static_cast<bool>( next_e->first - e->second -1  );
			if ( there_is_a_gap ){ gaps.push_back({e->second+1,next_e->first-1}); }
	}
	// Adding jumps
  // mid-first to all other mids (non-loops, then loops)
  int n = 1;
	for (auto e = edges.begin() + 1; e!= edges.end(); e++){
		int mid_e = static_cast<int>( (e->first + e->second) / 2);
		foldtree.add_edge(mid_first, mid_e, n);
		n++;
	}
  for (auto e = gaps.begin()  + 1; e!= gaps.end(); e++){
		int mid_e = static_cast<int>( (e->first + e->second) / 2);
		foldtree.add_edge(mid_first, mid_e, n);
		n++;
	}

  // Adding peptide edges
  // mid to both ways (non-loops, then loops)
  for (auto e = edges.begin(); e!= edges.end(); e++){
		int mid_e = static_cast<int>( (e->first + e->second) / 2);
		foldtree.add_edge(mid_e, e->first , core::kinematics::Edge::PEPTIDE);
		foldtree.add_edge(mid_e, e->second, core::kinematics::Edge::PEPTIDE);
	}
	for (auto e = gaps.begin(); e!= gaps.end(); e++){
		int mid_e = static_cast<int>( (e->first + e->second) / 2);
		foldtree.add_edge(mid_e, e->first , core::kinematics::Edge::PEPTIDE);
		foldtree.add_edge(mid_e, e->second, core::kinematics::Edge::PEPTIDE);
	}

return foldtree;

} // END of function : fold_tree_from_edges

/// @brief Apply the mover
void
BootCampMover::apply( core::pose::Pose& pose ){

	// Setting up score functions and movers
	core::scoring::ScoreFunctionCOP sfxn = core::scoring::get_score_function();
  protocols::moves::MonteCarloOP mover = utility::pointer::make_shared<protocols::moves::MonteCarlo>( pose, *sfxn, 0.5 );

	// Pose -> ss
	core::scoring::dssp::Dssp dssp = core::scoring::dssp::Dssp( pose );
  std::string secondary_structure = dssp.get_dssp_secstruct();
	// ss -> edges
  utility::vector1< std::pair< core::Size, core::Size > > edges_nonloop = 
  identify_secondary_structure_spans( secondary_structure ); 
	// edges -> foldtree
	core::kinematics::FoldTree ft = fold_tree_from_edges( edges_nonloop ); 
	// set ft to pose
	pose.fold_tree( ft );
  // purturb phi-psi angles
  double N = static_cast<double>( pose.size() );
  core::Real random_uni = numeric::random::uniform();
  core::Size random_res_num = static_cast< core::Size > ( random_uni * N + 1 );
  core::Real random_pert1 = numeric::random::gaussian();
  core::Real random_pert2 = numeric::random::gaussian();
  core::Real orig_phi = pose.phi( random_res_num );
  core::Real orig_psi = pose.psi( random_res_num );
  pose.set_phi( random_res_num, orig_phi + random_pert1 );
  pose.set_psi( random_res_num, orig_psi + random_pert2 );
  // monte carlo boltzmann minimizer
  mover->boltzmann( pose );
	

}
////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
BootCampMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
BootCampMover::parse_my_tag(
	utility::tag::TagCOP ,
	basic::datacache::DataMap&
) {

}
void BootCampMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Lab6 Mover", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
BootCampMover::fresh_instance() const
{
	return utility::pointer::make_shared< BootCampMover >();
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
BootCampMover::clone() const
{
	return utility::pointer::make_shared< BootCampMover >( *this );
}

std::string BootCampMover::get_name() const {
	return mover_name();
}

std::string BootCampMover::mover_name() {
	return "BootCampMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
BootCampMoverCreator::create_mover() const
{
	return utility::pointer::make_shared< BootCampMover >();
}

std::string
BootCampMoverCreator::keyname() const
{
	return BootCampMover::mover_name();
}

void BootCampMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	BootCampMover::provide_xml_schema( xsd );
}

/// @brief This mover is unpublished.  It returns Noora Azadvari as its author.
void
BootCampMover::provide_citation_info(basic::citation_manager::CitationCollectionList & citations ) const {
	citations.add(
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		"BootCampMover", basic::citation_manager::CitedModuleType::Mover,
		"Noora Azadvari",
		"TODO: institution",
		"eazadvar@uoregon.edu",
		"Wrote the BootCampMover."
		)
	);
}


////////////////////////////////////////////////////////////////////////////////
	/// private methods ///
	///////////////////////


std::ostream &
operator<<( std::ostream & os, BootCampMover const & mover )
{
	mover.show(os);
	return os;
}


} //bootcamp
} //protocols
