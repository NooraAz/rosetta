// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_O2StarMover
/// @brief O2Stars an RNA residue from a chain terminus.
/// @detailed
/// @author Rhiju Das

#include <protocols/swa/monte_carlo/RNA_O2StarMover.hh>

// libRosetta headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/VariantType.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>


static numeric::random::RandomGenerator RG(845099111);  // <- Magic number, do not change it!

using namespace core;
using core::Real;

//////////////////////////////////////////////////////////////////////////
// Removes one residue from a 5' or 3' chain terminus, and appropriately
// updates the pose full_model_info object.
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.swa.monte_carlo.o2star_mover" ) ;

namespace protocols {
namespace swa {
namespace monte_carlo {


  //////////////////////////////////////////////////////////////////////////
  //constructor!
	RNA_O2StarMover::RNA_O2StarMover( scoring::ScoreFunctionOP scorefxn ,
																		bool const sample_all_o2star,
																		Real const sample_range_small,
																		Real const sample_range_large ):
		scorefxn_( scorefxn ),
		sample_all_o2star_( sample_all_o2star ),
		sample_range_small_( sample_range_small ),
		sample_range_large_( sample_range_large )
  {}

  //////////////////////////////////////////////////////////////////////////
  //destructor
  RNA_O2StarMover::~RNA_O2StarMover()
  {}

  //////////////////////////////////////////////////////////////////////////
  void
  RNA_O2StarMover::apply( core::pose::Pose & pose )
	{
		std::string move_type = "";
		apply( pose, move_type );
	}

	//////////////////////////////////////////////////////////////////////
  void
  RNA_O2StarMover::apply( core::pose::Pose & pose, std::string & move_type )
	{

		using namespace core::pose::full_model_info;

		Size o2star_res( 0 );

		FullModelInfo & full_model_info = nonconst_full_model_info_from_pose( pose );
		utility::vector1< Size > moving_res_list = full_model_info.moving_res_list();

		if ( sample_all_o2star_ ){
			o2star_res = get_random_o2star_residue( pose );
		} else {
			// warning -- following might lead to weird 'hysteresis' effects since it picks
			// o2star to sample based on what's near moving residue.
			( *scorefxn_ )( pose ); //score it first to get energy graph.
			o2star_res = get_random_o2star_residue_near_moving_residue( pose, moving_res_list );
		}

		if ( o2star_res > 0 ) {
			Real const random_number2 = RG.uniform();
			if ( random_number2  < 0.5 ){
				move_type = "sml-o2star";
				sample_near_o2star_torsion( pose, o2star_res, sample_range_small_);
			} else {
				move_type = "lrg-o2star";
				sample_near_o2star_torsion( pose, o2star_res, sample_range_large_);
			}
		}

	}

	//////////////////////////////////////////////////////////////////////////////////////
	std::string
	RNA_O2StarMover::get_name() const {
		return "RNA_O2StarMover";
	}


	//////////////////////////////////////////////////////////////////////////////////////
	void
	RNA_O2StarMover::sample_near_o2star_torsion( pose::Pose & pose, Size const moving_res, Real const sample_range){
		id::TorsionID o2star_torsion_id( moving_res, id::CHI, 4 );
		Real o2star_torsion = pose.torsion( o2star_torsion_id ); //get
		o2star_torsion += RG.gaussian() * sample_range; //sample_near
		pose.set_torsion( o2star_torsion_id, o2star_torsion ); // apply
	}

	//////////////////////////////////////////////////////////////////////////////////////
	Size
	RNA_O2StarMover::get_random_o2star_residue( pose::Pose & pose ){
		// pick at random from whole pose -- a quick initial stab.
		Size const o2star_num = int( pose.total_residue() * RG.uniform() ) + 1;
		return o2star_num;
	}

	//////////////////////////////////////////////////////////////////////////////////////
	// This could be made smarter -- could go over nucleoside *and* suite.
	Size
	RNA_O2StarMover::get_random_o2star_residue_near_moving_residue( pose::Pose & pose, utility::vector1< Size > const moving_res_list ){

		// should be better -- actually look at o2star's that might be engaged in interactions with moving nucleoside

		// Following requires that pose is already scored once before entering the mover.
		//  Is there a way to avoid this? just ensure energy_graph is filled?
		utility::vector1< bool > residues_allowed_to_be_packed( pose.total_residue(), false );
		scoring::EnergyGraph const & energy_graph( pose.energies().energy_graph() );

		Distance DIST_CUTOFF( 4.0 );

		for (Size k = 1; k <= moving_res_list.size(); k++ ){
			Size const i = moving_res_list[ k ];

			for( graph::Graph::EdgeListConstIter
						 iter = energy_graph.get_node( i )->const_edge_list_begin();
					 iter != energy_graph.get_node( i )->const_edge_list_end();
					 ++iter ){

				Size j( (*iter)->get_other_ind( i ) );

				// check for potential interaction of o2* of this new residue and any atom in moving residue.
				Vector const & o2star_other = pose.residue( j ).xyz( " O2'" );
				for ( Size n = 1; n <= pose.residue( i ).natoms(); n++ ){
					if ( ( pose.residue( i ).xyz( n ) - o2star_other ).length() < DIST_CUTOFF ) {
						residues_allowed_to_be_packed[ j ] = true;
						break;
					}
				}

				// check for potential interaction of o2* of moving residue and any atom in this new residue
				if (residues_allowed_to_be_packed[ i ]) continue;

				Vector const & o2star_i = pose.residue( i ).xyz( " O2'" );
				for ( Size n = 1; n <= pose.residue( j ).natoms(); n++ ){
					if ( ( pose.residue( j ).xyz( n ) - o2star_i ).length() < DIST_CUTOFF ) {
						residues_allowed_to_be_packed[ i ] = true;
						break;
					}
				}

			}
		}

		utility::vector1< Size > res_list;
		for ( Size n = 1; n <= pose.total_residue(); n++ ) {
			if ( residues_allowed_to_be_packed[ n ] ) {
				res_list.push_back( n );
			}
		}
		if (res_list.size()==0) return 0; //nothing to move!

		Size const o2star_idx_within_res_list = int(  res_list.size() * RG.uniform() ) + 1;
		Size const o2star_num = res_list[ o2star_idx_within_res_list ];
		return o2star_num;
	}

}
}
}
