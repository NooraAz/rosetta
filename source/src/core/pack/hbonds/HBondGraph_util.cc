// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pack/hbonds/HBondGraph_util.cc
/// @brief A collections of methods that are useful for dealing with HBondGraphs
/// @author Jack Maguire, jack@med.unc.edu

#include <core/pack/hbonds/HBondGraph_util.hh>
#include <core/pack/hbonds/MCHBNetInteractionGraph.hh>
#include <core/scoring/hbonds/graph/HBondGraph.hh>
#include <core/scoring/hbonds/graph/AtomLevelHBondGraph.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/interaction_graph/AminoAcidNeighborSparseMatrix.hh>
#include <core/conformation/Residue.hh>
#include <utility/graph/Graph.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <list>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSet_.hh>
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace hbonds {

using namespace core::scoring::hbonds::graph;

scoring::hbonds::graph::HBondGraphOP create_init_and_create_edges_for_hbond_graph(
	rotamer_set::RotamerSetsOP rotamer_sets,
	core::scoring::ScoreFunction const & sfxn,
	core::pose::Pose const & pose,
	core::Real hydrogen_bond_threshold,
	core::Real clash_threshold
){
	scoring::hbonds::graph::HBondGraphOP hbond_graph = create_and_init_hbond_graph( * rotamer_sets );

	MCHBNetInteractionGraphOP ig( new MCHBNetInteractionGraph ( hbond_graph, rotamer_sets, hydrogen_bond_threshold, clash_threshold ) );
	utility::graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, sfxn, rotamer_sets->task() );
	rotamer_sets->precompute_two_body_energies( pose, sfxn, packer_neighbor_graph, ig, true );

	//If you are running with symmetry, you still need to score one-body interactions

	return hbond_graph;
}

scoring::hbonds::graph::AtomLevelHBondGraphOP create_init_and_create_edges_for_atom_level_hbond_graph(
	rotamer_set::RotamerSetsOP rotamer_sets,
	core::scoring::ScoreFunction const & sfxn,
	core::pose::Pose const & pose,
	core::Real hydrogen_bond_threshold,
	core::Real clash_threshold
){
	scoring::hbonds::graph::AtomLevelHBondGraphOP hbond_graph = create_and_init_atom_level_hbond_graph( * rotamer_sets );

	MCHBNetInteractionGraphOP ig( new MCHBNetInteractionGraph ( hbond_graph, rotamer_sets, hydrogen_bond_threshold, clash_threshold ) );
	utility::graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, sfxn, rotamer_sets->task() );
	rotamer_sets->precompute_two_body_energies( pose, sfxn, packer_neighbor_graph, ig, true );

	//If you are running with symmetry, you still need to score one-body interactions

	determine_atom_level_edge_info_for_all_edges(
		* hbond_graph,
		* rotamer_sets,
		* core::scoring::hbonds::HBondDatabase::get_database(),
		pose.energies().tenA_neighbor_graph(),
		pose
	);

	utility::vector1< bool > include_all_resids( pose.size(), true );
	determine_atom_level_node_info_for_all_nodes(
		* hbond_graph,
		* rotamer_sets,
		include_all_resids,
		false
	);

	return hbond_graph;
}


void init_node_info( AbstractHBondGraph & graph, rotamer_set::RotamerSets const & rotamer_sets ){
	core::Size const nrot = rotamer_sets.nrotamers();

	for ( core::Size rot = 1; rot <= nrot; ++rot ) {
		HBondNode * node = graph.get_hbondnode( rot );
		core::Size const mres = rotamer_sets.moltenres_for_rotamer( rot );
		debug_assert( mres );
		node->set_moltenres( mres );
		node->set_local_rotamer_id( rotamer_sets.rotid_on_moltenresidue( rot ) );
	}//for rot
}


void find_hbonds_in_residue_pair(
	core::conformation::Residue const & resA,
	core::conformation::Residue const & resB,
	core::scoring::hbonds::HBondDatabase const & database,
	utility::graph::Graph const & tenA_neighbor_graph,
	core::scoring::hbonds::HBondSet & set
){
	unsigned int const residA = resA.seqpos();
	unsigned short int const num_nbrsA =
		tenA_neighbor_graph.get_node( residA )->num_neighbors_counting_self();
	unsigned short int const num_nbrsB =
		tenA_neighbor_graph.get_node( resB.seqpos() )->num_neighbors_counting_self();

	core::scoring::hbonds::identify_hbonds_1way(
		database,
		resA,
		resB,
		num_nbrsA,
		num_nbrsB,
		false, //bool const evaluate_derivative,
		false, //bool const exclude_don_bb,
		false, //bool const exclude_don_bsc,
		false, //bool const exclude_acc_scb,
		false, //bool const exclude_acc_sc,
		// output
		set //HBondSet & hbond_set,
		//Real ssdep_weight_factor = 1.0,
		//bool bond_near_wat = false
	);

	core::scoring::hbonds::identify_hbonds_1way(
		database,
		resB,
		resA,
		num_nbrsB,
		num_nbrsA,
		false, //bool const evaluate_derivative,
		false, //bool const exclude_don_bb,
		false, //bool const exclude_don_bsc,
		false, //bool const exclude_acc_scb,
		false, //bool const exclude_acc_sc,
		// output
		set //HBondSet & hbond_set,
		//Real ssdep_weight_factor = 1.0,
		//bool bond_near_wat = false
	);

}

void determine_atom_level_edge_info_for_all_edges(
	scoring::hbonds::graph::AtomLevelHBondGraph & hb_graph,
	rotamer_set::RotamerSets const & rotamer_sets,
	core::scoring::hbonds::HBondDatabase const & database,
	utility::graph::Graph const & tenA_neighbor_graph,
	core::pose::Pose const & pose,
	bool skip_edges_with_degree_zero,
	core::Real hbond_energy_threshold_for_satisfaction,
	core::conformation::symmetry::SymmetryInfoCOP symm_info
){
	for ( utility::graph::EdgeListIterator it = hb_graph.edge_list_begin();
			it != hb_graph.edge_list_end(); ++it ) {
		AtomLevelHBondEdge * edge = static_cast< AtomLevelHBondEdge * >( *it );

		if ( skip_edges_with_degree_zero ) {
			if (
					hb_graph.get_node( edge->get_first_node_ind()  )->num_edges() == 1 &&
					hb_graph.get_node( edge->get_second_node_ind() )->num_edges() == 1
					) {
				continue;
			}
		}

		determine_atom_level_edge_info(
			* edge,
			rotamer_sets,
			database,
			tenA_neighbor_graph,
			pose,
			hbond_energy_threshold_for_satisfaction,
			symm_info
		);
	}
}

Size
get_symm_ind_res(
	core::pose::Pose const & pose,
	Size const resid,
	core::conformation::symmetry::SymmetryInfoCOP symm_info,
	std::map< char, std::pair< core::Size, core::Size > > chain_bounds
) {
	Size resi_ind( resid );
	if ( resid > symm_info->num_independent_residues() ) {
		char resi_chain = pose.chain( resid );
		if ( symm_info->get_num_components() > 1 ) {
			std::map< char, std::pair< Size, Size > > const & component_bounds = symm_info->get_component_bounds();
			char resi_comp = symm_info->get_component_of_residue( resid );
			resi_ind = resid - chain_bounds[ resi_chain ].first + component_bounds.at( resi_comp ).first;
		} else {
			resi_ind = resid - chain_bounds[ resi_chain ].first + 1;
		}
	}
	return resi_ind;
}


void determine_atom_level_edge_info(
	AtomLevelHBondEdge & hb_edge,
	rotamer_set::RotamerSets const & rotamer_sets,
	core::scoring::hbonds::HBondDatabase const & database,
	utility::graph::Graph const & tenA_neighbor_graph,
	core::pose::Pose const & pose,
	core::Real hbond_energy_threshold_for_satisfaction,
	core::conformation::symmetry::SymmetryInfoCOP symm_info
){
	core::scoring::hbonds::HBondSet set;

	core::conformation::Residue const & resA = * rotamer_sets.rotamer( hb_edge.get_first_node_ind() );
	core::conformation::Residue const & resB = * rotamer_sets.rotamer( hb_edge.get_second_node_ind() );
	core::Size const residA = resA.seqpos();
	core::Size const residB = resB.seqpos();

	if ( symm_info ) {

		//Stolen from HBNet.cc
		std::map< char, std::pair< core::Size, core::Size > > chain_bounds;
		for ( Size ic = 1; ic <= pose.conformation().num_chains(); ++ic ) {
			Size ic_begin = pose.conformation().chain_begin( ic );
			Size ic_end = pose.conformation().chain_end( ic );
			char chain = pose.chain( ic_begin );
			chain_bounds[ chain ].first = ic_begin;
			chain_bounds[ chain ].second = ic_end;
		}

		std::list< core::Size > resids_for_A_clones;
		std::list< core::Size > resids_for_B_clones;
		for ( core::Size resid = 1; resid <= pose.size(); ++resid ) {
			if ( resid == residA || resid == residB ) continue;

			core::Size const resid_ind = get_symm_ind_res( pose, resid, symm_info, chain_bounds );
			if ( resid_ind == residA ) {
				resids_for_A_clones.push_back( resid );
			}
			if ( resid_ind == residB ) {
				resids_for_B_clones.push_back( resid );
			}
		}

		if ( residA != residB ) {
			find_hbonds_in_residue_pair( resA, resB, database, tenA_neighbor_graph, set );
		}

		//We do not need to do an all-to-all comparison. We just need to iterate over one set of symmetric resids. If A is not symmetric, iterate over B
		if ( resids_for_A_clones.size() ) {
			core::pack::rotamer_set::symmetry::SymmetricRotamerSet_COP sym_set = utility::pointer::dynamic_pointer_cast< core::pack::rotamer_set::symmetry::SymmetricRotamerSet_ const > ( rotamer_sets.rotamer_set_for_residue( residA ) );
			runtime_assert( sym_set );
			for ( core::Size sympos : resids_for_A_clones ) {
				core::conformation::ResidueOP resA_ii = sym_set->orient_rotamer_to_symmetric_partner( pose, resA, sympos );
				find_hbonds_in_residue_pair( * resA_ii, resB, database, tenA_neighbor_graph, set );
			}
		} else if ( resids_for_B_clones.size() ) {
			core::pack::rotamer_set::symmetry::SymmetricRotamerSet_COP sym_set = utility::pointer::dynamic_pointer_cast< core::pack::rotamer_set::symmetry::SymmetricRotamerSet_ const > ( rotamer_sets.rotamer_set_for_residue( residB ) );
			runtime_assert( sym_set );
			for ( core::Size sympos : resids_for_B_clones ) {
				core::conformation::ResidueOP resB_ii = sym_set-> orient_rotamer_to_symmetric_partner( pose, resB, sympos );
				find_hbonds_in_residue_pair( resA, * resB_ii, database, tenA_neighbor_graph, set );
			}
		}

	} else {
		find_hbonds_in_residue_pair( resA, resB, database, tenA_neighbor_graph, set );
	}

	unsigned short int const num_hbonds = set.nhbonds();
	for ( unsigned short int ii = 1; ii <= num_hbonds; ++ii ) {
		core::scoring::hbonds::HBond const & hbond = set.hbond( ii );
		if ( hbond.energy() > hbond_energy_threshold_for_satisfaction ) continue;

		bool const first_node_is_donor = hbond.don_res() == residA;
		unsigned short int const Hatm = hbond.don_hatm();
		unsigned short int const Datm = ( first_node_is_donor ? resA.atom_base( Hatm ) : resB.atom_base( Hatm ) );
		unsigned short int const Aatm = hbond.acc_atm();

		if ( false ) { //TODO bridging_waters
			//LKAtomLevelHBondEdge * lk_edge = static_cast< LKAtomLevelHBondEdge * >( lk_edge );
			//LKHBondInfo info = ...
			//lk_edge->register_hbond( info, first_node_is_donor, Hatm, Datm, Aatm );
		} else {
			hb_edge.register_hbond( first_node_is_donor, Hatm, Datm, Aatm );
		}
	}// for all hbonds
}

void determine_atom_level_node_info_for_all_nodes(
	scoring::hbonds::graph::AtomLevelHBondGraph & hb_graph,
	rotamer_set::RotamerSets const & rotamer_sets,
	utility::vector1< bool > const & include_these_resids,
	bool skip_nodes_with_no_edges
){
	core::Size const num_nodes = hb_graph.num_nodes();
	for ( unsigned int ii=1; ii <= num_nodes; ++ii ) {
		AtomLevelHBondNode * al_node = static_cast< AtomLevelHBondNode * >( hb_graph.get_node( ii ) );
		if ( skip_nodes_with_no_edges && al_node->num_edges() == 0 ) {
			continue;
		}

		determine_atom_level_node_info( * al_node, rotamer_sets, include_these_resids );
	}
}

void determine_atom_level_node_info(
	AtomLevelHBondNode & al_node,
	rotamer_set::RotamerSets const & rotamer_sets,
	utility::vector1< bool > const & include_these_resids
){
	core::Size const rot = al_node.global_rotamer_id();
	core::Size const mres = rotamer_sets.moltenres_for_rotamer( rot );

	core::conformation::ResidueCOP rotamer = rotamer_sets.rotamer( rot );
	debug_assert( rotamer->seqpos() == rotamer_sets.moltenres_2_resid( mres ) );
	if ( ! include_these_resids[ rotamer->seqpos() ] ) return;

	////////
	//HEAVY
	unsigned short int const nheavyatoms = rotamer->nheavyatoms();
	for ( unsigned short int ii = rotamer->first_sidechain_atom(); ii <= nheavyatoms; ++ii ) {

		bool const is_don = rotamer->heavyatom_is_an_acceptor( ii );
		bool const is_acc = rotamer->heavyatom_has_polar_hydrogens( ii );
		if ( is_don || is_acc ) {
			al_node.add_polar_atom(
				ii,
				rotamer->xyz( ii ),
				false, //is_hydrogen
				is_don,
				is_acc,
				rotamer->atom_type( ii ).name() == "OH" //is_hydroxyl
			);
		}
	}//for heavy atoms

	//////
	//HPOL
	unsigned short int const natoms = rotamer->natoms();
	for ( unsigned short int jj = nheavyatoms + 1; jj <= natoms; ++jj ) {
		if ( rotamer->atom_is_polar_hydrogen( jj ) ) {
			al_node.add_polar_atom(
				jj,
				rotamer->xyz( jj ),
				true, //is_hydrogen
				false,
				false,
				rotamer->atom_type( rotamer->atom_base( jj ) ).name() == "OH" //is_hydroxyl
			);
		}
	}

	if ( false ) { //TODO bridging waters
		//LKAtomLevelHBondNode * lkal_node = static_cast< LKAtomLevelHBondNode * >( al_node );
	}
}



void find_satisfying_interactions_with_background(
	AtomLevelHBondNode & node,
	rotamer_set::RotamerSets const & rotamer_sets,
	utility::graph::Graph const & packer_neighbor_graph,
	core::pose::Pose const & poly_ala_pose,
	core::Real hbond_energy_threshold_for_satisfaction
) {

	core::conformation::ResidueCOP const rotamer = rotamer_sets.rotamer( node.global_rotamer_id() );
	unsigned int const resid = rotamer->seqpos();
	debug_assert( resid == rotamer_sets.moltenres_2_resid( node.moltenres() ) );

	utility::graph::Node const * const resid_nbr_node = packer_neighbor_graph.get_node( resid );
	core::scoring::hbonds::HBondSet set;
	auto const & database = * core::scoring::hbonds::HBondDatabase::get_database( set.hbond_options().params_database_tag() );

	//////////////////
	// FIND ALL HBONDS
	for ( utility::graph::EdgeListConstIterator it = resid_nbr_node->const_edge_list_begin(),
			end = resid_nbr_node->const_edge_list_end();
			it != end;
			++it ) {

		unsigned int const other_resid = (*it)->get_other_ind( resid );
		core::conformation::Residue const & other_residue = poly_ala_pose.residue( other_resid );
		find_hbonds_in_residue_pair( *rotamer, other_residue, database, poly_ala_pose.energies().tenA_neighbor_graph(), set );

	}//for all neighbors


	///////////////////////////////////////
	// REMOVE ATOM INFO FOR SATISFIED ATOMS
	//
	// For every rotamer in the hbond graph,
	// determine which atoms are hydrogen bonding
	// with the background and remove them from
	// the container of unsatisfied atoms
	unsigned short int const num_hbonds = set.nhbonds();
	for ( unsigned short int ii = 1; ii <= num_hbonds; ++ii ) {
		core::scoring::hbonds::HBond const & hbond = set.hbond( ii );
		if ( hbond.energy() > hbond_energy_threshold_for_satisfaction ) continue;

		if ( hbond.don_res() == resid ) {
			unsigned short int const Hatm = hbond.don_hatm();
			unsigned short int const Datm = rotamer->atom_base( Hatm );
			node.remove_atom_info_stable( Hatm );
			node.remove_atom_info_stable( Datm );
		} else {
			debug_assert( hbond.acc_res() == resid );
			unsigned short int const Aatm = hbond.acc_atm();
			node.remove_atom_info_stable( Aatm );
		}//Acceptor
	}// for all hbonds
}

void delete_edges_with_degree_zero( scoring::hbonds::graph::AbstractHBondGraph & hb_graph ){
	for ( utility::graph::EdgeListIterator it = hb_graph.edge_list_begin(), end = hb_graph.edge_list_end();
			it != end; ) {

		utility::graph::Edge * edge = *it;

		if ( hb_graph.get_node( edge->get_first_node_ind() )->num_edges() == 1 &&
				hb_graph.get_node( edge->get_second_node_ind() )->num_edges() == 1
				) {
			++it;//very important to increment before deleting the edge that "it" is pointing to
			hb_graph.delete_edge( edge );
		} else {
			++it;
		}
	}
}


} //hbonds
} //pack
} //core