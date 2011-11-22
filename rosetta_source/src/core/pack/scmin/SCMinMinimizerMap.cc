// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/scmin/SCMinMinimizerMap.cc
/// @brief  Class for identifying the sidechain DOFs in the AtomTree which are free during
///         any particular call to the minimizer.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Unit headers
#include <core/pack/scmin/SCMinMinimizerMap.hh>

// Package Headers
// AUTO-REMOVED #include <core/pack/rotamer_set/RotamerSet.fwd.hh>
// AUTO-REMOVED #include <core/pack/rotamer_set/RotamerSets.fwd.hh>
#include <core/pack/scmin/AtomTreeCollection.hh>

// Project Headers
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/optimization/DOF_Node.hh>

//#include <core/pose/Pose.fwd.hh>
//#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <core/scoring/DerivVectorPair.hh>
#include <utility/options/BooleanVectorOption.hh>


namespace core {
namespace pack {
namespace scmin {

SCMinMinimizerMap::SCMinMinimizerMap() :
	nactive_residues_( 0 ),
	dof_mask_( 1, false ),
	focused_residue_( 0 ),
	chi_start_for_focused_residue_( 0 ),
	nchi_added_for_focused_residue_( 0 ),
	n_active_dof_nodes_( 0 )
{}

SCMinMinimizerMap::~SCMinMinimizerMap() {}

void SCMinMinimizerMap::set_total_residue( Size total_residue )
{
	reset_dof_nodes();
	nactive_residues_ = 0;
	active_residues_.resize( total_residue );
	active_residue_index_for_res_.resize( total_residue );
	chi_start_for_active_residue_.resize( total_residue );
	active_residue_atom_to_dofnode_index_.resize( total_residue );
	atom_derivatives_.resize( total_residue );

	domain_map_.dimension( total_residue );
	atcs_for_residues_.resize( total_residue );
	if ( dof_nodes_.size() < total_residue ) dof_nodes_.reserve( total_residue ); // just a guess
	std::fill( active_residues_.begin(), active_residues_.end(), 0 );
	std::fill( active_residue_index_for_res_.begin(), active_residue_index_for_res_.end(), 0 );
	std::fill( chi_start_for_active_residue_.begin(), chi_start_for_active_residue_.end(), 0 );
	//std::fill( domain_map_.begin(), domain_map_.end(), 1 );
	for ( Size ii = 1; ii <= domain_map_.size(); ++ii ) domain_map_( ii ) = 1;
}

/// @brief Disable the minimization for all residues.  Ammortized O(1).
void SCMinMinimizerMap::clear_active_chi()
{
	for ( Size ii = 1; ii <= nactive_residues_; ++ii ) {
		domain_map_( active_residues_[ ii ] ) = 1;
		atcs_for_residues_[ active_residues_[ ii ]] = 0;
		active_residue_index_for_res_[ active_residues_[ ii ] ] = 0;

		active_residues_[ ii ] = 0;
		std::fill( active_residue_atom_to_dofnode_index_[ ii ].begin(), active_residue_atom_to_dofnode_index_[ ii ].end(), 0 );
		chi_start_for_active_residue_[ ii ] = 0;
	}
	nactive_residues_ = 0;
	reset_dof_nodes();
}

/// @details This should be called at most once per residue between calls to "clear_active_chi"
void SCMinMinimizerMap::activate_residue_chi( Size resindex )
{
	assert( domain_map_( resindex ) == 1 ); // activate_residue_chi should not have already been called.
	assert( active_residue_index_for_res_[ resindex ] == 0 ); // activate_residue_chi should not have already been called.

	domain_map_( resindex ) = 0;
	active_residues_[ ++nactive_residues_ ] = resindex;
	active_residue_index_for_res_[ resindex ] = nactive_residues_;

}

/// @brief Invoked during the depth-first traversal through the AtomTree.  The AtomTree
/// is indicating that a particular torsion is dependent on another torsion.  Record
/// that fact.
void
SCMinMinimizerMap::add_torsion(
	DOF_ID const & new_torsion,
	DOF_ID const & parent
)
{
	DOF_ID new_torsion_corrected( id::AtomID( new_torsion.atomno(), focused_residue_ ), new_torsion.type() );

	DOF_NodeOP pnode( 0 );
	if ( parent.valid() ) {
		Size chiind = atoms_representing_chis_[ parent.atomno() ];
		assert( chiind != 0 );
		Size const pind = chi_start_for_focused_residue_ + chiind - 1;
		assert( pind != 0 && pind <= n_active_dof_nodes_ );
		pnode = dof_nodes_[ pind ];
	}

	++nchi_added_for_focused_residue_;
	assert( atoms_representing_chis_[ new_torsion.atomno() ] == 0 );
	atoms_representing_chis_[ new_torsion.atomno() ] = nchi_added_for_focused_residue_;
	if ( n_active_dof_nodes_ >= dof_nodes_.size() ) {
		dof_nodes_.push_back( new optimization::DOF_Node( new_torsion_corrected, pnode ) );
		++n_active_dof_nodes_;
	} else {
		++n_active_dof_nodes_;
		dof_nodes_[ n_active_dof_nodes_ ]->set_id( new_torsion_corrected );
		dof_nodes_[ n_active_dof_nodes_ ]->set_parent( pnode );
	}
}

/// @brief Invoked during the depth-first traversal through the AtomTree; the atom
/// tree is indicating that a given atom is controlled by a particular DOF.  Record
/// that fact.
void
SCMinMinimizerMap::add_atom(
	AtomID const & atom_id,
	DOF_ID const & dof_id
)
{
	if ( dof_id.valid() ) {
		Size const chi = atoms_representing_chis_[ dof_id.atomno() ];
		assert( chi != 0 );
		Size const dofind = chi_start_for_focused_residue_ + chi - 1;
		assert( dofind <= n_active_dof_nodes_ );
		/// add the atom, but give it the correct residue id, since the input atom_id will be 1.
		dof_nodes_[ dofind ]->add_atom( id::AtomID( atom_id.atomno(), focused_residue_ ) );
	}
}

/// @brief Traverse the atom trees in preparation for minimization to tie together all the
/// DOFs and the atoms they control.
void
SCMinMinimizerMap::setup( AtomTreeCollectionOP trees )
{
	reset_dof_nodes();
	atom_tree_collection_ = trees;
	for ( Size ii = 1; ii <= nactive_residues_; ++ii ) {
		Size iiresid = active_residues_[ ii ];
		//assert( atcs_for_residues_[ iiresid ] == 0 );
		atcs_for_residues_[ iiresid ] = trees->residue_atomtree_collection_op( iiresid );
		active_residue_atom_to_dofnode_index_[ ii ].resize( atcs_for_residues_[ iiresid ]->active_residue().natoms(), 0 );
		focused_residue_ = iiresid;

		conformation::Residue const & iires( atcs_for_residues_[ iiresid ]->active_residue() );

		/// mark the DOFs in the DOF_ID_Mask for the chi in this residue as being free
		assert( dof_mask_.n_atom( 1 ) == atoms_representing_chis_.size() );
		if ( dof_mask_.n_atom( 1 ) < iires.natoms() ) {
			dof_mask_.resize( 1, iires.natoms(), false );
			atoms_representing_chis_.resize( iires.natoms(), 0 );
		}
		for ( Size jj = 1; jj <= iires.nchi(); ++jj ) {
			dof_mask_[ id::DOF_ID( id::AtomID( iires.chi_atoms( jj )[ 4 ], 1 ), id::PHI ) ] = true;
		}
		chi_start_for_focused_residue_ = n_active_dof_nodes_ + 1;
		nchi_added_for_focused_residue_ = 0;

		// descend through the atom tree for this residue.
		// In this call, add_atom() and add_torsion() wil be invoked on this MinimizerMap.
		DOF_ID invalid;
		atcs_for_residues_[ iiresid ]->active_atom_tree().root()->setup_min_map( invalid, dof_mask_, *this );

		// now mark the DOFs in the DOF_ID_Mask for the chi in this residue as being fixed
		for ( Size jj = 1; jj <= iires.nchi(); ++jj ) {
			dof_mask_[ id::DOF_ID( id::AtomID( iires.chi_atoms( jj )[ 4 ], 1 ), id::PHI ) ] = false;
		}
		std::fill( atoms_representing_chis_.begin(), atoms_representing_chis_.end(), 0 ); // overkill

	}
}

/// @brief Convenience lookup -- turns over the request to the AtomTreeCollection
conformation::Residue const &
SCMinMinimizerMap::residue( Size seqpos ) const
{
	return atcs_for_residues_[ seqpos ]->active_residue();
}

void
SCMinMinimizerMap::starting_dofs( optimization::Multivec & chi ) const
{
	chi.resize( n_active_dof_nodes_ );
	for ( Size ii = 1; ii <= n_active_dof_nodes_; ++ii ) {
		DOF_Node const & iinode( * dof_nodes_[ ii ] );
		Size const iirsd = iinode.rsd();
		Size const iichi = atcs_for_residues_[ iirsd ]->active_restype().last_controlling_chi( iinode.atomno() );
		assert( atcs_for_residues_[ iirsd ]->active_restype().chi_atoms( iichi )[ 4 ] == (Size) iinode.atomno() );
		chi[ ii ] = atcs_for_residues_[ iirsd ]->active_residue().chi( iichi );
	}
}

void
SCMinMinimizerMap::assign_dofs_to_mobile_residues( optimization::Multivec const & chi )
{
	assert( chi.size() == n_active_dof_nodes_ );
	for ( Size ii = 1; ii <= n_active_dof_nodes_; ++ii ) {
		DOF_Node const & iinode( * dof_nodes_[ ii ] );
		Size const iirsd = iinode.rsd();
		Size const iichi = atcs_for_residues_[ iirsd ]->active_restype().last_controlling_chi( iinode.atomno() );
		assert( atcs_for_residues_[ iirsd ]->active_restype().chi_atoms( iichi )[ 4 ] == (Size) iinode.atomno() );
		atcs_for_residues_[ iirsd ]->set_chi( iichi, chi[ ii ] );
	}
	for ( Size ii = 1; ii <= nactive_residues_; ++ii ) {
		atcs_for_residues_[ active_residues_[ ii ] ]->update_residue();
	}

}

SCMinMinimizerMap::DOF_Node const &
SCMinMinimizerMap::dof_node_for_chi( Size resid, Size chiid ) const
{
	Size atno = atcs_for_residues_[ resid ]->active_residue().chi_atoms( chiid )[ 4 ];
	Size actid = active_residue_index_for_res_[ resid ];
	Size dofnode_ind = active_residue_atom_to_dofnode_index_[ actid ][ atno ];
	return * dof_nodes_[ dofnode_ind ];
}

id::TorsionID
SCMinMinimizerMap::tor_for_dof( DOF_ID const & dofid ) const
{
	Size const rsd( dofid.rsd() );
	Size const chi( residue( rsd ).type().last_controlling_chi( dofid.atomno() ) );
	id::TorsionID torid( rsd, id::CHI, chi );
	return torid;
}

kinematics::tree::Atom const &
SCMinMinimizerMap::atom( AtomID const & atid ) const
{
	return atcs_for_residues_[ atid.rsd() ]->active_atom_tree().atom( id::AtomID( atid.atomno(), 1 ) );
}

void SCMinMinimizerMap::zero_atom_derivative_vectors()
{
	for ( Size ii = 1; ii <= nactive_residues_; ++ii ) {
		Size iiresid = active_residues_[ ii ];
		for ( Size jj = 1, jjend = atom_derivatives_[ iiresid ].size(); jj <= jjend; ++jj ) {
			atom_derivatives_[ iiresid ][ jj ].f1() = 0.0;
			atom_derivatives_[ iiresid ][ jj ].f2() = 0.0;
		}
	}
}

/// @details super simple -- no need for a sort (nor is there a need in the optimization::MinimizerMap for that matter).
void SCMinMinimizerMap::link_torsion_vectors()
{
	/// by construction, DOFs are added such that index( parent ) < index( child ).
	/// (i.e. we have a partial order of DOFs and DOF_Nodes are added in a monotonically increasing manner in that partial order)
	for ( Size ii = n_active_dof_nodes_; ii >= 1; --ii ) {
		dof_nodes_[ ii ]->link_vectors();
	}
}

void SCMinMinimizerMap::set_natoms_for_residue( Size resid, Size natoms )
{
	if ( atom_derivatives_[ resid ].size() < natoms ) {
		atom_derivatives_[ resid ].resize( natoms );
	}
}

void SCMinMinimizerMap::reset_dof_nodes()
{
	for ( Size ii = 1; ii <= n_active_dof_nodes_; ++ii ) {
		dof_nodes_[ ii ]->set_parent( 0 );
		dof_nodes_[ ii ]->clear_atoms();
	}
	n_active_dof_nodes_ = 0;
}



} // namespace scmin
} // namespace pack
} // namespace core

