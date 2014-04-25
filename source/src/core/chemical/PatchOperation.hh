// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/PatchOperation.hh
/// @brief  Polymorphic classes representing the contents of a residue-type patch file
/// @author Phil Bradley

#ifndef INCLUDED_core_chemical_PatchOperation_hh
#define INCLUDED_core_chemical_PatchOperation_hh

// Unit headers
#include <core/chemical/PatchOperation.fwd.hh>

// Package headers
#include <core/chemical/ResidueType.fwd.hh>

// Basic header
#include <basic/Tracer.hh>

// Utility header
#include <utility/vector1.hh>


namespace core {
namespace chemical {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief  A single operation that needs to be applied in a residue patch
class PatchOperation : public utility::pointer::ReferenceCount {
public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~PatchOperation();

	/// @brief Returns true to signal failure, false to indicate success.
	virtual
	bool
	apply( ResidueType & rsd ) const = 0;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief delete an atom
class DeleteAtom : public PatchOperation {
public:

	/// @brief constructor
	DeleteAtom( std::string const & atom_name_in );

	/// @brief delete an atom from ResidueType rsd
	bool
	apply( ResidueType & rsd ) const;

private:
	/// name of the atom to be deleted
	std::string atom_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set an atom as backbone heavy atom
class SetBackboneHeavyatom : public PatchOperation {
public:
	/// @brief constructor
	SetBackboneHeavyatom( std::string const & atom_name_in );

	/// set an atom in ResidueType rsd as backbone heavy atom
	bool
	apply( ResidueType & rsd ) const;

private:
	// name of the atom to be set
	std::string atom_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set an atom as polymer connection
class SetPolymerConnectAtom : public PatchOperation {
public:

	/// @brief constructor the type of connection is "LOWER" or "UPPER"
	SetPolymerConnectAtom( std::string const & atom_name_in, std::string const & upper_lower_in );

	/// @brief set an atom in ResidueType rsd as a polymer connection atom
	bool
	apply( ResidueType & rsd ) const;

private:
	// "NONE" to delete the connection by setting its atomno to ZERO
	std::string atom_name_;
	// -1 for lower connection, 1 for upper connection
	int upper_lower_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class AddConnect : public PatchOperation {
public:

	// c-tor
	AddConnect( std::string const & connect_atom,
		Real const phi, Real const theta, Real const d,
		std::string const & parent_atom,
		std::string const & angle_atom,
		std::string const & torsion_atom
	);

	/// add a property
	bool
	apply( ResidueType & rsd ) const;

private:
	std::string const connect_atom_;
	Real const phi_;
	Real const theta_;
	Real const d_;
	std::string const  parent_atom_;
	std::string const   angle_atom_;
	std::string const torsion_atom_;

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a property to ResidueType
class AddProperty : public PatchOperation {
public:

	/// @brief constructor
	AddProperty( std::string const & property_in );

	/// @brief add a property
	bool
	apply( ResidueType & rsd ) const;

private:
	/// property to be added
	std::string property_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief delete a property from ResidueType
///    Added by Andy M. Chen in June 2009
///    This is needed for deleting properties, which occurs in certain PTM's (e.g. methylation)
class DeleteProperty : public PatchOperation {
public:

	/// @brief constructor
	DeleteProperty( std::string const & property_in );

	/// @brief delete a property
	bool
	apply( ResidueType & rsd ) const;

private:
	// property to be added
	std::string property_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief   Add a chi angle to ResidueType.
/// @author  Added by Andy M. Chen in June 2009
/// @details This is needed for PTMs, which often result in one or more extra chi angles.
class AddChi : public PatchOperation {
public:
	/// @brief Constructor for when the chi index is specified.
	AddChi(Size const & chino_in,
			std::string const & atom1_in,
			std::string const & atom2_in,
			std::string const & atom3_in,
			std::string const & atom4_in);

	/// @brief Constructor for when the chi index is not specified.
	AddChi(std::string const & atom1_in,
			std::string const & atom2_in,
			std::string const & atom3_in,
			std::string const & atom4_in);

	/// @brief Add a chi angle.
	bool apply(ResidueType & rsd) const;

private:
	bool no_index_;  // indicates that no chi index is provided, and the new chi will be added to the end of the list
	Size chino_;

	// atoms defining the added chi angle
	std::string atom1_;
	std::string atom2_;
	std::string atom3_;
	std::string atom4_;
};


class AddProtonChi : public PatchOperation {
public:

	/// constructor
	AddProtonChi(
		Size const & chino_in,
		utility::vector1<core::Real> const & samples,
		utility::vector1<core::Real> const & extrasamples
	);

	/// add a proton chi angle
	bool
	apply( ResidueType & rsd ) const;

private:
	/// atoms between which a chi angle is added
	Size chino_;
	utility::vector1<core::Real> samples_;
	utility::vector1<core::Real> extrasamples_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Redefine a chi angle
///    Added by Andy M. Chen in June 2009
///    This is needed for certain PTMs
class RedefineChi : public PatchOperation {
public:
	/// constructor
	RedefineChi(
		Size const & chino_in,
		std::string const & atom1_in,
		std::string const & atom2_in,
		std::string const & atom3_in,
		std::string const & atom4_in
	);

	/// redefine a chi angle
	bool
	apply( ResidueType & rsd ) const;

private:
	/// atoms between which a chi angle is added
	Size chino_;
	std::string atom1_;
	std::string atom2_;
	std::string atom3_;
	std::string atom4_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief   Add a rotamer sample to a chi angle of the ResidueType.
/// @author  Added by Andy M. Chen in June 2009
/// @details This is needed for PTMs.
class AddChiRotamer : public PatchOperation {
public:
	/// @brief Constructor for when the chi index is specified
	AddChiRotamer(Size const & chino_in, Angle const & mean_in, Angle const & sdev_in);

	/// @brief Constructor for when the chi index is not specified
	AddChiRotamer(Angle const & mean_in, Angle const & sdev_in);

	/// @brief Add a rotamer sample.
	bool apply(ResidueType & rsd) const;

private:
	bool no_index_;  // indicates that no chi index is provided, and the rotamer sample will be added to the last chi
	Size chino_;
	Real mean_;
	Real sdev_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add an atom to ResidueType
#ifdef WIN32
class AddAtomWIN32 : public PatchOperation {
#else
class AddAtom : public PatchOperation {
#endif
public:

	/// constructor
#ifdef WIN32
	AddAtomWIN32(
#else
	AddAtom(
#endif
		std::string const & atom_name_in,
		std::string const & atom_type_name_in,
		std::string const & mm_atom_type_name_in,
		Real const charge
	);

	/// add an atom
	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_name_;
	std::string atom_type_name_;
	std::string mm_atom_type_name_;
	Real const charge_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a bond to ResidueType
class AddBond : public PatchOperation {
public:

	/// constructor
	AddBond(
		std::string const & atom1_in,
		std::string const & atom2_in
	);

	/// add a bond
	bool
	apply( ResidueType & rsd ) const;

private:
	/// atoms between which a bond is added
	std::string atom1_;
	std::string atom2_;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set an atom's charge
class SetAtomicCharge : public PatchOperation {
public:

	/// constructor
	SetAtomicCharge(
		std::string const & atom_name_in,
		Real const charge_in
	);

	/// set an atom's charge
	bool
	apply( ResidueType & rsd ) const;

private:
	/// atom's name
	std::string atom_name_;
	/// atom's charge
	Real charge_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set atom's chemical type
class SetAtomType : public PatchOperation {
public:
	/// constructor
	SetAtomType(
		std::string const & atom_name_in,
		std::string const & atom_type_name_in
	);

	/// set atom's chemical type
	bool
	apply( ResidueType & rsd ) const;

private:
	/// atom's name
	std::string atom_name_;
	/// atom's type name
	std::string atom_type_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set atom's chemical type
class SetIO_String : public PatchOperation {
public:
	/// constructor
	SetIO_String(
		std::string const & name3,
		char const name1
	);

	/// set atom's chemical type
	bool
	apply( ResidueType & rsd ) const;

private:
	std::string const name3_;
	char const name1_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the interchangeability_group string for a ResidueType
class SetInterchangeabilityGroup_String : public PatchOperation {
public:
	/// constructor
	SetInterchangeabilityGroup_String(
		std::string const & intgrp
	);

	bool
	apply( ResidueType & rsd ) const;

private:
	std::string const intgrp_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Append a string to the existing interchangeability_group string for a ResidueType
class AppendInterchangeabilityGroup_String : public PatchOperation {
public:
	/// constructor
	AppendInterchangeabilityGroup_String(
		std::string const & intgrp_addendum
	);

	bool
	apply( ResidueType & rsd ) const;

private:
	std::string const intgrp_addendum_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set atom's MM chemical type
class SetMMAtomType : public PatchOperation {
public:
	/// constructor
	SetMMAtomType(
		std::string const & atom_name_in,
		std::string const & mm_atom_type_name_in
	);

	/// set atom's chemical type
	bool
	apply( ResidueType & rsd ) const;

private:
	/// atom's name
	std::string atom_name_;
	/// atom's type name
	std::string mm_atom_type_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set an atom's AtomICoord
class SetICoor : public PatchOperation {
public:
	/// constructor
	SetICoor(
		std::string const & atom_in,
		Real const phi_in,
		Real const theta_in,
		Real const d_in,
		std::string const & stub1_in,
		std::string const & stub2_in,
		std::string const & stub3_in
	);

	/// set an atom's AtomICoord
	bool
	apply( ResidueType & rsd ) const;

private:
	/// atom's name
	std::string atom_;
	/// AtomICoord data
	Real phi_;
	Real theta_;
	Real d_;
	std::string stub1_;
	std::string stub2_;
	std::string stub3_;
};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a mainchain atom before the first mainchain atom
class PrependMainchainAtom : public PatchOperation {
public:
	/// @brief constructor
	PrependMainchainAtom( std::string const & atom_name_in );

	/// @brief set an atom to be the first mainchain atom
	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a mainchain atom after the last mainchain atom
class AppendMainchainAtom : public PatchOperation {
public:
	/// @brief constructor
	AppendMainchainAtom( std::string const & atom_name_in );

	/// @brief set an atom to be the last mainchain atom
	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the residue neighbor atom
class SetNbrAtom : public PatchOperation {
public:
	/// @brief constructor
	SetNbrAtom( std::string const & atom_name_in );

	/// @brief set the residue neighbor atom
	bool
	apply( ResidueType & rsd ) const;

private:
	std::string atom_name_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the residue neighbor radius
class SetNbrRadius : public PatchOperation {
public:
	/// @brief constructor
	SetNbrRadius( Real const & radius );

	/// @brief set the residue neighbor atom
	bool
	apply( ResidueType & rsd ) const;

private:
	Real radius_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Set orient atom selection mode.
class SetOrientAtom : public PatchOperation {
public:
	SetOrientAtom(bool force_nbr_atom_orient);

	bool
	apply( ResidueType & rsd ) const;

private:
	bool force_nbr_atom_orient_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the path to a rotamer library for an NCAA that is not in dunbrack
class NCAARotLibPath : public PatchOperation {
public:
	/// @brief constructor
	NCAARotLibPath( std::string const & path_in );

	/// @brief set the NCAA rotamer library path in the residue type
	bool
	apply( ResidueType & rsd ) const;

private:
	std::string path_;
};

/// @brief  Virtual constructor, returns 0 if no match
PatchOperationOP
patch_operation_from_patch_file_line( std::string const & line );

} // chemical
} // core

#endif
