// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/PatchOperation.hh
/// @brief  Polymorphic classes representing the contents of a residue-type patch file
/// @author Phil Bradley

#ifndef INCLUDED_core_chemical_PatchOperation_hh
#define INCLUDED_core_chemical_PatchOperation_hh

// Unit headers
#include <core/chemical/PatchOperation.fwd.hh>

// Package headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/MutableResidueType.fwd.hh>
#include <core/chemical/VariantType.hh>

// Utility header
#include <utility/vector1.hh>
#include <utility/VirtualBase.hh>

#include <map>

#include <core/chemical/ResidueProperty.hh> // AUTO IWYU For ResidueProperty, NO_PROPERTY

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief  A single operation that needs to be applied in a residue patch
class PatchOperation : public utility::VirtualBase {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from VirtualBase
	~PatchOperation() override;

	/// @brief Returns the name of the patch operation.  Useful for debugging.
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	virtual
	std::string name() const = 0;

	/// @brief Returns true to signal failure, false to indicate success.
	virtual
	bool
	apply( MutableResidueType & rsd ) const = 0;

	/// @brief   Which atom(s), if any, is/are added. Used for fast matching of ResidueType/Patches to PDB residues.
	/// @details This includes both atoms and all possible aliases for those atoms.
	virtual utility::vector1< std::string > adds_atoms() { return utility::vector1< std::string >(); }

	/// @brief Which atom, if any, is deleted. Used for fast matching of ResidueType/Patches to PDB residues.
	virtual
	std::string
	deletes_atom(){ return ""; }

	/// @brief Which property, if any, is added.
	virtual
	std::string
	adds_property() const { return ""; }

	/// @brief Which property, if any, is added.
	/// @details This returns an enum value.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	virtual
	ResidueProperty
	adds_property_enum() const { return NO_PROPERTY; }

	/// @brief Which property, if any, is deleted.
	virtual
	std::string
	deletes_property() const { return ""; }

	/// @brief Which property, if any, is deleted.
	/// @details This returns an enum value.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	virtual
	ResidueProperty
	deletes_property_enum() const { return NO_PROPERTY; }

	/// @brief Which variant, if any, is deleted.
	virtual
	std::string
	deletes_variant() const { return ""; }

	/// @brief Which variant, if any, is deleted, by enum.
	virtual
	VariantType
	deletes_variant_enum() const { return NO_VARIANT; }

	/// @brief Generates a new aa
	virtual
	bool
	may_change_aa(){ return false; }

	/// @brief Can this case change connections for the atom on the residue?
	/// @details - Be a little careful, as the passed atom name string may not have the same
	/// whitespace padding as any internal atom name.
	virtual
	bool
	changes_connections_on( ResidueType const &, std::string const & /*atom*/ ) const { return false; }

	/// @brief Generates name3.
	virtual
	std::string
	generates_name3(){ return ""; }

	/// @brief Generates interchangeability_group.
	virtual
	std::string
	generates_interchangeability_group(){ return ""; }

	/// @brief Generates base residue -- legacy for D_AA -- do not use otherwise.
	virtual
	bool
	generates_base_residue_type() const { return false; }

	/// @brief Special -- does this apply to 'minimal', placeholder types? Generally true, unless updating aa or name3.
	virtual
	bool
	applies_to_placeholder() const { return false; }

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief delete an atom
class DeleteAtom : public PatchOperation {
public:

	/// @brief constructor
	DeleteAtom( std::string const & atom_name_in );

	/// @brief delete an atom from ResidueType rsd
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("DeleteAtom").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

	std::string
	deletes_atom() override{ return atom_name_; }

private:
	/// name of the atom to be deleted
	std::string atom_name_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DeleteAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set an atom as backbone heavy atom
class SetBackboneHeavyatom : public PatchOperation {
public:
	/// @brief constructor
	SetBackboneHeavyatom( std::string const & atom_name_in );

	/// @brief Return the name of this PatchOperation ("SetBackboneHeavyatom").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

	/// set an atom in ResidueType rsd as backbone heavy atom
	bool
	apply( MutableResidueType & rsd ) const override;

private:
	// name of the atom to be set
	std::string atom_name_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetBackboneHeavyatom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set an atom as polymer connection
class SetPolymerConnectAtom : public PatchOperation {
public:

	/// @brief constructor the type of connection is "LOWER" or "UPPER"
	SetPolymerConnectAtom( std::string const & atom_name_in, std::string const & upper_lower_in );

	/// @brief set an atom in ResidueType rsd as a polymer connection atom
	bool
	apply( MutableResidueType & rsd ) const override;

	bool
	changes_connections_on( ResidueType const & rsd_type, std::string const & atom ) const override;

	/// @brief Return the name of this PatchOperation ("SetPolymerConnectAtom").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	// "NONE" to delete the connection by setting its atomno to ZERO
	std::string atom_name_;
	// -1 for lower connection, 1 for upper connection
	int upper_lower_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetPolymerConnectAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( MutableResidueType & rsd ) const override;

	bool
	changes_connections_on( ResidueType const & rsd_type, std::string const & atom ) const override;

	/// @brief Return the name of this PatchOperation ("AddConnect").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string connect_atom_;
	Real phi_;
	Real theta_;
	Real d_;
	std::string  parent_atom_;
	std::string   angle_atom_;
	std::string torsion_atom_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddConnect();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a property to ResidueType
class AddProperty : public PatchOperation {
public:

	/// @brief constructor
	AddProperty( std::string const & property_in );

	/// @brief add a property
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("AddProperty").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

	/// @brief Which property, if any, is added.
	std::string
	adds_property() const override { return property_; }

	/// @brief Which property, if any, is added.
	/// @details This returns an enum value.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	ResidueProperty
	adds_property_enum() const override { return property_enum_; }

private:
	/// @brief Name of property to be added.
	std::string property_;

	/// @brief Enum for property to be added.
	/// @details Will be NO_PROPERTY if this is a custom property type.
	ResidueProperty property_enum_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddProperty();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief delete a property from ResidueType
///    Added by Andy M. Chen in June 2009
///    This is needed for deleting properties, which occurs in certain PTM's (e.g. methylation)
///    Rewritten by Vikram K. Mulligan on 25 Aug. 2016 to use enums wherever possible for speed.
class DeleteProperty : public PatchOperation {
public:

	/// @brief constructor
	DeleteProperty( std::string const & property_in );

	/// @brief delete a property
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("DeleteProperty").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

	/// @brief Which property, if any, is deleted.
	///
	std::string
	deletes_property() const override { return property_name_; }

	/// @brief Which property, if any, is deleted.
	/// @details Returns NO_PROPERTY if this PatchOperation deletes a custom, on-the-fly property.
	ResidueProperty
	deletes_property_enum() const override { return property_; }


private:
	/// @brief Name (string) of property to be deleted.
	std::string property_name_;

	/// @brief Property to be deleted (enum).  Will be NO_PROPERTY if
	/// this deletes a custom (on-the-fly) property.
	ResidueProperty property_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DeleteProperty();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////
/// @brief   A patch operation for deleting a VariantType from a ResidueType.
/// @author  Labonte <JWLabonte@jhu.edu>
/// @author  Vikram K. Mulligan (vmullig@uw.edu) -- modified to primarily use enums instead of strings.
class DeleteVariantType : public PatchOperation {
public:
	// Constructor
	DeleteVariantType( std::string const & variant_in );

	// Constructor
	DeleteVariantType( VariantType const variant_in );

	/// @brief  Apply this patch to the given ResidueType.
	bool apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("DeleteVariantType").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

	/// @brief Which variant, if any, is deleted.
	std::string
	deletes_variant() const override { return variant_str_; }

	/// @brief Which variant, if any, is deleted, by enum.
	VariantType
	deletes_variant_enum() const override { return variant_; }

private:
	std::string variant_str_;
	VariantType variant_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DeleteVariantType();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	bool apply(MutableResidueType & rsd) const override;

	/// @brief Return the name of this PatchOperation ("AddChi").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	bool no_index_;  // indicates that no chi index is provided, and the new chi will be added to the end of the list
	Size chino_;

	// atoms defining the added chi angle
	std::string atom1_;
	std::string atom2_;
	std::string atom3_;
	std::string atom4_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddChi();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("AddProtonChi").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	/// atoms between which a chi angle is added
	Size chino_;
	utility::vector1<core::Real> samples_;
	utility::vector1<core::Real> extrasamples_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddProtonChi();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("RedefineChi").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	/// atoms between which a chi angle is added
	Size chino_;
	std::string atom1_;
	std::string atom2_;
	std::string atom3_;
	std::string atom4_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	RedefineChi();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Delete terminal chi angle
///    Added by Andrew M. Watkins in April 2015
class DeleteTerminalChi : public PatchOperation {
public:
	/// constructor
	DeleteTerminalChi() {};

	/// redefine a chi angle
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("DeleteTerminalChi").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Delete a metal binding atom
///    Added by Andrew M. Watkins in April 2015
class DeleteMetalbindingAtom : public PatchOperation {
public:
	/// constructor
	DeleteMetalbindingAtom( std::string const & atom_name );

	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("DeleteMetalbindingAtom").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom_name_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DeleteMetalbindingAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Delete an act coord atom
///    Added by Andrew M. Watkins in April 2015
class DeleteActCoordAtom : public PatchOperation {
public:
	/// constructor
	DeleteActCoordAtom( std::string const & atom_name );

	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("DeleteActCoordAtom").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom_name_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DeleteActCoordAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	bool apply(MutableResidueType & rsd) const override;

	/// @brief Return the name of this PatchOperation ("AddChiRotamer").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	bool no_index_;  // indicates that no chi index is provided, and the rotamer sample will be added to the last chi
	Size chino_;
	Real mean_;
	Real sdev_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddChiRotamer();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////
/// @brief   A patch operation for clearing all rotamer bins from the chi of a ResidueType.
/// @note    This is useful if one has redefined a chi.
/// @author  Labonte <JWLabonte@jhu.edu>
class ClearChiRotamers : public PatchOperation {
public:
	// Constructor
	ClearChiRotamers( core::uint const chi_no_in );

	/// @brief  Apply this patch to the given ResidueType.
	bool apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("ClearChiRotamers").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	core::uint chi_no_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ClearChiRotamers();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add an atom to ResidueType
class AddAtom : public PatchOperation {
public:

	/// constructor
	AddAtom(
		std::string const & atom_name_in,
		std::string const & atom_type_name_in,
		std::string const & mm_atom_type_name_in,
		Real const charge
	);

	/// @brief Return the name of this PatchOperation ("AddAtom").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

	/// add an atom
	bool
	apply( MutableResidueType & rsd ) const override;

	utility::vector1< std::string > adds_atoms() override { return utility::vector1< std::string >( 1, atom_name_ ); }

private:
	std::string atom_name_;
	std::string atom_type_name_;
	std::string mm_atom_type_name_;
	Real charge_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////
/// @brief   A patch operation for adding an atom alias to a ResidueType.
/// @note    See residue_io.cc for a description of atom aliases.
/// @remarks Atom aliases were graciously added to Rosetta by Rocco Moretti.
/// @author  Labonte <JWLabonte@jhu.edu>
class AddAtomAlias : public PatchOperation {
public:
	// Constructor
	AddAtomAlias(
		std::string const & rosetta_atom_name_in,
		utility::vector1< std::string > const & aliases_in,
		std::string const & preferred_alias_in );

	/// @brief  Apply this patch to the given ResidueType.
	bool apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("AddAtomAlias").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string name() const override;

	/// @brief   Return a list of all atom names that this operation adds (as aliases).
	/// @details All of the aliases for an atom must be offered as options to the ResidueTypeFinder.
	utility::vector1< std::string > adds_atoms() override { return aliases_; }

private:
	std::string rosetta_atom_name_;
	utility::vector1< std::string > aliases_;
	std::string preferred_alias_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddAtomAlias();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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

	/// @brief Return the name of this PatchOperation ("AddBond").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

	/// add a bond
	bool
	apply( MutableResidueType & rsd ) const override;

private:
	/// atoms between which a bond is added
	std::string atom1_;
	std::string atom2_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddBond();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////
/// @brief   A patch operation for adding a specific type of bond to a ResidueType.
/// @note    See residue_io.cc for a description of bond types.
/// @author  Labonte <JWLabonte@jhu.edu>
class AddBondType : public PatchOperation {
public:
	// Constructor
	AddBondType(
		std::string const & atom1_in,
		std::string const & atom2_in,
		std::string const & bond_type_in );

	/// @brief  Apply this patch to the given ResidueType.
	bool apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("AddBondType").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom1_;
	std::string atom2_;
	std::string bond_type_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddBondType();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////
/// @brief   A patch operation for changing the bond type of a given bond.
/// @note    See residue_io.cc for a description of bond types.
/// @author  Labonte <JWLabonte@jhu.edu>
class ChangeBondType : public PatchOperation {
public:
	// Constructor
	ChangeBondType(
		std::string const & atom1_in,
		std::string const & atom2_in,
		std::string const & old_bond_type_in,
		std::string const & new_bond_type_in );

	/// @brief  Apply this patch to the given ResidueType.
	bool apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("ChangeBondType").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom1_;
	std::string atom2_;
	std::string old_bond_type_;
	std::string new_bond_type_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ChangeBondType();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("SetAtomicCharge").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	/// atom's name
	std::string atom_name_;
	/// atom's charge
	Real charge_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetAtomicCharge();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////
/// @brief   A patch operation for setting the formal charge of a ResidueType's atom.
/// @author  Labonte <JWLabonte@jhu.edu>
class SetFormalCharge : public PatchOperation {
public:
	// Constructor
	SetFormalCharge( std::string const & atom_name_in, core::SSize charge_in );

	/// @brief  Apply this patch to the given ResidueType.
	bool apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("SetFormalCharge").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom_name_;
	core::SSize charge_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetFormalCharge();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////
/// @brief   A patch operation for setting the net formal charge of a whole ResidueType.
/// @author  Vikram K. Mulligan (vmullig@uw.edu)
class SetNetFormalCharge : public PatchOperation {
public:

	/// @brief Constructor
	/// @details Note that this deliberately takes a signed int.
	SetNetFormalCharge( signed int const charge_in );

	/// @brief  Apply this patch to the given ResidueType.
	bool apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("SetNetFormalCharge").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	signed int charge_;

#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetNetFormalCharge();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("SetAtomType").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	/// atom's name
	std::string atom_name_;
	/// atom's type name
	std::string atom_type_name_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetAtomType();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set residue's aa
class Set_AA : public PatchOperation {
public:
	/// constructor
	Set_AA(
		std::string const & aa
	) : aa_( aa ) {}

	/// set atom's chemical type
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("Set_AA").
	std::string
	name() const override { return "Set_AA"; }

	/// @brief Generates a new aa
	bool
	may_change_aa() override { return true; }

	bool
	applies_to_placeholder() const override { return true; }

private:
	std::string aa_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	Set_AA();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set residue's name1 and name3
class SetIO_String : public PatchOperation {
public:
	/// constructor
	SetIO_String(
		std::string const & name3,
		char const name1
	);

	/// set atom's chemical type
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("SetIO_String").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

	/// @brief Generates name3.
	std::string
	generates_name3() override{ return name3_; }


	bool
	applies_to_placeholder() const override { return true; }

private:
	std::string name3_;
	char name1_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetIO_String();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("SetInterchangeabilityGroup_String").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

	std::string
	generates_interchangeability_group() override{ return intgrp_; }

private:
	std::string intgrp_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetInterchangeabilityGroup_String();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("SetMMAtomType").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	/// atom's name
	std::string atom_name_;
	/// atom's type name
	std::string mm_atom_type_name_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetMMAtomType();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

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

	/// @brief Return the name of this PatchOperation ("SetICoor").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

	/// set an atom's AtomICoord
	bool
	apply( MutableResidueType & rsd ) const override;

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
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetICoor();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Change the parent, grandparent, or great-grandparent of an atom
enum Ancestor { anc_parent, anc_grandparent, anc_greatgrandparent };
class ChangeAncestory : public PatchOperation {
public:
	/// constructor
	ChangeAncestory(
		std::string const & target_atom,
		Ancestor which_ancestor,
		std::string const & ancestor_name
	);

	/// @brief change the ancestory, but leave the icoors intact.
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("ChangeAncestory").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	/// atom's name
	std::string atom_;
	Ancestor which_ancestor_;
	std::string ancestor_name_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ChangeAncestory();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////
/// @brief   A patch operation for resetting the length of a bond within a ResidueType.
/// @note    This is useful for when an atom is rehybridized within a patch file.
/// @author  Labonte <JWLabonte@jhu.edu>
class ResetBondLength : public PatchOperation {
public:
	// Constructor
	ResetBondLength( std::string const & atm_in, core::Distance d_in );

	/// @brief  Apply this patch to the given ResidueType.
	bool apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("ResetBondLength").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atm_;
	core::Distance d_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ResetBondLength();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a mainchain atom before the first mainchain atom
class PrependMainchainAtom : public PatchOperation {
public:
	/// @brief constructor
	PrependMainchainAtom( std::string const & atom_name_in );

	/// @brief set an atom to be the first mainchain atom
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("PrependMainchainAtom").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom_name_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	PrependMainchainAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a mainchain atom after the last mainchain atom
class AppendMainchainAtom : public PatchOperation {
public:
	/// @brief constructor
	AppendMainchainAtom( std::string const & atom_name_in );

	/// @brief set an atom to be the last mainchain atom
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("AppendMainchainAtom").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom_name_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AppendMainchainAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace a mainchain atom
class ReplaceMainchainAtom : public PatchOperation {
public:
	/// @brief constructor
	ReplaceMainchainAtom( std::string const & target, std::string const & new_atom );

	/// @brief set an atom to be the last mainchain atom
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("ReplaceMainchainAtom").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string target_;
	std::string new_atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceMainchainAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the residue neighbor atom
class SetNbrAtom : public PatchOperation {
public:
	/// @brief constructor
	SetNbrAtom( std::string const & atom_name_in );

	/// @brief set the residue neighbor atom
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("SetNbrAtom").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom_name_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetNbrAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the residue neighbor radius
class SetNbrRadius : public PatchOperation {
public:
	/// @brief constructor
	SetNbrRadius( Real const & radius );

	/// @brief set the residue neighbor atom
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("SetNbrRadius").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	Real radius_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetNbrRadius();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the residue neighbor radius
class SetAllAtomsRepulsive : public PatchOperation {
public:
	/// @brief constructor
	SetAllAtomsRepulsive() {};

	/// @brief set the residue neighbor atom
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("SetAllAtomsRepulsive").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Set orient atom selection mode.
class SetOrientAtom : public PatchOperation {
public:
	SetOrientAtom(bool force_nbr_atom_orient);

	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("SetOrientAtom").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	bool force_nbr_atom_orient_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetOrientAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Remove existing rotamer specifications (of any type).
/// @author Vikram K. Mulligan (vmullig@uw.edu)
class RemoveRotamerSpecifications : public PatchOperation {
public:
	/// @brief Constructor.
	///
	RemoveRotamerSpecifications();

	/// @brief Strip all RotamerSpecifications from the ResidueType.
	///
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("RemoveRotamerSpecifications").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Set the filenames for RamaPrePro scoring tables.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
class RamaPreproFilename : public PatchOperation {
public:

	/// @brief Constructor.
	///
	RamaPreproFilename( std::string const &non_prepro_file, std::string const &prepro_file );

	/// @brief Set the RamaPrepro library paths in the residue type.
	///
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("RamaPreproFilename").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:

	/// @brief The rama table file to use for scoring residues that do NOT precede proline.
	///
	std::string non_prepro_file_;

	/// @brief The rama table file to use for scoring residues that DO precede proline.
	///
	std::string prepro_file_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Set the residue name for RamaPrePro scoring tables.
/// @details This is the name in the scoring table AND the reference string used to look up the table.  Should be
/// unique.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
class RamaPreproResname : public PatchOperation {
public:

	/// @brief Constructor.
	///
	RamaPreproResname( std::string const &resname_in );

	/// @brief Set the RamaPrepro reference string in the residue type.
	///
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("RamaPreproResname").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:

	/// @brief The name listed for the amino acid type in the scoring table.
	/// @details Must be unique.
	std::string resname_;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set the path to a rotamer library for an NCAA that is not in dunbrack
class NCAARotLibPath : public PatchOperation {
public:
	/// @brief constructor
	NCAARotLibPath( std::string const & path_in );

	/// @brief set the NCAA rotamer library path in the residue type
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("NCAARotLibPath").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string path_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	NCAARotLibPath();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Set the mainchain torsion indices that a noncanonical rotamer library depends upon.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
class NCAARotLibBBTorsions : public PatchOperation {
public:
	/// @brief Constructor.
	NCAARotLibBBTorsions( utility::vector1< core::Size > const &torsions_in );

	/// @brief Set the mainchain torsion indices that a noncanonical rotamer library depends upon.
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("NCAARotLibBBTorsions").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	utility::vector1< core::Size > indices_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	NCAARotLibBBTorsions();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Set the number of rotamer bins per chi for an NCAA that is not in dunbrack.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
class NCAARotLibNumRotamerBins : public PatchOperation {
public:
	/// @brief Constructor
	///
	NCAARotLibNumRotamerBins( utility::vector1< core::Size > const &binsizes_in );

	/// @brief Set the number of rotamer bins per chi for an NCAA.
	///
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("NCAARotLibNumRotamerBins").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	utility::vector1 < core::Size > binsizes_;
};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Add a connection to the residue's sulfur and make a virtual proton to track the position of the connection atom
class ConnectSulfurAndMakeVirtualProton : public PatchOperation {
public:
	/// @brief constructor
	ConnectSulfurAndMakeVirtualProton() {};

	bool
	apply( MutableResidueType & rsd ) const override;

	bool
	changes_connections_on( ResidueType const & rsd_type, std::string const & atom ) const override;

	/// @brief Return the name of this PatchOperation ("ConnectSulfurAndMakeVirtualProton").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Alter the base name.
/// @details This creates a new base residue type, clearing anything after the colon in the name.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
class SetBaseName : public PatchOperation {
public:
	/// @brief Deafult constructor -- defaulted.
	SetBaseName() = default;

	/// @brief Constructor.
	SetBaseName( std::string const & new_base_name );

	/// @brief Generates a new aa
	bool
	may_change_aa() override{ return true; }

	/// @brief Apply this patch to generate a new base residue type.
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("SetBaseName").
	std::string
	name() const override;

	/// @brief This patch operaton DOES result in a new base residue type.
	bool
	generates_base_residue_type() const override;

private:

	std::string new_base_name_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Execute chiral flip (primarily: at CA)
class ChiralFlipNaming : public PatchOperation {
public:
	/// @brief constructor
	ChiralFlipNaming() {}// std::string const atom1, std::string const atom2 ): atom1_( atom1 ), atom2_( atom2 ) {};

	bool
	applies_to_placeholder() const override { return true; }

	/// @brief Generates a new aa
	bool
	may_change_aa() override{ return true; }

	/// @brief set the NCAA rotamer library path in the residue type
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("ChiralFlipNaming").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

	/// @brief This patch operaton DOES result in a new base residue type.
	bool
	generates_base_residue_type() const override;

private:

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Execute chiral flip (primarily: at CA)
class ChiralFlipAtoms : public PatchOperation {
public:
	/// @brief constructor
	ChiralFlipAtoms() {}// std::string const atom1, std::string const atom2 ): atom1_( atom1 ), atom2_( atom2 ) {};
	/// @brief set the NCAA rotamer library path in the residue type

	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("ChiralFlipAtoms").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with trifluoromethyl
class ReplaceProtonWithTrifluoromethyl : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithTrifluoromethyl( std::string const & atom ): atom_( atom ) {};

	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("ReplaceProtonWithTrifluoromethyl").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceProtonWithTrifluoromethyl();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with methyl
class ReplaceProtonWithMethyl : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithMethyl( std::string const & atom ): atom_( atom ) {};

	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("ReplaceProtonWithMethyl").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceProtonWithMethyl();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with methoxy
class ReplaceProtonWithMethoxy: public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithMethoxy( std::string const & atom ): atom_( atom ) {};

	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("ReplaceProtonWithMethoxy").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceProtonWithMethoxy();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with ethyl
class ReplaceProtonWithEthyl : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithEthyl( std::string const & atom ): atom_( atom ) {};

	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("ReplaceProtonWithEthyl").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceProtonWithEthyl();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with chlorine
class ReplaceProtonWithChlorine : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithChlorine( std::string const & atom ): atom_( atom ) {};

	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("ReplaceProtonWithChlorine").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceProtonWithChlorine();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with fluorine
class ReplaceProtonWithFluorine : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithFluorine( std::string const & atom ): atom_( atom ) {};

	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("ReplaceProtonWithFluorine").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceProtonWithFluorine();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with bromine
class ReplaceProtonWithBromine : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithBromine( std::string const & atom ): atom_( atom ) {};

	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("ReplaceProtonWithBromine").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceProtonWithBromine();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with iodine
class ReplaceProtonWithIodine : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithIodine( std::string const & atom ): atom_( atom ) {};

	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("ReplaceProtonWithIodine").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceProtonWithIodine();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief replace proton with hydroxyl
class ReplaceProtonWithHydroxyl : public PatchOperation {
public:
	/// @brief constructor
	ReplaceProtonWithHydroxyl( std::string const & atom ): atom_( atom ) {};

	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("ReplaceProtonWithHydroxyl").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	ReplaceProtonWithHydroxyl();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a connect and tracking virt to the atom
class AddConnectAndTrackingVirt : public PatchOperation {
public:
	/// @brief constructor
	AddConnectAndTrackingVirt( std::string const & atom ): atom_( atom ) {};

	bool
	apply( MutableResidueType & rsd ) const override;

	bool
	changes_connections_on( ResidueType const & rsd_type, std::string const & atom ) const override;

	/// @brief Return the name of this PatchOperation ("AddConnectAndTrackingVirt").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddConnectAndTrackingVirt();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief add a connect to the atom, delete child proton
class AddConnectDeleteChildProton : public PatchOperation {
public:
	/// @brief constructor
	AddConnectDeleteChildProton( std::string const & atom ): atom_( atom ) {};

	bool
	apply( MutableResidueType & rsd ) const override;

	bool
	changes_connections_on( ResidueType const & rsd_type, std::string const & atom ) const override;

	/// @brief Return the name of this PatchOperation ("AddConnectDeleteChildProton").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	AddConnectDeleteChildProton();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief delete child proton
class DeleteChildProton : public PatchOperation {
public:
	/// @brief constructor
	DeleteChildProton( std::string const & atom ): atom_( atom ) {};

	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("DeleteChildProton").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	std::string atom_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	DeleteChildProton();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief virtualize all
class VirtualizeAll: public PatchOperation {
public:
	/// @brief constructor
	VirtualizeAll() {};

	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("VirtualizeAll").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief virtualize sidechain
class VirtualizeSidechain: public PatchOperation {
public:
	/// @brief constructor
	VirtualizeSidechain() {};

	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("VirtualizeSidechain").
	/// @author Andrew M. Watkins (amw579@stanford.edu).
	std::string
	name() const override;

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set virtual shadow atoms
class SetVirtualShadow : public PatchOperation {
public:
	/// constructor
	SetVirtualShadow(
		std::string const & shadower_,
		std::string const & shadowee_
	);

	/// set atom's chemical type
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("SetVirtualShadow").
	/// @author Vikram K. Mulligan (vmullig@uw.edu).
	std::string
	name() const override;

private:
	/// atom's name
	std::string shadower_;
	/// atom's type name
	std::string shadowee_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetVirtualShadow();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief rename atom
class RenameAtom : public PatchOperation {
public:
	/// constructor
	RenameAtom(
		std::string const & old_name,
		std::string const & new_name
	): old_name_( old_name ), new_name_( new_name ) {}

	/// set atom's chemical type
	bool
	apply( MutableResidueType & rsd ) const override;

	/// @brief Return the name of this PatchOperation ("RenameAtom").
	std::string
	name() const override { return "RenameAtom"; }

private:
	std::string old_name_;
	std::string new_name_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	RenameAtom();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief set an atom as this residue's disulfide forming atom
class SetDisulfideAtomName : public PatchOperation {
public:
	/// @brief constructor
	SetDisulfideAtomName( std::string const & atom_name_in );

	/// @brief Return the name of this PatchOperation ("SetDisulfideAtomName").
	/// @author Andy Watkins (amw579@stanford.edu)
	std::string
	name() const override;

	/// set an atom in ResidueType rsd as backbone heavy atom
	bool
	apply( MutableResidueType & rsd ) const override;

private:
	// name of the atom to be set
	std::string atom_name_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;
	SetDisulfideAtomName();

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // chemical
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_chemical_PatchOperation )
#endif // SERIALIZATION


#endif
