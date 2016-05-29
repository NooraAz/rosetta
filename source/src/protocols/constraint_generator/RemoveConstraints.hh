// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/constraint_generator/RemoveConstraints.hh
/// @brief Removes constraints generated by a ConstraintGenerator
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_constraint_generator_RemoveConstraints_hh
#define INCLUDED_protocols_constraint_generator_RemoveConstraints_hh

// Unit headers
#include <protocols/constraint_generator/RemoveConstraints.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/constraint_generator/ConstraintGenerator.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/excn/EXCN_Base.hh>

namespace protocols {
namespace constraint_generator {

///@brief Removes constraints generated by a ConstraintGenerator
class RemoveConstraints : public protocols::moves::Mover {

public:
	RemoveConstraints();
	RemoveConstraints( ConstraintGeneratorCOPs const & generators );

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~RemoveConstraints();

	static std::string
	class_name() { return "RemoveConstraints"; }

	virtual void
	apply( core::pose::Pose & pose );

public:
	std::string
	get_name() const;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	/// @brief required in the context of the parser/scripting scheme
	virtual moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const;

public:
	void
	add_generator( protocols::constraint_generator::ConstraintGeneratorCOP generator );

private:
	ConstraintGeneratorCOPs generators_;

};

class EXCN_RemoveCstsFailed : public utility::excn::EXCN_Base {
public:
	EXCN_RemoveCstsFailed():
		utility::excn::EXCN_Base()
	{}
	virtual void show( std::ostream & os ) const { os << "Remodel constraints somehow got lost along the way" << std::endl; }
};

} //protocols
} //constraint_generator

#endif //protocols/constraint_generator_RemoveConstraints_hh