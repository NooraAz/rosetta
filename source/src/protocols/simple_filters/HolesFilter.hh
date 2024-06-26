// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/simple_filters/HolesFilter.hh
/// @brief header file for HolesFilter class.
/// @details
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_simple_filters_HolesFilter_hh
#define INCLUDED_protocols_simple_filters_HolesFilter_hh

// Unit Headers
#include <protocols/simple_filters/HolesFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Utility headers

// Parser headers
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>



//// C++ headers

namespace protocols {
namespace simple_filters {

class HolesFilter : public protocols::filters::Filter {
public:

	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef std::string String;

	typedef utility::tag::TagCOP TagCOP;
	typedef basic::datacache::DataMap DataMap;


public:// constructor/destructor


	// @brief default constructor
	HolesFilter();

	// @brief constructor with arguments
	//HolesFilter( Real const & ss );

	// @brief copy constructor
	HolesFilter( HolesFilter const & rval );

	~HolesFilter() override= default;


public:// virtual constructor


	// @brief make clone
	filters::FilterOP clone() const override { return utility::pointer::make_shared< HolesFilter >( *this ); }

	// @brief make fresh instance
	filters::FilterOP fresh_instance() const override { return utility::pointer::make_shared< HolesFilter >(); }


public:// mutator


	// @brief
	void filtered_value( Real const & ss );


public:// accessor


	// @brief get name of this filter


public:// parser

	void parse_my_tag( TagCOP tag,
		basic::datacache::DataMap & data
	) override;


public:// virtual main operation


	// @brief returns true if the given pose passes the filter, false otherwise.
	// In this case, the test is whether the give pose is the topology we want.
	bool apply( Pose const & pose ) const override;

	/// @brief
	Real report_sm( Pose const & pose ) const override;

	/// @brief used to report score
	void report( std::ostream & out, Pose const & pose ) const override;

	/// @brief
	Real compute( Pose const & pose ) const;

	std::string
	name() const override;

	static
	std::string
	class_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void
	set_residue_selector( core::select::residue_selector::ResidueSelectorOP residue_selector ) {
		residue_selector_ = residue_selector;
	}

	void
	set_threshold( Real threshold ) {
		filtered_value_ = threshold;
	}

	void
	set_exclude_bb_atoms( bool exclude_bb_atoms ) {
		exclude_bb_atoms_ = exclude_bb_atoms;
	}

	void
	set_normalize_per_atom( bool normalize_per_atom ) {
		normalize_per_atom_ = normalize_per_atom;
	}

	void
	set_normalize_per_residue( bool normalize_per_residue ) {
		normalize_per_residue_ = normalize_per_residue;
	}


private:


	Real filtered_value_;

	String cmd_;

	// options for resisdue selector case
	bool exclude_bb_atoms_;
	bool normalize_per_atom_;
	bool normalize_per_residue_;

	core::select::residue_selector::ResidueSelectorOP residue_selector_;
};

} // filters
} // protocols

#endif
