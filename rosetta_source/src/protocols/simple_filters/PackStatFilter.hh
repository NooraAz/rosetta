// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/simple_filters/PackStatFilter.hh
/// @brief header file for PackStatFilter class.
/// @detailed
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_simple_filters_PackStatFilter_hh
#define INCLUDED_protocols_simple_filters_PackStatFilter_hh

// Unit Headers
#include <protocols/simple_filters/PackStatFilter.fwd.hh>

// Package Headers
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>

// Utility headers

// Parser headers
#include <protocols/moves/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

#include <utility/vector1.hh>


// C++ headers

namespace protocols {
namespace simple_filters {

class PackStatFilter : public protocols::filters::Filter {
public:

	typedef protocols::filters::Filter Super;
	typedef protocols::filters::Filter Filter;
	typedef protocols::filters::FilterOP FilterOP;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;

	typedef utility::tag::TagPtr TagPtr;
	typedef protocols::filters::Filters_map Filters_map;
	typedef protocols::moves::DataMap DataMap;
	typedef protocols::moves::Movers_map Movers_map;


public:// constructor/destructor


	/// @brief default constructor
	PackStatFilter();

	/// @brief constructor with arguments
	PackStatFilter( Real const & ss );

	/// @brief copy constructor
	PackStatFilter( PackStatFilter const & rval );

	virtual ~PackStatFilter(){}


public:// virtual constructor


	/// @brief make clone
	virtual filters::FilterOP clone() const { return new PackStatFilter( *this ); }

	/// @brief make fresh instance
	virtual filters::FilterOP fresh_instance() const {	return new PackStatFilter(); }


public:// mutator


	/// @brief set filtered secondary structure
	void filtered_score( Real const & ss );


public:// accessor


	/// @brief get name of this filter
	virtual std::string name() const { return "PackStatFilter"; }



public:// parser

	virtual void parse_my_tag( TagPtr const tag,
														 DataMap &,
														 filters::Filters_map const &,
														 Movers_map const &,
														 Pose const & );


public:// virtual main operation


	/// @brief returns true if the given pose passes the filter, false otherwise.
	/// In this case, the test is whether the give pose is the topology we want.
	virtual bool apply( Pose const & pose ) const;

	/// @brief
	virtual Real report_sm( Pose const & pose ) const;

	/// @brief calc packstat score
	Real compute( Pose const & pose ) const;

	/// @brief determine which chain to evaluate, 0 is all
	core::Size chain_;
	
	/// @brief how many times to repeat evaluation
	core::Size repeats_;
	
private:

	Real filtered_score_;

};

} // filters
} // protocols

#endif
