// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/Filter.cc
/// @brief
/// @detailed
///	  Contains currently:
///
///
/// @author Florian Richter, Sarel Fleishman (sarelf@uw.edu)

// Unit Headers
#include <protocols/filters/Filter.hh>

// Project Headers
// AUTO-REMOVED #include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

// Boost Headers
#include <boost/foreach.hpp>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <core/pose/util.hh>

#define foreach BOOST_FOREACH


static basic::Tracer TR("protocols.filters.Filter");

namespace protocols {
namespace filters {

#ifdef USELUA
void lregister_Filter( lua_State * lstate ) {
	luabind::module(lstate, "protocols")
	[
		luabind::namespace_("filters")
		[
			luabind::class_<Filter>("Filter")
				.def("apply", ( void (Filter::*)( core::io::serialization::PipeMap & )) &Filter::apply)
				.def("score", &Filter::score)
		]
	];
}
#endif

using namespace core;
typedef std::pair< std::string const, FilterCOP > StringFilter_pair;
typedef utility::tag::TagPtr TagPtr;

bool
FilterCollection::apply( core::pose::Pose const & pose ) const
{
	foreach(protocols::filters::FilterCOP filter, filters_){
		if( ! filter->apply( pose ) ){
			return false;
		}
	}

	return true;
}

void
FilterCollection::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	foreach(protocols::filters::FilterCOP filter, filters_){
		filter->report( out, pose );
	}
}

Filter::Filter()
	: utility::pointer::ReferenceCount(),
		type_( "UNDEFINED TYPE" ),
		scorename_("defaultscorename")
{}

Filter::Filter( std::string const & type )
	: utility::pointer::ReferenceCount(),
		type_( type ),
		scorename_("defaultscorename")
{}

Filter::Filter( Filter const & init )
	:	utility::pointer::ReferenceCount(),
		type_( init.type_ ),
		scorename_("defaultscorename"),
		user_defined_name_( init.user_defined_name_ )
{}

Filter::~Filter() {}

void
Filter::parse_my_tag(
	TagPtr const,
	moves::DataMap &,
	Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const & )
{}

} // filters
} // protocols
