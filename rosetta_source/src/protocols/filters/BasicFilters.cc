// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/BasicFilters.cc
/// @brief
/// @detailed
///	  Contains currently:
///
///
/// @author Florian Richter, Sarel Fleishman (sarelf@uw.edu), Rocco Moretti (rmoretti@u.washington.edu)

// Unit Headers
#include <protocols/filters/BasicFilters.hh>
#include <protocols/filters/BasicFilterCreators.hh>

// AUTO-REMOVED #include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <protocols/moves/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

// Package Headers

// Project Headers
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

// ObjexxFCL Headers

// Utility headers
#include <numeric/random/random.hh>

// Boost Headers
#include <boost/foreach.hpp>

#include <protocols/jobdist/Jobs.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#define foreach BOOST_FOREACH

//// C++ headers
static basic::Tracer TR("protocols.filters.Filter");
static numeric::random::RandomGenerator RG( 140789 ); // <- Magic number, do not change it!!!

namespace protocols {
namespace filters {

using namespace core;
typedef std::pair< std::string const, FilterCOP > StringFilter_pair;
typedef utility::tag::TagPtr TagPtr;
typedef core::pose::Pose Pose;

////////////////////////////////////////////////////////////////////////////////////////////////////

StochasticFilter::StochasticFilter() : Filter( "Stochastic" ) {}
StochasticFilter::~StochasticFilter() {}

StochasticFilter::StochasticFilter( core::Real const confidence )
	: Filter( "Stochastic" ), confidence_( confidence )
{}

bool
StochasticFilter::apply( Pose const & ) const
{
	if( confidence_ >= 0.999 ) return true;

	core::Real const random_number( RG.uniform() );
	if( random_number <= confidence_ ) {
		TR<<"stochastic filter returning false"<<std::endl;
		return false;
	}
	TR<<"stochastic filter returning true"<<std::endl;
	return true;
}


FilterOP
StochasticFilter::clone() const
{
	return new StochasticFilter( *this );
}

FilterOP
StochasticFilter::fresh_instance() const
{
	return new StochasticFilter();
}

void
StochasticFilter::parse_my_tag(
	TagPtr const tag,
	moves::DataMap &,
	Filters_map const &,
	moves::Movers_map const &,
	Pose const & )
{
	confidence_ = tag->getOption< core::Real >( "confidence", 1.0 );
	TR<<"stochastic filter with confidence "<<confidence_<<std::endl;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// @brief Used to define a compound logical statement involving other filters with
// AND, OR and XOR
CompoundFilter::CompoundFilter() : Filter( "CompoundStatement" ) {}
CompoundFilter::~CompoundFilter() {}

CompoundFilter::CompoundFilter( CompoundStatement const & compound_statement ) :
	Filter( "CompoundStatement" ),
	compound_statement_( compound_statement )
{}

bool
CompoundFilter::apply( Pose const & pose ) const
{
	bool const value( compute( pose ) );

	TR<<"Compound logical statement is "<<value<<"."<<std::endl;
	return( value );
}

FilterOP
CompoundFilter::clone() const
{
	return new CompoundFilter( *this );
}

FilterOP
CompoundFilter::fresh_instance() const
{
	return new CompoundFilter();
}

void
CompoundFilter::report( std::ostream & out, Pose const & pose ) const
{
	if( compound_statement_.size() == 2 ){
	//special case for filters that are defined with a confidence value. In that case, we want to report the value of the filter regardless of the stochastic filter
		bool confidence( false );
		CompoundStatement::const_iterator non_stochastic_filter;
		for( CompoundStatement::const_iterator it=compound_statement_.begin(); it!=compound_statement_.end(); ++it ){
			if( it->first->get_type() == "Stochastic" ) confidence = true;
			else non_stochastic_filter = it;
		}
		if( confidence ) non_stochastic_filter->first->report( out, pose );
	}
	bool const value( compute( pose ) );

	out<<"Compound filter returns: "<<value<<'\n';
}

core::Real
CompoundFilter::report_sm( Pose const & pose ) const
{
	if( compound_statement_.size() == 2 ){
	//special case for filters that are defined with a confidence value. In that case, we want to report the value of the filter regardless of the stochastic filter
		bool confidence( false );
		CompoundStatement::const_iterator non_stochastic_filter;
		for( CompoundStatement::const_iterator it=compound_statement_.begin(); it!=compound_statement_.end(); ++it ){
			if( it->first->get_type() == "Stochastic" ) confidence = true;
			else non_stochastic_filter = it;
		}
		if( confidence ) return( non_stochastic_filter->first->report_sm( pose ) );
	}
	bool const value( compute( pose ) );
	return( value );
}

bool
CompoundFilter::compute( Pose const & pose ) const
{
	bool value( true );

	for( CompoundStatement::const_iterator it=compound_statement_.begin(); it!=compound_statement_.end(); ++it ) {
		if( it - compound_statement_.begin() == 0 ){ // ignore first logical operation
			value = it->first->apply( pose );
			continue;
		}
		switch( it->second  ) {
			case ( AND ) : value = value && it->first->apply( pose ); break;
			case ( OR  ) : value = value || it->first->apply( pose ); break;
			case ( XOR ) : value = value ^ it->first->apply( pose ); break;
			case ( NOR ) : value = value || !it->first->apply( pose ); break;
			case (NAND ) : value = value && !it->first->apply( pose ); break;
		}
	}
	return( value );
}

void
CompoundFilter::clear()
{
	compound_statement_.clear();
}

CompoundFilter::iterator
CompoundFilter::begin()
{
	return( compound_statement_.begin() );
}
CompoundFilter::const_iterator
CompoundFilter::begin() const
{
	return( compound_statement_.begin() );
}

CompoundFilter::iterator
CompoundFilter::end()
{
	return( compound_statement_.end() );
}

CompoundFilter::const_iterator
CompoundFilter::end() const
{
	return( compound_statement_.end() );
}

/// @details call the compound statement's constituent filters' set_resid
void
CompoundFilter::set_resid( core::Size const resid )
{
	for( iterator it( compound_statement_.begin() ); it!=compound_statement_.end(); ++it )
		protocols::moves::modify_ResId_based_object( it->first, resid );
}

void
CompoundFilter::parse_my_tag(
	TagPtr const tag,
	moves::DataMap &,
	Filters_map const & filters,
	moves::Movers_map const &,
	Pose const & )
{
	TR<<"CompoundStatement"<<std::endl;

	foreach(TagPtr cmp_tag_ptr, tag->getTags() ){
		std::string const operation( cmp_tag_ptr->getName() );
		std::pair< FilterOP, boolean_operations > filter_pair;
		if( operation == "AND" ) filter_pair.second = AND;
		else if( operation == "OR" ) filter_pair.second = OR;
		else if( operation == "XOR" ) filter_pair.second = XOR;
		else if( operation == "NOR" ) filter_pair.second = NOR;
		else if( operation =="NAND" ) filter_pair.second = NAND;
		else {
			utility_exit_with_message( "Error: Boolean operation in tag is undefined." );
		}
		std::string const filter_name( cmp_tag_ptr->getOption<std::string>( "filter_name" ) );

		Filters_map::const_iterator find_filter( filters.find( filter_name ));
		bool const filter_found( find_filter!=filters.end() );
		if( filter_found )
			filter_pair.first = find_filter->second->clone();
		else {
			TR<<"***WARNING WARNING! Filter defined for CompoundStatement not found in filter_list!!!! Defaulting to truefilter***"<<std::endl;
			filter_pair.first = new filters::TrueFilter;
		}
		runtime_assert( filter_found );
		compound_statement_.push_back( filter_pair );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// @brief Used to combine multiple seperate filters into a single filter value
	CombinedFilter::CombinedFilter() : Filter( "CombinedValue" ), threshold_(0.0) {}
CombinedFilter::~CombinedFilter() {}

bool
CombinedFilter::apply( core::pose::Pose const & pose ) const
{
	bool const value( compute( pose ) );

	TR<<"CombinedFilter value is "<<value<<", threshold is "<< threshold_ << "."<<std::endl;
	return( value <= threshold_ );
}

FilterOP
CombinedFilter::clone() const
{
	return new CombinedFilter( *this );
}

FilterOP
CombinedFilter::fresh_instance() const
{
	return new CombinedFilter();
}

void
CombinedFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	Real value( compute( pose ) );

	out<<"Combined filter returns: "<<value<<'\n';
}

core::Real
CombinedFilter::report_sm( core::pose::Pose const & pose ) const
{
	return compute(pose);
}

core::Real
CombinedFilter::compute( core::pose::Pose const & pose ) const
{
	core::Real value( 0.0 );

	foreach(FilterWeightPair fw_pair, filterlist_){
		value += fw_pair.second * fw_pair.first->report_sm( pose );
	}
	return( value );
}

void
CombinedFilter::parse_my_tag(
	TagPtr const tag,
	moves::DataMap &,
	Filters_map const & filters,
	moves::Movers_map const &,
	Pose const & )
{
	threshold_ = tag->getOption<core::Real>( "threshold", 0.0 );
	utility::vector1< TagPtr > const sub_tags( tag->getTags() );
	foreach(TagPtr tag_ptr, sub_tags){
		core::Real weight(1.0);
		if (tag_ptr->hasOption("factor") ) {
			weight = tag_ptr->getOption<core::Real>( "factor" );
		} else if (tag_ptr->hasOption("temp") ) {
			weight = 1.0 / tag_ptr->getOption<core::Real>( "temp" );
		}

		FilterOP filter;
		std::string const filter_name( tag_ptr->getOption<std::string>( "filter_name" ) );
		Filters_map::const_iterator find_filter( filters.find( filter_name ));
		bool const filter_found( find_filter!=filters.end() );
		if( filter_found ) {
			filter = find_filter->second->clone();
		}
		else {
			TR.Warning<<"***WARNING WARNING! Filter " << filter_name << " defined for CombinedValue not found in filter_list!!!! ***"<<std::endl;
			utility_exit_with_message("Filter "+filter_name+" not found in filter list.");
		}
		filterlist_.push_back( FilterWeightPair(filter, weight) );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Apply a sub-mover prior to calculating a filter value
MoveBeforeFilter::MoveBeforeFilter() : Filter( "MoveBeforeFilter" ) {}
MoveBeforeFilter::MoveBeforeFilter(moves::MoverOP mover, FilterCOP filter) : Filter( "MoveBeforeFilter" ), submover_(mover), subfilter_(filter) {}
MoveBeforeFilter::~MoveBeforeFilter() {}

bool
MoveBeforeFilter::apply( core::pose::Pose const & pose ) const
{
	core::pose::Pose p_mod(pose);
	submover_->apply(p_mod);
	return subfilter_->apply(p_mod);
}

FilterOP
MoveBeforeFilter::clone() const
{
	return new MoveBeforeFilter( *this );
}

FilterOP
MoveBeforeFilter::fresh_instance() const
{
	return new MoveBeforeFilter();
}

void
MoveBeforeFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	core::pose::Pose p_mod(pose);
	submover_->apply(p_mod);
	subfilter_->report(out, p_mod);
}

core::Real
MoveBeforeFilter::report_sm( core::pose::Pose const & pose ) const
{
	core::pose::Pose p_mod(pose);
	submover_->apply(p_mod);
	return subfilter_->report_sm(p_mod);
}

void
MoveBeforeFilter::parse_my_tag(
	TagPtr const tag,
	moves::DataMap &,
	Filters_map const & filters,
	moves::Movers_map const & movers,
	Pose const & )
{
	std::string mover_name("");
	std::string filter_name("");
	if ( tag->hasOption("mover") ) mover_name = tag->getOption< std::string >( "mover" );
	if ( tag->hasOption("mover_name") ) mover_name = tag->getOption< std::string >( "mover_name" );
	if ( tag->hasOption("filter") ) filter_name = tag->getOption< std::string >( "filter" );
	if ( tag->hasOption("filter_name") ) filter_name = tag->getOption< std::string >( "filter_name" );

	moves::Movers_map::const_iterator  find_mover ( movers.find( mover_name ));
  Filters_map::const_iterator find_filter( filters.find( filter_name ));

  if( find_mover == movers.end() ) {
    TR.Error << "ERROR !! mover '"<<mover_name<<"' not found in map: \n" << tag << std::endl;
    runtime_assert( find_mover != movers.end() );
  }
  if( find_filter == filters.end() ) {
    TR.Error << "ERROR !! filter '"<<filter_name<<"' not found in map: \n" << tag << std::endl;
    runtime_assert( find_filter != filters.end() );
  }
	submover_ = find_mover->second;
	subfilter_ = find_filter->second;

	TR << "Setting MoveBeforeFilter for mover '"<<mover_name<<"' and filter '"<<filter_name<<"'"<<std::endl;

}


////////////////////////////////////////////////////////////////////////////////////////////////////
// @brief FilterCreator methods

FilterOP
TrueFilterCreator::create_filter() const { return new TrueFilter; }

std::string
TrueFilterCreator::keyname() const { return "TrueFilter"; }

FilterOP
FalseFilterCreator::create_filter() const { return new FalseFilter; }

std::string
FalseFilterCreator::keyname() const { return "FalseFilter"; }

FilterOP
StochasticFilterCreator::create_filter() const { return new StochasticFilter; }

std::string
StochasticFilterCreator::keyname() const { return "Stochastic"; }

FilterOP
CompoundFilterCreator::create_filter() const { return new CompoundFilter; }

std::string
CompoundFilterCreator::keyname() const { return "CompoundStatement"; }

FilterOP
CombinedFilterCreator::create_filter() const { return new CombinedFilter; }

std::string
CombinedFilterCreator::keyname() const { return "CombinedValue"; }

FilterOP
MoveBeforeFilterCreator::create_filter() const { return new MoveBeforeFilter; }

std::string
MoveBeforeFilterCreator::keyname() const { return "MoveBeforeFilter"; }

} // filters
} // protocols
