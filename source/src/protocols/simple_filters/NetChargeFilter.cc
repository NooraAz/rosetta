// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/NetChargeFilter.cc
/// @brief
/// @author Dave La (davela@u.washington.edu)
// Project Headers

#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/types.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <map>
#include <numeric/random/random.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_filters/ScoreTypeFilter.hh>
#include <protocols/simple_filters/NetChargeFilter.hh>
#include <protocols/simple_filters/NetChargeFilterCreator.hh>
#include <protocols/simple_moves/ddG.hh>
#include <string>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <sstream>


namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.NetChargeFilter" );

protocols::filters::FilterOP
NetChargeFilterCreator::create_filter() const { return new NetChargeFilter; }

std::string
NetChargeFilterCreator::keyname() const { return "NetCharge"; }

NetChargeFilter::~NetChargeFilter(){}

void
NetChargeFilter::parse_my_tag( utility::tag::TagPtr const tag, moves::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	chain_ = tag->getOption<core::Size>( "chain", 0 );
	net_charge_max_ = tag->getOption<int>( "max", 100 );
	net_charge_min_ = tag->getOption<int>( "min", -100 );

	TR<<"Net charge will be caculated for chain " << chain_ << " with maximum cutoff " << net_charge_max_ << " and minimum cutoff " << net_charge_min_ << "." << std::endl;
}

bool
NetChargeFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const net_charge( compute( pose ) );
	TR<<"Net Charge: "<<net_charge<<". " ;
	bool const status_max = (net_charge <= net_charge_max_) ? (true) : (false);
	bool const status_min = (net_charge >= net_charge_min_) ? (true) : (false);
	bool const status = (status_max && status_min) ? (true) : (false);
	if( status_max && status_min ) {
			TR << "passing." << std::endl;
	}
	else {
			TR << "failing." << std::endl;
  }
	return status;
}

void
NetChargeFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const net_charge( compute( pose ) );
	out<<"Net Charge: "<< net_charge <<'\n';
}

core::Real
NetChargeFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const net_charge( compute( pose ) );
	return( net_charge );
}

core::Real
NetChargeFilter::compute( core::pose::Pose const & pose ) const {
	core::pose::Pose copy_pose = pose;

	int net_charge = 0;

	for ( core::Size i=1; i <= pose.total_residue(); ++i ) {
	
		core::Size const chain = copy_pose.residue( i ).chain();

		// Skip if current residue is not part of the chain specified. 
		// Otherwise, the default chain=0 means consider all chains.
		if (chain_ != 0) {
				if (chain != chain_) continue;
		}

		std::string arg_res ("ARG");
		std::string lys_res ("LYS");
			
		std::string asp_res ("ASP");
		std::string glu_res ("GLU");

		std::stringstream out;
		std::string cur_res;

		out << pose.aa(i);
		cur_res = out.str();

		if (arg_res.compare(cur_res) == 0) {
				TR << "AA:  +1  " << cur_res << " " << i << std::endl;
				net_charge++;
		}
		else if (lys_res.compare(cur_res) == 0) {
				TR << "AA:  +1  " << cur_res << " " << i << std::endl;
				net_charge++;
		}
		else if (asp_res.compare(cur_res) == 0) {
				TR << "AA:  -1  " << cur_res << " " << i << std::endl;
				net_charge--;
		}
		else if (glu_res.compare(cur_res) == 0) {
				TR << "AA:  -1  " << cur_res << " " << i << std::endl;
				net_charge--;
		}

	}

	TR << "The net charge is: " << net_charge << std::endl;

	return( net_charge );
}


}
}
