// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	devel/buns/BuriedUnsatPolarsFinder2.cc
/// @brief
/// @details
/// @author	Kevin Houlihan (khouli@unc.edu)
/// @author	Bryan Der

// Unit headers

#include <devel/buns/BuriedUnsatisfiedPolarsCalculator2.hh>

//Project Headers
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <devel/vardist_solaccess/VarSolDRotamerDots.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <core/conformation/Residue.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <numeric/conversions.hh>
#include <basic/options/keys/bunsat_calc2.OptionKeys.gen.hh>


// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
#include <basic/MetricValue.hh>

#include <cassert>

#include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>

static basic::Tracer TR("devel.buns.BuriedUnsatisfiedPolarsCalculator2");


namespace devel {
namespace buns {

using namespace core;
using namespace core::chemical;
using namespace core::conformation;
using namespace core::pose;
using namespace basic::options;
using namespace core::pose::metrics;
using namespace protocols::toolbox::pose_metric_calculators;
using namespace utility;
typedef numeric::xyzVector< core::Real > Vector;

BuriedUnsatisfiedPolarsCalculator2::BuriedUnsatisfiedPolarsCalculator2(
	std::string weak_bunsat_calc
) :
	all_bur_unsat_polars_( 0 ),
	special_region_bur_unsat_polars_(0),
	name_of_weak_bunsat_calc_( weak_bunsat_calc ),
	layered_sasa_(option[OptionKeys::bunsat_calc2::layered_sasa]),
	generous_hbonds_(option[OptionKeys::bunsat_calc2::generous_hbonds]),
	sasa_burial_cutoff_(option[OptionKeys::bunsat_calc2::sasa_burial_cutoff]),
	AHD_cutoff_(option[OptionKeys::bunsat_calc2::AHD_cutoff]),
	dist_cutoff_(option[OptionKeys::bunsat_calc2::dist_cutoff]),
	hxl_dist_cutoff_(option[OptionKeys::bunsat_calc2::hxl_dist_cutoff]),
	sulph_dist_cutoff_(option[OptionKeys::bunsat_calc2::sulph_dist_cutoff]),
	metal_dist_cutoff_(option[OptionKeys::bunsat_calc2::metal_dist_cutoff])
{
	atom_bur_unsat_.clear();
	residue_bur_unsat_polars_.clear();
	special_region_.clear();
	assert_calculators();
}
			
			
BuriedUnsatisfiedPolarsCalculator2::BuriedUnsatisfiedPolarsCalculator2(
	std::string weak_bunsat_calc,
	std::set< core::Size > const & special_region
) :
	all_bur_unsat_polars_(0),
	special_region_bur_unsat_polars_(0),
	name_of_weak_bunsat_calc_( weak_bunsat_calc ),
	special_region_( special_region ),
	layered_sasa_(option[OptionKeys::bunsat_calc2::layered_sasa]),
	generous_hbonds_(option[OptionKeys::bunsat_calc2::generous_hbonds]),
	sasa_burial_cutoff_(option[OptionKeys::bunsat_calc2::sasa_burial_cutoff]),
	AHD_cutoff_(option[OptionKeys::bunsat_calc2::AHD_cutoff]),
	dist_cutoff_(option[OptionKeys::bunsat_calc2::dist_cutoff]),
	hxl_dist_cutoff_(option[OptionKeys::bunsat_calc2::hxl_dist_cutoff]),
	sulph_dist_cutoff_(option[OptionKeys::bunsat_calc2::sulph_dist_cutoff]),
	metal_dist_cutoff_(option[OptionKeys::bunsat_calc2::metal_dist_cutoff])
{
	atom_bur_unsat_.clear();
	residue_bur_unsat_polars_.clear();
	assert_calculators();
}
			
			
void
BuriedUnsatisfiedPolarsCalculator2::assert_calculators() {
	if( !CalculatorFactory::Instance().check_calculator_exists( name_of_weak_bunsat_calc_ ) ) {
		if( name_of_weak_bunsat_calc_ != "default" ) TR <<
				"Attention: couldn't find the specified buried unsat calculator ( " <<
				name_of_weak_bunsat_calc_ << " ), instantiating default one." << std::endl;
		name_of_weak_bunsat_calc_ = "bunsat_calc2_default_weak_bunsat_calc";
		
		std::string sasa_calc_name("sasa_calc_name");
		if( !CalculatorFactory::Instance().check_calculator_exists( sasa_calc_name ) ) {
			if(layered_sasa_) {
				TR << "Registering VarSolDist SASA Calculator" << std::endl;
				CalculatorFactory::Instance().register_calculator(sasa_calc_name,
						new devel::vardist_solaccess::VarSolDistSasaCalculator());
			}
			else {
				TR << "Registering SASA Calculator" << std::endl;
				CalculatorFactory::Instance().register_calculator(sasa_calc_name,
						new pose::metrics::simple_calculators::SasaCalculatorLegacy());
			}
		}
		std::string num_hbonds_calc_name("num_hbonds_calc_name");
		if (!CalculatorFactory::Instance().check_calculator_exists( num_hbonds_calc_name ) ) {
			CalculatorFactory::Instance().register_calculator(num_hbonds_calc_name,
					new NumberHBondsCalculator());
		}
		if( !CalculatorFactory::Instance().check_calculator_exists( name_of_weak_bunsat_calc_ ) ){
			TR << "Registering new basic buried unsat calculator." << std::endl;
			CalculatorFactory::Instance().register_calculator( name_of_weak_bunsat_calc_,
					new protocols::toolbox::pose_metric_calculators::
					BuriedUnsatisfiedPolarsCalculator(sasa_calc_name, num_hbonds_calc_name,
					sasa_burial_cutoff_) );
		}
	}
}
			

void
BuriedUnsatisfiedPolarsCalculator2::lookup(
	std::string const & key,
	basic::MetricValueBase * valptr
) const
{
	
	if ( key == "all_bur_unsat_polars" ) {
		basic::check_cast( valptr, &all_bur_unsat_polars_,
				"all_bur_unsat_polars expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( all_bur_unsat_polars_ );

	} else if ( key == "special_region_bur_unsat_polars" ) {
		basic::check_cast( valptr, &special_region_bur_unsat_polars_,
				"special_region_bur_unsat_polars expects to return a Size" );
		(static_cast<basic::MetricValue<Size> *>(valptr))->set( special_region_bur_unsat_polars_ );
		
	} else if ( key == "atom_bur_unsat" ) {
		basic::check_cast( valptr, &atom_bur_unsat_,
				"atom_bur_unsat expects to return a id::AtomID_Map< bool >" );
		(static_cast<basic::MetricValue<id::AtomID_Map< bool > > *>(valptr))->set( atom_bur_unsat_ );
		
	} else if ( key == "residue_bur_unsat_polars" ) {
		basic::check_cast( valptr, &residue_bur_unsat_polars_,
				"residue_bur_unsat_polars expects to return a utility::vector1< Size >" );
		(static_cast<basic::MetricValue<utility::vector1< Size > > *>(valptr))->set( residue_bur_unsat_polars_ );
	
	} else {
		basic::Error() << "BuriedUnsatisfiedPolarsCalculator2 cannot compute the requested metric " << key << std::endl;
		utility_exit();
	}
	
} //lookup
			
			
			
std::string
BuriedUnsatisfiedPolarsCalculator2::print( std::string const & key ) const
{
	if ( key == "all_bur_unsat_polars" ) {
		return utility::to_string( all_bur_unsat_polars_ );
	} else if ( key == "special_region_bur_unsat_polars" ) {
		return utility::to_string( special_region_bur_unsat_polars_ );
	} else if ( key == "residue_bur_unsat_polars" ) {
		return utility::to_string( residue_bur_unsat_polars_ );
	}
	
	basic::Error() << "BuriedUnsatisfiedPolarsCalculator2 cannot compute metric " << key << std::endl;
	utility_exit();
	return "";
	
} //print
			
			
void
BuriedUnsatisfiedPolarsCalculator2::recompute( Pose const & pose )
{
	TR << "Recomputing buried unsats" << std::endl;
	this->show();

	all_bur_unsat_polars_ = 0;
	special_region_bur_unsat_polars_ = 0;
	
	if( pose.total_residue() != residue_bur_unsat_polars_.size() ) {
		residue_bur_unsat_polars_.resize( pose.total_residue() );
		atom_bur_unsat_.resize( pose.total_residue() );
	}
	
	if(generous_hbonds_)
		generous_hbond();
	
	TR << "Running basic buried unsat calc" << std::endl;
	basic::MetricValue< id::AtomID_Map< bool > > bunsat_atomid_map;
	pose.metric(name_of_weak_bunsat_calc_, "atom_bur_unsat", bunsat_atomid_map);
	
	//id::AtomID_Map< bool > bunsat_thorough_atomid_map(bunsat_atomid_map.value());
	//bunsats_thorough_check(pose, bunsat_thorough_atomid_map);
	atom_bur_unsat_ = bunsat_atomid_map.value();
	
	TR << "Validating buried unsats" << std::endl;
	
	bunsats_thorough_check(pose, atom_bur_unsat_);
	
} //recompute

void
BuriedUnsatisfiedPolarsCalculator2::generous_hbond()
const
{
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
	scoring::methods::EnergyMethodOptionsOP emopts(
			new scoring::methods::EnergyMethodOptions(
			scorefxn->energy_method_options() ) );
	emopts->hbond_options().decompose_bb_hb_into_pair_energies(true);
	emopts->hbond_options().use_hb_env_dep(false);
	emopts->hbond_options().bb_donor_acceptor_check(false);
	scorefxn->set_energy_method_options(*emopts);
}

void
BuriedUnsatisfiedPolarsCalculator2::bunsats_thorough_check(
	pose::Pose const & pose,
	id::AtomID_Map< bool > & bunsat_thorough_atomid_map
)
{
	for (Size res = 1; res <= bunsat_thorough_atomid_map.n_residue(); ++res) {
		
		residue_bur_unsat_polars_[res] = 0;
		
		for (Size atm = 1; atm <= bunsat_thorough_atomid_map.n_atom(res); ++atm) {

			std::string res_debug = pose.pdb_info()->pose2pdb(res);
			std::string atm_debug = pose.residue(res).atom_name(static_cast< int >(atm));
			
			// if unsat by weak bunsat calc, do extra checks
			if (bunsat_thorough_atomid_map(res, atm)) {
				//TR << "Considering candidate buried unsat, residue " << res_debug << ", atom " << atm_debug << "." << std::endl;
				id::AtomID bunsat_candidate_atom_id(atm, res);
				bunsat_thorough_atomid_map.set(bunsat_candidate_atom_id,
						single_bunsat_thorough_check(pose, bunsat_candidate_atom_id));
				//if (!bunsat_thorough_atomid_map(res, atm)) TR << "Rejected a buried unsat, residue " << res_debug << ", atom " << atm_debug << "." << std::endl;
			}
			// if still bunsat, add to counts
			if (bunsat_thorough_atomid_map(res, atm)) {
				TR << "Validated a buried unsat, rosetta resi " << res << ", pdb residue " << res_debug << ", atom " << atm_debug << "." << std::endl;
				residue_bur_unsat_polars_[res]++;
				if( special_region_.find( res ) != special_region_.end() )
						special_region_bur_unsat_polars_++;
			}
		}
	}
}

/*
 return true if bunsat_candidate_atom_id in pose passes additional checks,
 otherwise return false
 */
bool
BuriedUnsatisfiedPolarsCalculator2::single_bunsat_thorough_check(
	pose::Pose const & pose,
	id::AtomID const & bunsat_candidate_atom_id
)
{
	Size const bunsat_resi = bunsat_candidate_atom_id.rsd();
	
	//TR << "Calculating neighbors for resi " << bunsat_resi << std::endl;
	
	NeighborsByDistanceCalculator nbr_calc(bunsat_resi);
	basic::MetricValue< std::set< Size > > neighbors;
	nbr_calc.get("neighbors", neighbors, pose);
	//pose.metric(nbr_calc_name, "neighbors", neighbors);

	Residue const & bunsat_rsd(pose.residue(bunsat_resi));
	Size const bunsat_atom_num = bunsat_candidate_atom_id.atomno();
	
	Vector bunsat_xyz = bunsat_rsd.atom(bunsat_atom_num).xyz();
	
	bool bunsat_is_donor(bunsat_rsd.atom_type(
			static_cast< int >(bunsat_atom_num)).is_donor());
	
	bool bunsat_is_acceptor(bunsat_rsd.atom_type(
			static_cast< int >(bunsat_atom_num)).is_acceptor());
	
	bool pass_checks = true;
	
	if (bunsat_is_donor) {
		for(std::set<Size>::const_iterator test_res_it = neighbors.value().begin();
				test_res_it != neighbors.value().end(); ++test_res_it) {
			Size const test_resi = *test_res_it;
			
			std::string test_res_debug = pose.pdb_info()->pose2pdb(test_resi);
			//TR << "checking if donates to " << test_res_debug << std::endl;
			
			if (!bunsat_donor_nbr_residue_check(
					pose, bunsat_candidate_atom_id, bunsat_rsd, bunsat_xyz, test_resi)) {
				pass_checks = false;
				break;
			}
		}
	}
	if (bunsat_is_acceptor) {
		for(std::set<Size>::const_iterator test_res_it = neighbors.value().begin();
				test_res_it != neighbors.value().end(); ++test_res_it) {
			Size const test_resi = *test_res_it;
			if (!bunsat_acc_nbr_residue_check(
					pose, bunsat_candidate_atom_id, bunsat_rsd, bunsat_xyz, test_resi)) {
				pass_checks = false;
				break;
			}
		}
	}
	
	if (pass_checks) all_bur_unsat_polars_++;
	return pass_checks;
}

/*
 return false if a plausible hbond between bunsat_candidate_atom_id as a donor
 and test_resi as an acceptor is found, otherwise return true
 */
bool
BuriedUnsatisfiedPolarsCalculator2::bunsat_donor_nbr_residue_check(
	pose::Pose const & pose,
	id::AtomID const & bunsat_candidate_atom_id,
	Residue const & bunsat_rsd,
	Vector const & bunsat_xyz,
	Size const test_resi
) const
{
	Size const bunsat_resi = bunsat_candidate_atom_id.rsd();
	Size const bunsat_atom_num = bunsat_candidate_atom_id.atomno();
	std::string const bunsat_atom_name =
		bunsat_rsd.atom_name(static_cast< int >(bunsat_atom_num));
	
	Residue const & test_rsd(pose.residue(test_resi));
	
	AtomIndices const & test_accpt_pos = test_rsd.accpt_pos();
	for(vector1< Size >::const_iterator test_atom_it = test_accpt_pos.begin();
		test_atom_it < test_accpt_pos.end(); ++test_atom_it)
	{
		Size const test_atom_num = *test_atom_it;
		
		// exclude neighboring backbone-backbone
		std::string const test_atom_name(
				test_rsd.atom_name(static_cast< int >(test_atom_num)));
		if (adjacent_bbbb_check(bunsat_resi, bunsat_atom_name, test_resi, test_atom_name)) {
			continue;
		}
		
		// exclude self sidechain-sidechain
		if (self_scsc(bunsat_rsd, bunsat_resi, bunsat_atom_num, test_rsd, test_resi, test_atom_num)) {
			continue;
		}
		
		Vector test_xyz = test_rsd.atom(test_atom_num).xyz();
		// check for sulfur-bonds
		if (sulphur_bond_check(test_rsd, test_atom_num, bunsat_xyz, test_xyz)) {
			//TR << "HBond donating to a sulphur atom detected" << std::endl;
			return false;
		}
		
		if (test_rsd.atom_type(static_cast< int >(test_atom_num)).is_acceptor()){
			//if(bunsat_xyz.distance(test_xyz) < 3.3) {
			//if (check_AHD_angle(pose, bunsat_candidate_atom_id, id::AtomID( test_atom_num, test_resi)))
			//TR << "Checking if candidate sites donates to acceptor" << std::endl;
			if (don_geom_check(pose, bunsat_rsd, bunsat_atom_num, bunsat_xyz, test_xyz)) {
				//TR << "HBond donating to a standard acceptor detected" << std::endl;
				return false;
			}
		}
		
	}
	
	return true;
}

/*
 return false if a plausible hbond between bunsat_candidate_atom_id as an acceptor
 and test_resi as an donor is found, otherwise return true
 */
bool
BuriedUnsatisfiedPolarsCalculator2::bunsat_acc_nbr_residue_check(
	pose::Pose const & pose,
	id::AtomID const & bunsat_candidate_atom_id,
	Residue const & bunsat_rsd,
	Vector const & bunsat_xyz,
	Size const & test_resi
) const
{
	Size const bunsat_resi = bunsat_candidate_atom_id.rsd();
	Size const bunsat_atom_num = bunsat_candidate_atom_id.atomno();
	std::string const bunsat_atom_name =
	bunsat_rsd.atom_name(static_cast< int >(bunsat_atom_num));
	
	Residue const & test_rsd(pose.residue(test_resi));
	
	AtomIndices const & test_hpos_polar = test_rsd.Hpos_polar();
	for(vector1< Size>::const_iterator test_atom_it = test_hpos_polar.begin();
			test_atom_it < test_hpos_polar.end(); ++test_atom_it) {
		Size const test_atom_num = *test_atom_it;
		// exclude neighboring backbone-backbone
		std::string const test_atom_name(
				test_rsd.atom_name(static_cast< int >(test_atom_num)));
		if (adjacent_bbbb_check(bunsat_resi, bunsat_atom_name, test_resi, test_atom_name)) {
			continue;
		}
		
		// exclude self sidechain-sidechain
		if (self_scsc(bunsat_rsd, bunsat_resi, bunsat_atom_num, test_rsd, test_resi, test_atom_num)) {
			continue;
		}
		
		Vector test_xyz = test_rsd.atom(test_atom_num).xyz();
		
		// check if coordinating metal
		if (metal_check(test_rsd, bunsat_xyz, test_xyz)) {
			//TR << "HBond accepting from a metal ion detected" << std::endl;
			return false;
		}
		if (test_rsd.atom_type(static_cast< int >(test_atom_num)).is_polar_hydrogen()) {
			//if (check_AHD_angle(pose, bunsat_candidate_atom_id, id::AtomID( test_atom_num, test_resi)))
			if (acc_geom_check(pose, bunsat_xyz, test_rsd, test_atom_num, test_xyz)) {
				//TR << "HBond accepting from a standard donor detected" << std::endl;
				return false;
			}
		}
	}
	
	return true;
}

/*
 return true if bunsat_xyz within h-bond distance to test_xyz and test_rsd is a metal ion,
 otherwise return false
 */
bool
BuriedUnsatisfiedPolarsCalculator2::metal_check(
	Residue const & test_rsd,
	Vector const & bunsat_xyz,
	Vector const & test_xyz
) const
{
	std::string test_rsd_name = test_rsd.name();
	if(test_rsd_name == "CA" || test_rsd_name == "MG" || test_rsd_name == "ZN" ||
		   test_rsd_name == "FE" || test_rsd_name == "MN" || test_rsd_name == "NA") {
		return (bunsat_xyz.distance(test_xyz) < metal_dist_cutoff_); //metal_dist_cutoff
	} else {
		return false;
	}
}

/*
 return true if bunsat_atom_name on bunsat_resi and test_atom_name on test_resi
 are adjacent backbone atoms, otherwise false
 */
bool
BuriedUnsatisfiedPolarsCalculator2::adjacent_bbbb_check(
	Size const & bunsat_resi,
	std::string const & bunsat_atom_name,
	Size const & test_resi,
	std::string const & test_atom_name
) const
{
	bool adjacent = abs(static_cast< int >(bunsat_resi - test_resi)) <= 1;
	bool bunsat_bb = (bunsat_atom_name == "H" || bunsat_atom_name == "O");
	bool test_bb = (test_atom_name == "H" && test_atom_name == "O");
	return (adjacent && bunsat_bb && test_bb);
}

/*
 return true if bunsat_atom_num on bunsat_resi and test_atom_num on test_resi
 are both SC atoms in the same residue, otherwise false
 */
bool
BuriedUnsatisfiedPolarsCalculator2::self_scsc(
	Residue const & bunsat_rsd,
	Size const & bunsat_resi,
	Size const & bunsat_atom_num,
	Residue const & test_rsd,
	Size const & test_resi,
	Size const & test_atom_num
) const
{
	bool same_res = bunsat_resi == test_resi;
	bool bunsat_sc = bunsat_atom_num >= bunsat_rsd.first_sidechain_atom();
	bool test_sc = test_atom_num >= test_rsd.first_sidechain_atom();
	return (same_res && bunsat_sc && test_sc);
}

/*
 return true if bunsat_xyz is within hbond distance of test_xyz and test_xyz
 is a sulphur atom
 */
bool
BuriedUnsatisfiedPolarsCalculator2::sulphur_bond_check(
	Residue const & test_rsd,
	Size const & test_atom_num,
	Vector const & bunsat_xyz,
	Vector const & test_xyz
) const
{
	if(test_rsd.atom_type(static_cast< int >(test_atom_num)).element() == "S") {
		return(bunsat_xyz.distance(test_xyz) < sulph_dist_cutoff_); //suphur_dist_cutoff
	} else {
		return false;
	}
}

/*
 return true if bunsat_atom_num in bunsat_rsd as an hbond donor to an acceptor
 at test_xyz has acceptable hbond geometry, otherwise return false
 */
bool
BuriedUnsatisfiedPolarsCalculator2::don_geom_check(
	core::pose::Pose const &,
	Residue const & bunsat_rsd,
	Size const & bunsat_atom_num,
	Vector const & bunsat_xyz,
	Vector const & test_xyz
) const
{
	for (Size bunsat_H_atom_num = bunsat_rsd.attached_H_begin(static_cast< int >(bunsat_atom_num));
			 bunsat_H_atom_num <= bunsat_rsd.attached_H_end(static_cast< int >(bunsat_atom_num));
			 ++bunsat_H_atom_num) {
		//TR << "checking site attached H for donor geom" << std::endl;
		Vector bunsat_H_xyz = bunsat_rsd.atom(bunsat_H_atom_num).xyz();
		Real AHD_angle = numeric::conversions::degrees(
				angle_of(bunsat_xyz, bunsat_H_xyz, test_xyz));
		Real dist = bunsat_H_xyz.distance(test_xyz);
		//TR << "bunsat_xyz = " << bunsat_xyz[0] << "," << bunsat_xyz[1] << "," << bunsat_xyz[2] << "; bunsat_H_xyz = " << bunsat_H_xyz[0] << "," << bunsat_H_xyz[1] << "," << bunsat_H_xyz[2] << "; test_xyz = " << test_xyz[0] << "," << test_xyz[1] << "," << test_xyz[2] << std::endl;
		//TR << "AH dist = " << dist << ", AHD angle = " << AHD_angle << std::endl;
		//if(dist >= option[OptionKeys::bunsat_calc2::dist_cutoff]) TR << "failed dist check" << std::endl;
		//if (AHD_angle <= option[OptionKeys::bunsat_calc2::AHD_cutoff]) TR << "failed AHD check" << std::endl;
		if (dist < dist_cutoff_ && AHD_angle > AHD_cutoff_) { //dist_cutoff and AHD_cutoff
			return true;
		} else if ((bunsat_rsd.name() == "SER" || bunsat_rsd.name3() == "THR" || bunsat_rsd.name3() == "TYR")
				   && bunsat_atom_num >= bunsat_rsd.first_sidechain_atom())
		{
			/*
			 HXL hydrogen placement is ambiguous, for serine or threonine do
			 a second check based only on heavy atom distance
			 */
			Real dist = bunsat_xyz.distance(test_xyz);
			return dist < hxl_dist_cutoff_; // hydroxyl_dist_cutoff
		}
	}
	return false;
}
	

/// @brief return true if bunsat_atom_num in bunsat_rsd as an hbond acceptor to a donor
/// at test_xyz has acceptable hbond geometry, otherwise return false
bool
BuriedUnsatisfiedPolarsCalculator2::acc_geom_check(
	core::pose::Pose const &,
	Vector const & bunsat_xyz,
	Residue const & test_rsd,
	Size const & test_atom_num,
	Vector const & test_xyz
) const
{
	Vector test_base_xyz = test_rsd.atom(
			(test_rsd.atom_base(static_cast< int >(test_atom_num)))).xyz();
	Real AHD_angle = numeric::conversions::degrees(
			angle_of(bunsat_xyz, test_xyz, test_base_xyz));
	Real dist = bunsat_xyz.distance(test_xyz);
	//dist_cutoff and AHD_cutoff
	if (dist < dist_cutoff_
			&& AHD_angle > AHD_cutoff_) {
		return true;
	} else {
		/*
		 HXL hydrogen placement is ambiguous, for serine or threonine do
		 a second check based only on heavy atom distance
		 */
		if ((test_rsd.name() == "SER" || test_rsd.name3() == "THR" || test_rsd.name3() == "TYR")
				&& test_atom_num >= test_rsd.first_sidechain_atom()) {
			Real dist = bunsat_xyz.distance(test_base_xyz);
			return dist < hxl_dist_cutoff_; // hydroxyl_dist_cutoff
		} else {
			return false;
		}
	}
}

void
BuriedUnsatisfiedPolarsCalculator2::show() {
	std::stringstream sstream;
	sstream << "Buns2calc with name_of_weak_bunsat_calc_ = "
		<< name_of_weak_bunsat_calc_ << ", layered_sasa_ = "
		<< layered_sasa_ << ", generous_hbonds_ = "
		<< generous_hbonds_ << ", sasa_burial_cutoff_ = "
		<< sasa_burial_cutoff_ << ", AHD_cutoff_ = "
		<< AHD_cutoff_ << ", dist_cutoff_ = "
		<< dist_cutoff_ << ".";
	TR << sstream.str() << std::endl;
}

} // namespace buns
} // namespace devel
