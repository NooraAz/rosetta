// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/farna/RNA_FragmentMonteCarloOptions.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/farna/RNA_FragmentMonteCarloOptions.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/database/open.hh>
#include <utility/file/file_sys_util.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.farna.RNA_FragmentMonteCarloOptions" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

namespace protocols {
namespace farna {

//Constructor
RNA_FragmentMonteCarloOptions::RNA_FragmentMonteCarloOptions():
	rounds_( 10 ),
	monte_carlo_cycles_( 0 ), /* will be reset later */
	user_defined_cycles_( false ), /* will change to true if set_monte_carlo_cycles() is called */
	dump_pdb_( false ),
	minimize_structure_( false ),
	relax_structure_( false ),
	temperature_( 2.0 ),
	ignore_secstruct_( false ),
	chainbreak_weight_( -1.0 ), /* use farna/rna_lores.wts number unless user specified. -1.0 is never really used. */
	linear_chainbreak_weight_( -1.0 ),  /* use farna/rna_lores.wts number unless user specified. -1.0 is never really used. */
	close_loops_( true ),
	close_loops_in_last_round_( true ),
	close_loops_after_each_move_( false ),
	allow_bulge_( false ),
	allow_consecutive_bulges_( false ),
	jump_change_frequency_( 0.1 ),
	autofilter_( false ),
	autofilter_score_quantile_( 0.20 ),
	titrate_stack_bonus_( true ),
	move_first_rigid_body_( false ),
	root_at_first_rigid_body_( false ),
	suppress_bp_constraint_( 1.0 ),
	filter_lores_base_pairs_( false ),
	filter_lores_base_pairs_early_( false ),
	filter_chain_closure_( true ),
	filter_chain_closure_distance_( 6.0 ), /* in Angstroms. This is pretty loose!*/
	filter_chain_closure_halfway_( true ),
	staged_constraints_( false ),
	vary_bond_geometry_( false ),
	simple_rmsd_cutoff_relax_( false ),
	refine_from_silent_( false ),
	refine_pose_( false ),
	bps_moves_( false ),
	minimizer_use_coordinate_constraints_( false ),
	minimize_bps_( false ),
	// following is odd, but note that core::scoring::rna::chemical_shift machinery also checks global options system.
	use_chem_shift_data_( option[ OptionKeys::score::rna_chemical_shift_exp_data].user() ),
	verbose_( true ),
	all_rna_fragments_file_( basic::database::full_name("sampling/rna/RICHARDSON_RNA09.torsions") ),
	rna_params_file_( "" ),
	jump_library_file_( basic::database::full_name("sampling/rna/1jj2_RNA_jump_library.dat" ) )
{}

//Destructor
RNA_FragmentMonteCarloOptions::~RNA_FragmentMonteCarloOptions()
{}

/// @brief copy constructor
RNA_FragmentMonteCarloOptions::RNA_FragmentMonteCarloOptions( RNA_FragmentMonteCarloOptions const & src ) :
	ResourceOptions( src )
{
	*this = src;
}

/// @brief clone the options
RNA_FragmentMonteCarloOptionsOP
RNA_FragmentMonteCarloOptions::clone() const
{
	return RNA_FragmentMonteCarloOptionsOP( new RNA_FragmentMonteCarloOptions( *this ) );
}

///////////////////////////////////////////////////////////////////
void
RNA_FragmentMonteCarloOptions::initialize_from_command_line() {

	if ( option[ rna::farna::cycles ].user() ) {
		set_monte_carlo_cycles( option[ rna::farna::cycles ]() );
		set_user_defined_cycles( true );
	}
	minimize_structure_ = option[ rna::farna::minimize_rna ]();
	relax_structure_ = option[ rna::farna::relax_rna ]();
	allow_bulge_ = option[ rna::farna::allow_bulge ]();
	set_temperature( option[ rna::farna::temperature ] );
	set_ignore_secstruct( option[ rna::farna::ignore_secstruct ] );
	set_jump_change_frequency( option[ rna::farna::jump_change_frequency ] );
	set_close_loops( option[ rna::farna::close_loops] );
	set_close_loops_after_each_move( option[ rna::farna::close_loops_after_each_move ] );
	set_dump_pdb( option[ rna::farna::dump ] ) ;
	set_staged_constraints( option[ rna::farna::staged_constraints ] ) ;
	set_filter_lores_base_pairs(  option[ rna::farna::filter_lores_base_pairs] );
	set_filter_lores_base_pairs_early(  option[ rna::farna::filter_lores_base_pairs_early] );
	set_suppress_bp_constraint(  option[ rna::farna::suppress_bp_constraint]() );
	set_filter_chain_closure(  option[ rna::farna::filter_chain_closure ]() );
	set_filter_chain_closure_distance(  option[ rna::farna::filter_chain_closure_distance ]() );
	set_filter_chain_closure_halfway(  option[ rna::farna::filter_chain_closure_halfway ]() );
	set_vary_bond_geometry(  option[ rna::farna::vary_geometry ] );
	set_move_first_rigid_body(  option[ rna::farna::move_first_rigid_body ] );
	set_root_at_first_rigid_body(  option[ rna::farna::root_at_first_rigid_body ] );
	set_autofilter( option[ rna::farna::autofilter ] );
	set_bps_moves( option[ rna::farna::bps_moves ] );
	set_minimizer_use_coordinate_constraints( option[ rna::farna::minimizer_use_coordinate_constraints ]() );
	set_minimize_bps( option[ rna::farna::minimize_bps ]() );

	set_allow_consecutive_bulges( option[ rna::farna::allow_consecutive_bulges ]() ) ;
	set_allowed_bulge_res( option[ rna::farna::allowed_bulge_res ]() ) ;
	set_extra_minimize_res( option[ rna::farna::extra_minimize_res ]() ) ;
	set_extra_minimize_chi_res( option[ rna::farna::extra_minimize_chi_res ]() ) ;

	set_simple_rmsd_cutoff_relax( option[ rna::farna::simple_relax ] );

	std::string const in_path = option[ in::path::path ]()[1];
	if ( option[ rna::farna::vall_torsions ].user() ) {
		// check in database first
		std::string vall_torsions_file = basic::database::full_name("/sampling/rna/" + option[ rna::farna::vall_torsions ]() );
		if ( !utility::file::file_exists( vall_torsions_file ) && !utility::file::file_exists( vall_torsions_file + ".gz" ) )  vall_torsions_file = in_path + option[ rna::farna::vall_torsions ]();
		set_vall_torsions_file( vall_torsions_file );
	}
	if ( option[ rna::farna::use_1jj2_torsions ]() ) set_vall_torsions_file( basic::database::full_name("sampling/rna/1jj2.torsions") );
	if ( option[ rna::farna::jump_library_file ].user() ) set_jump_library_file( in_path + option[ rna::farna::jump_library_file] );
	if ( option[ rna::farna::params_file ].user() ) set_rna_params_file( in_path + option[ rna::farna::params_file ] );

	if ( option[ rna::farna::rna_lores_linear_chainbreak_weight ].user() ) set_linear_chainbreak_weight( option[ rna::farna::rna_lores_linear_chainbreak_weight ]() );
	if ( option[ rna::farna::rna_lores_chainbreak_weight ].user() ) set_chainbreak_weight( option[ rna::farna::rna_lores_chainbreak_weight ]() );

	if ( filter_lores_base_pairs_early_ ) set_filter_lores_base_pairs( true );
}

} //farna
} //protocols
