// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_DeNovo_Protocol.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_rna_RNA_DeNovoProtocol_HH
#define INCLUDED_protocols_rna_RNA_DeNovoProtocol_HH

#include <core/types.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/farna/RNA_FragmentMonteCarlo.fwd.hh>
#include <protocols/farna/RNA_StructureParameters.hh>
#include <protocols/farna/RNA_DeNovoProtocolOptions.fwd.hh>

//Oooh.
#include <ObjexxFCL/FArray1D.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#include <list>

#include <core/io/silent/silent.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/moves/MonteCarlo.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace farna {

/// @brief The RNA de novo structure modeling protocol
class RNA_DeNovoProtocol: public protocols::moves::Mover {
public:

	/// @brief Construct the protocol object given
	/// the RNA fragment library to use.
	RNA_DeNovoProtocol( RNA_DeNovoProtocolOptionsCOP options = 0 );

	~RNA_DeNovoProtocol();

	/// @brief Clone this object
	virtual protocols::moves::MoverOP clone() const;

	/// @brief Apply the RNA denovo modeling protocol to the input pose
	void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

	virtual void show(std::ostream & output=std::cout) const;

	void
	output_to_silent_file( core::pose::Pose & pose, std::string const & silent_file, std::string const & out_file_tag,
		bool const score_only = false ) const;

	void
	align_and_output_to_silent_file( core::pose::Pose & pose, std::string const & silent_file, std::string const & out_file_tag ) const;

	void
	set_refine_pose_list( utility::vector1<core::pose::PoseOP> const & setting ) { refine_pose_list_ = setting; };

private:

	void
	initialize_lores_silent_file();

	void
	initialize_constraints( core::pose::Pose & pose );

	void
	initialize_scorefxn( core::pose::Pose & pose );

	void
	initialize_tag_is_done();

	void
	output_silent_struct( core::io::silent::SilentStruct & s,
		core::io::silent::SilentFileData & silent_file_data,
		std::string const & silent_file,
		core::pose::Pose & pose,
		std::string const out_file_tag,
		bool const score_only = false ) const;

	void
	add_number_base_pairs( core::pose::Pose const & pose_input, core::io::silent::SilentStruct & s ) const;

	void
	add_number_native_base_pairs( core::pose::Pose & pose, core::io::silent::SilentStruct & s ) const;

	void
	check_for_loop_modeling_case( std::map< core::id::AtomID, core::id::AtomID > & atom_id_map, core::pose::Pose const & pose ) const;

	void
	calc_rmsds( core::io::silent::SilentStruct & s, core::pose::Pose & pose, std::string const & out_file_tag ) const;

	void
	add_chem_shift_info(core::io::silent::SilentStruct & silent_struct, core::pose::Pose const & const_pose) const;

private:

	RNA_DeNovoProtocolOptionsCOP options_;

	std::string lores_silent_file_;

	RNA_FragmentMonteCarloOP rna_fragment_monte_carlo_;

	std::map< std::string, bool > tag_is_done_;

	core::scoring::ScoreFunctionOP denovo_scorefxn_;
	core::scoring::ScoreFunctionOP hires_scorefxn_;

	utility::vector1<core::pose::PoseOP> refine_pose_list_;

}; // class RNA_DeNovoProtocol

std::ostream &operator<< ( std::ostream &os, RNA_DeNovoProtocol const &mover );

} //farna
} //protocols

#endif
