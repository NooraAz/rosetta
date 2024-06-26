// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/pose_length_moves/InsertResMover.hh
/// @brief inserts ideal residues into pose. Useful for extending helices
///
/// @author TJ Brunette tjbrunette@gmail.com

#ifndef INCLUDED_protocols_pose_length_moves_InsertResMover_hh
#define INCLUDED_protocols_pose_length_moves_InsertResMover_hh


#include <protocols/moves/Mover.hh>

#include <protocols/pose_length_moves/InsertResMover.fwd.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>

// C++ Headers
#include <string>
// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.fwd.hh>
#include <numeric/xyzVector.hh>


namespace protocols {
namespace pose_length_moves {


class InsertResMover : public protocols::moves::Mover {
public:
	InsertResMover();
	numeric::xyzVector<core::Real> center_of_mass(core::pose::Pose const & pose);
	void extendRegion(core::pose::PoseOP poseOP, core::Size chain_id, core::Size length);
	void apply( Pose & pose ) override;
	moves::MoverOP clone() const override { return utility::pointer::make_shared< InsertResMover >( *this ); }
	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & datamap ) override;
	core::pose::PoseOP get_additional_output() override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	utility::vector1<core::pose::PoseOP> posesToOutput_;
	utility::vector1<bool> posesOutputed_;
	std::string chain_;
	core::Size residue_;
	bool grow_toward_Nterm_;
	bool ideal_;
	core::Size lowAddRes_;
	core::Size highAddRes_;
	std::string resType_;
	bool useInputAngles_;
	core::Real phi_;
	core::Real psi_;
	core::Real omega_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::Size steal_angles_from_res_;
};


} // pose_length_moves
} // protocols

#endif
