// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/design_opt/PointMutationCalculator.hh
/// @author Chris King (chrisk1@uw.edu)

#ifndef INCLUDED_protocols_design_opt_PointMutationCalculator_hh
#define INCLUDED_protocols_design_opt_PointMutationCalculator_hh
#include <protocols/design_opt/PointMutationCalculator.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/AA.hh>
#include <utility/tag/Tag.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <utility/vector1.hh>

#ifdef PYROSETTA
	#include <protocols/filters/Filter.hh>
#endif


namespace protocols {
namespace design_opt{

class PointMutationCalculator : public utility::pointer::ReferenceCount
{
public:
	typedef core::pose::Pose Pose;
public:
	PointMutationCalculator();
	PointMutationCalculator(
		core::pack::task::TaskFactoryOP task_factory,
		core::scoring::ScoreFunctionOP scorefxn,
		utility::vector1< protocols::filters::FilterOP > filters,
		protocols::moves::MoverOP relax_mover,
		bool dump_pdb = false,
		std::string sample_type = "low"
	);
	PointMutationCalculator(
		core::pack::task::TaskFactoryOP task_factory,
		core::scoring::ScoreFunctionOP scorefxn,
		protocols::filters::FilterOP filter,
		protocols::moves::MoverOP relax_mover,
		bool dump_pdb = false,
		std::string sample_type = "low"
	);
	virtual ~PointMutationCalculator();

	void calc_point_mut_filters( Pose const & start_pose,
			utility::vector1< std::pair< core::Size, utility::vector1< std::pair< core::chemical::AA, utility::vector1< core::Real > > > > > & seqpos_aa_vals_vec );
	void calc_point_mut_filters( Pose const & start_pose,
			utility::vector1< std::pair< core::Size, utility::vector1< std::pair< core::chemical::AA, core::Real > > > > & seqpos_aa_val_vec );
	protocols::design_opt::PointMutationCalculatorOP clone() const;

	core::pack::task::TaskFactoryOP task_factory() const;
	void task_factory( core::pack::task::TaskFactoryOP task_factory );
	core::scoring::ScoreFunctionOP scorefxn() const;
	void scorefxn( core::scoring::ScoreFunctionOP scorefxn );
	utility::vector1< protocols::filters::FilterOP > filters() const;
	void filters( utility::vector1< protocols::filters::FilterOP > filters );
	protocols::moves::MoverOP relax_mover() const;
	void relax_mover( protocols::moves::MoverOP relax_mover );
	bool dump_pdb() const;
	void dump_pdb( bool const dump_pdb );
	std::string sample_type() const;
	void sample_type( std::string const sample_type );
private:
	core::pack::task::TaskFactoryOP task_factory_;
	core::scoring::ScoreFunctionOP scorefxn_;
	utility::vector1< protocols::filters::FilterOP > filters_;
	protocols::moves::MoverOP relax_mover_;
	std::string sample_type_;
	core::Real flip_sign_;
	bool dump_pdb_;
};


} // design_opt 
} // protocols


#endif /*INCLUDED_protocols_design_opt_PointMutationCalculator_HH*/
