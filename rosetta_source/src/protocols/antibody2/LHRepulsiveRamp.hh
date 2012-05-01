// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/LHRepulsiveRamp.hh
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#ifndef INCLUDED_protocols_antibody2_LHRepulsiveRamp_hh
#define INCLUDED_protocols_antibody2_LHRepulsiveRamp_hh





#include <protocols/antibody2/LHRepulsiveRamp.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/antibody2/AntibodyInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/RepeatMover.fwd.hh>
#include <protocols/moves/PyMolMover.fwd.hh>
#include <protocols/docking/DockTaskFactory.fwd.hh>

#ifdef PYROSETTA
#include <protocols/moves/PyMolMover.hh>
#endif


using namespace core;
namespace protocols {
namespace antibody2 {
        
class LHRepulsiveRamp: public moves::Mover {
            
            
public:
    
    /// @brief default constructor
	LHRepulsiveRamp();
    
	/// @brief constructor with arguments
    
    LHRepulsiveRamp(AntibodyInfoOP antibody_in );
    
	LHRepulsiveRamp(AntibodyInfoOP antibody_in, bool camelid );
    
    LHRepulsiveRamp(AntibodyInfoOP antibody_in,
                    core::scoring::ScoreFunctionCOP dock_scorefxn,
                    core::scoring::ScoreFunctionCOP pack_scorefxn );
        
    virtual protocols::moves::MoverOP clone() const;
    
	/// @brief default destructor
	~LHRepulsiveRamp();
    
    void set_default();
    
    void set_dock_score_func(scoring::ScoreFunctionCOP dock_scorefxn ){
        dock_scorefxn_ = new core::scoring::ScoreFunction(*dock_scorefxn);
    }
    
    void set_pack_score_func(scoring::ScoreFunctionCOP pack_scorefxn){
        pack_scorefxn_ = new core::scoring::ScoreFunction(*pack_scorefxn);
    }
    
    virtual void apply( core::pose::Pose & pose );
    
    virtual std::string get_name() const;
    
    
    void set_task_factory(pack::task::TaskFactoryCOP tf);
    void set_move_map(kinematics::MoveMapCOP cdr_dock_map);
    
    core::Real set_rot_mag  (core::Real rot_mag)  {return rot_mag_  =rot_mag;  }
    core::Real set_trans_mag(core::Real trans_mag){return trans_mag_=trans_mag;}
    
    void turn_on_and_pass_the_pymol(moves::PyMolMoverOP pymol){
        use_pymol_diy_ = true;
        pymol_ = pymol;
    }
    
private:

    AntibodyInfoOP ab_info_;
    
    bool user_defined_;
    bool benchmark_;
    bool is_camelid_;
    core::Size rep_ramp_cycles_;
    core::Real rot_mag_;
    core::Real trans_mag_;
    core::Size num_repeats_;

    moves::PyMolMoverOP pymol_;
    bool use_pymol_diy_;
    
    scoring::ScoreFunctionOP dock_scorefxn_;
    scoring::ScoreFunctionOP pack_scorefxn_;
    
    void init(AntibodyInfoOP ab_info, bool camelid);
            
    void snugfit_MC_min(pose::Pose & pose, core::scoring::ScoreFunctionOP  temp_scorefxn);

    
	void repulsive_ramp( pose::Pose & pose_in, loops::Loops loops_in );
    

    
    
	//packer task
    pack::task::TaskFactoryOP tf_;
    kinematics::MoveMapOP cdr_dock_map_;
};
    
    
    
    
    
    
    
} // namespace antibody2
} // namespace protocols

#endif








