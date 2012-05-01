// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/LHRepulsiveRamp.cc
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#include <protocols/antibody2/LHRepulsiveRamp.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <protocols/antibody2/AntibodyInfo.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/import_pose/import_pose.hh>


#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>



#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/OptH.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>


#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
using namespace ObjexxFCL::fmt;

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/docking/SidechainMinMover.hh>
#include <protocols/moves/JumpOutMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/antibody2/AntibodyUtil.hh>
#include <protocols/moves/PyMolMover.hh>

#include <protocols/docking/DockMCMCycle.hh>
#include <protocols/docking/util.hh>






#include <core/chemical/VariantType.hh>
//JQX:: this header file took care of the "CUTPOINT_LOWER" options below

// TODO: make a general DockRepulsiveRamp mover without antibody info and implement into Docking


using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer TR("protocols.antibody2.LHRepulsiveRamp");


using namespace core;
namespace protocols {
namespace antibody2 {


    
// default constructor
LHRepulsiveRamp::LHRepulsiveRamp() : Mover() {

}


    
LHRepulsiveRamp::LHRepulsiveRamp(AntibodyInfoOP antibody_in) : Mover() {
    user_defined_ = false;
    
    init(antibody_in,false);
}
    
LHRepulsiveRamp::LHRepulsiveRamp(AntibodyInfoOP antibody_in, bool camelid) : Mover() {
    user_defined_ = false;
    
    init(antibody_in, camelid);
}
    
    
LHRepulsiveRamp::LHRepulsiveRamp( AntibodyInfoOP antibody_in,
                                   core::scoring::ScoreFunctionCOP dock_scorefxn,
                                   core::scoring::ScoreFunctionCOP pack_scorefxn ) : Mover(){
    user_defined_ = true;
    dock_scorefxn_ = new core::scoring::ScoreFunction(*dock_scorefxn);
    pack_scorefxn_ = new core::scoring::ScoreFunction(*pack_scorefxn);
    
    init(antibody_in, false);
}
    
    
// default destructor
LHRepulsiveRamp::~LHRepulsiveRamp() {}
    
//clone
protocols::moves::MoverOP LHRepulsiveRamp::clone() const {
    return( new LHRepulsiveRamp() );
}
    
    

    
    
void LHRepulsiveRamp::init(AntibodyInfoOP antibody_in, bool camelid ) 
{
    set_default();
    ab_info_ = antibody_in;
    is_camelid_ = camelid;
    
    tf_ = NULL;
    
}
    
    
void LHRepulsiveRamp::set_default(){
    benchmark_       = false;
    use_pymol_diy_   = false;

    rep_ramp_cycles_ = 3 ;
    rot_mag_         = 2.0 ;
    trans_mag_       = 0.1 ;
    num_repeats_     = 4;
    
    if(!user_defined_){
        dock_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
            dock_scorefxn_->set_weight( core::scoring::chainbreak, 1.0 );
            dock_scorefxn_->set_weight( core::scoring::overlap_chainbreak, 10./3. ); 
        pack_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "standard" );
    }

}
    
    
    
    
std::string LHRepulsiveRamp::get_name() const {
    return "LHRepulsiveRamp";
}

    
    

    
    
    
///////////////////////////////////////////////////////////////////////////
/// @begin repulsive_ramp
///
/// @brief ramping up the fullatom repulsive weight slowly to allow the
///        partners to relieve clashes and make way for each other
///
/// @detailed This routine is specially targetted to the coupled
///           optimization of docking partners and the loop region.  The
///           loop modelling & all previous  steps  involve mainly
///           centroid  mode .On switching  on fullatom mode, one is bound
///           to end up with clashes.To relieve the clashes, it is
///           essential to slowly  dial up the  repulsive weight instead of
///           turning it on to the maximum  value in one single step
///
/// @param[in] input pose which is assumed to have a docking fold tree
///
/// @global_read fa_rep : fullatom repulsive weight
///
/// @global_write fa_rep ( It is reset to the original value at the end )
///
/// @remarks A particular portion is  commented out,which can be
///          uncommented if one  uses a  low  resolution  homology  model.
///          Check details in the beginning of the commented out region
///
/// @references
///
/// @authors Aroop 07/13/2010    
///
/// @last_modified 07/13/2010
///////////////////////////////////////////////////////////////////////////
    
void LHRepulsiveRamp::apply( pose::Pose & pose ) {
    TR<<"start apply function ..."<<std::endl;
    
    // remove cutpoints variants for all cdrs
    // "true" forces removal of variants even from non-cutpoints
    loops::remove_cutpoint_variants( pose, true );
    using namespace core::chemical;
    for ( loops::Loops::const_iterator it = ab_info_->all_cdr_loops_.begin(),
         it_end = ab_info_->all_cdr_loops_.end();	it != it_end; ++it ) {
        core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, it->cut() );
        core::pose::add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER,it->cut()+1);
    }
    //TODO: JQX don't understand what the above is doing, why remove that?
    //      It seems it re-added all the cutpoints to all the loops, so that the chain-break scoring function can work.
    
    // add scores to map
    ( *dock_scorefxn_ )( pose );
    
    // dampen fa_rep weight
    core::Real rep_weight_max = dock_scorefxn_->get_weight( core::scoring::fa_rep );

    if( benchmark_ ) {
        rep_ramp_cycles_ = 1;
        num_repeats_ = 1;
    }
    

    
    core::Real rep_ramp_step = (rep_weight_max - 0.02) / core::Real(rep_ramp_cycles_-1);
    core::scoring::ScoreFunctionOP temp_dock_scorefxn =new core::scoring::ScoreFunction( *dock_scorefxn_);
    
    for ( Size i = 1; i <= rep_ramp_cycles_; i++ ) {
        core::Real rep_weight = 0.02 + rep_ramp_step * Real(i-1);
        TR<<"   repulsive ramp cycle "<<i<<":     rep_weight = "<<rep_weight<<std::endl;
        temp_dock_scorefxn->set_weight( core::scoring::fa_rep, rep_weight );
        

        docking::DockMCMCycleOP dockmcm_cyclemover = new docking::DockMCMCycle( ab_info_->LH_dock_jump(), temp_dock_scorefxn, pack_scorefxn_ );
        //TODO: print scoring function in apply and move "new" out
            dockmcm_cyclemover->set_rot_magnitude(rot_mag_);
            dockmcm_cyclemover->set_task_factory(tf_);
            dockmcm_cyclemover->set_move_map(cdr_dock_map_);
        
        for (Size j=1; j<=num_repeats_; j++) {
            dockmcm_cyclemover -> apply(pose);
            TR<<"       doing rb_mover_min_trial in the DockMCMCycle  ...   "<<j<<std::endl;
            if(use_pymol_diy_) pymol_->apply(pose);
        }
        dockmcm_cyclemover -> reset_cycle_index(); //JQX: only do the rb_mover_min_trial (index<7)
        dockmcm_cyclemover -> get_mc()->recover_low( pose ); //chose the pose having the lowest score

    }

    TR<<"finish apply function !!!"<<std::endl;

}
    

    
    
void LHRepulsiveRamp::set_task_factory(pack::task::TaskFactoryCOP tf){        
    tf_ = new pack::task::TaskFactory(*tf);
}
    
void LHRepulsiveRamp::set_move_map(kinematics::MoveMapCOP cdr_dock_map){
    cdr_dock_map_ = new kinematics::MoveMap(*cdr_dock_map);
}


} // namespace antibody2
} // namespace protocols





