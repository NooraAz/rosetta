// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/task_operations/JointSequenceOperation.cc
/// @brief set designable residues to those observed in a set of structures
/// @author Rocco Moretti, rmoretti@u.washington.edu


// Unit Headers
#include <protocols/toolbox/task_operations/JointSequenceOperation.hh>
#include <protocols/toolbox/task_operations/JointSequenceOperationCreator.hh>

// Project Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/import_pose/import_pose.hh>
#include <core/sequence/Sequence.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <protocols/jd2/util.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <string>

#include <utility/vector0.hh>


static basic::Tracer TR("protocols.toolbox.task_operations.JointSequenceOperation");

namespace protocols{
namespace toolbox{
namespace task_operations{

core::pack::task::operation::TaskOperationOP
JointSequenceOperationCreator::create_task_operation() const
{
	return new JointSequenceOperation;
}

/// @brief default constructor
JointSequenceOperation::JointSequenceOperation():
	TaskOperation(),
	use_current_pose_(true),
	use_natro_(false),
	ubr_(0),
	chain_(0)
{
}

/// @brief destructor
JointSequenceOperation::~JointSequenceOperation() {}

/// @brief clone
core::pack::task::operation::TaskOperationOP
JointSequenceOperation::clone() const {
	return new JointSequenceOperation( *this );
}

/// @brief all AA that have a higher probability in the seqprofile
/// than the native residue are allowed. probability also
/// needs to be higher than min_aa_probability_
void
JointSequenceOperation::apply( Pose const & pose, PackerTask & task ) const
{

	core::conformation::symmetry::SymmetryInfoCOP syminfo = NULL;
	if( core::pose::symmetry::is_symmetric(pose) ) syminfo = core::pose::symmetry::symmetry_info(pose);

	core::Size start = 1;
	core::Size end = pose.total_residue();

	if( chain_ > 0 ){ 
		start = pose.conformation().chain_begin( chain_ );
		end = 	pose.conformation().chain_end( chain_ );
	}

	core::Size seq_length = end - start + 1;

	for( std::vector<core::sequence::SequenceOP>::const_iterator iter(sequences_.begin()); iter != sequences_.end(); iter++ ) {
		if( (*iter)->length() != seq_length ) {
				std::string name("current pdb");
				if( pose.pdb_info() ) {
					name = pose.pdb_info()->name();
				}
				TR.Warning << "WARNING: considered sequence " << (*iter)->id() << " contains a different number of residues than " << name << std::endl;
		}
	}
	for( core::Size ii = start; ii <= end; ++ii){
		if( !pose.residue_type( ii ).is_protein() ) continue;
		utility::vector1< bool > allowed(core::chemical::num_canonical_aas, false);

		if(use_current_pose_) {
			if( pose.aa(ii) <= allowed.size() ) allowed[ pose.aa(ii) ] = true;
		}
		if( !core::pose::symmetry::is_symmetric(pose) || syminfo->chi_is_independent(ii) ) {
			for( std::vector<core::sequence::SequenceOP>::const_iterator iter(sequences_.begin()); iter != sequences_.end(); iter++ ) {
				if ( ii > (*iter)->length() ) continue; // ignore short references
				char aa( (*(*iter))[ ii ] );
				if( core::chemical::oneletter_code_specifies_aa(aa) ) {
					if(core::chemical::aa_from_oneletter_code(aa)<=allowed.size()) {
						allowed[ core::chemical::aa_from_oneletter_code(aa)  ] = true;
					}
				}
			}
		}
		task.nonconst_residue_task(ii).restrict_absent_canonical_aas( allowed );
	} //loop over all residues

	if( use_natro_ && ubr_) {
		task.append_rotamerset_operation( ubr_ );
	}
} // apply

void
JointSequenceOperation::parse_tag( TagPtr tag )
{
	use_current_pose( tag->getOption< bool >( "use_current", true ) );
	use_natro( tag->getOption< bool >( "use_natro", false ) );

	// specifiy a chain, if 0 use all
	chain_ = tag->getOption< core::Size >("chain", 0 );

	if( tag->getOption< bool >( "use_native", false )) {
		if( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
			add_native_pdb( basic::options::option[ basic::options::OptionKeys::in::file::native ] );
		}
		else {
			utility_exit_with_message("Native PDB not specified on command line.");
		}
	}

 if( tag->getOption< bool >( "use_starting_as_native", false )) {
    core::pose::PoseCOP pose (protocols::jd2::get_current_jobs_starting_pose());
    add_native_pose( pose );
  }

	if( tag->hasOption("filename") ){
		add_pdb( tag->getOption< String >( "filename" ) );
	}
	if( tag->hasOption("native") ){
		add_native_pdb( tag->getOption< String >( "native" ) );
	}
}

/// @brief Add the sequence from the given filename to the set of allowed aas.
void
JointSequenceOperation::add_pdb( std::string filename )
{
	
	core::pose::Pose new_pose;
	core::pose::Pose pose;

	core::import_pose::pose_from_pdb( new_pose, filename );
	if( chain_ == 0 ){
		TR << "taking only chain " << chain_ << std::endl;
		pose = new_pose.split_by_chain( chain_ ) ;
	}
	else
		pose = new_pose;
	
	add_pose( pose );
}

/// @brief Add the sequence from the given pose to the set of allowed aas.
void
JointSequenceOperation::add_pose( Pose const & pose )
{
	std::string name("unknown");
	if( pose.pdb_info() ) {
		name = pose.pdb_info()->name();
	}
	sequences_.push_back( new core::sequence::Sequence(pose.sequence(), name) );
}

/// @brief Add the sequence from the given filename to the set of allowed aas
/// and add the rotamers to the set of possible rotamers (if use_natro_ is set)
void
JointSequenceOperation::add_native_pdb( std::string filename ) {
	core::pose::PoseOP poseop(new core::pose::Pose);
	core::import_pose::pose_from_pdb( *poseop, filename );
	add_native_pose( poseop );
}

/// @brief Add the sequence from the given pose to the set of allowed aas
/// and add the rotamers to the set of possible rotamers
void
JointSequenceOperation::add_native_pose( core::pose::PoseCOP posecop ){ // PoseCOP needed by UnboundRot, unfortunately
	if( use_natro_ ) { // Deliberate check now to avoid keeping native poses around if we're never going to need them.
		ubr_->add_pose(posecop);
	}
	add_pose(*posecop);
}

/// @brief Should the current pose (pose supplied to apply) be used in addition to the other ones?
void
JointSequenceOperation::use_current_pose( bool ucp )
{
	use_current_pose_ = ucp;
}

/// @brief Should the rotamers for the native poses be used?
void
JointSequenceOperation::use_natro( bool unr ) {
	if( !use_natro_ && unr ) {
		ubr_ = new core::pack::rotamer_set::UnboundRotamersOperation;
	}
	if( use_natro_ && !unr ) {
		ubr_ = 0; // Allow owning pointer to garbage collect as necessary.
	}
	use_natro_ = unr;
}

/// @brief which chain should be considered
void
JointSequenceOperation::set_chain( core::Size chain){
	chain_ = chain;
}

} // TaskOperations
} // toolbox
} // protocols

