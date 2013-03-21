// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Add constraints to the current pose conformation.
/// @author Yifan Song

#ifndef INCLUDED_protocols_hybridization_HybridizeProtocol_hh
#define INCLUDED_protocols_hybridization_HybridizeProtocol_hh

#include <protocols/hybridization/HybridizeProtocol.fwd.hh>
#include <protocols/hybridization/FoldTreeHybridize.hh>

#include <protocols/moves/Mover.hh>

#include <protocols/loops/Loops.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <utility/file/FileName.hh>

namespace protocols {
namespace hybridization {

class HybridizeProtocol : public protocols::moves::Mover {

public:
	HybridizeProtocol();
	//HybridizeProtocol(std::string template_list_file);

	void init();

	void read_template_structures(utility::file::FileName template_list);
	void read_template_structures(utility::vector1 < utility::file::FileName > const & template_filenames);

	void add_template(
		std::string template_fn,
		std::string cst_fn,
		std::string symmdef_file = "NULL",
		core::Real weight = 1.,
		core::Real domain_assembly_weight = 0.,
		core::Size cluster_id = 1,
		utility::vector1<core::Size> cst_reses = utility::vector1<core::Size>(0) );

	void add_template(
		core::pose::PoseOP template_pose,
		std::string cst_fn,
		std::string symmdef_file = "NULL",
		core::Real weight = 1.,
		core::Real domain_assembly_weight = 0.,
		core::Size cluster_id = 1,
		utility::vector1<core::Size> cst_reses = utility::vector1<core::Size>(0),
		std::string filename="default" );

	void pick_starting_template(core::Size & initial_template_index,
		core::Size & initial_template_index_icluster,
		utility::vector1 < core::Size > & template_index_icluster,
		utility::vector1 < core::pose::PoseOP > & templates_icluster,
		utility::vector1 < core::Real > & weights_icluster,
		utility::vector1 < protocols::loops::Loops > & template_chunks_icluster,
		utility::vector1 < protocols::loops::Loops > & template_contigs_icluster);

	utility::vector1 <Loops>
	expand_domains_to_full_length(utility::vector1 < utility::vector1 < Loops > > all_domains, Size ref_domains_index, Size n_residues);

	void
	align_by_domain(utility::vector1<core::pose::PoseOP> & poses, utility::vector1 < Loops > domains, core::pose::PoseOP & ref_pose);

	void
	align_by_domain(core::pose::Pose & pose, core::pose::Pose const & ref_pose, utility::vector1 <Loops> domains);

	//fpd optionally do not hybridize in stage 1
	void
	initialize_and_sample_loops(
		core::pose::Pose &pose,
		core::pose::PoseOP chosen_templ,
		protocols::loops::Loops template_contigs_icluster,
		core::scoring::ScoreFunctionOP scorefxn);

	// check fragments ... if they do not exist dynamically allocate them
	void check_and_create_fragments( Pose & );

	virtual void apply( Pose & );
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	virtual void
	parse_my_tag( TagPtr const, DataMap &, Filters_map const &, Movers_map const &, Pose const & );

	// set options
	void set_batch_relax( core::Size newval ) { batch_relax_ = newval; }
	void add_big_fragments( core::fragment::FragSetOP newval ) { fragments_big_.push_back(newval); }
	void add_small_fragments( core::fragment::FragSetOP newval ) { fragments_small_.push_back(newval); }
	void set_stage1_scorefxn( core::scoring::ScoreFunctionOP newval ) { stage1_scorefxn_ = newval; }
	void set_stage2_scorefxn( core::scoring::ScoreFunctionOP newval ) { stage2_scorefxn_ = newval; }
	void set_fullatom_scorefxn( core::scoring::ScoreFunctionOP newval ) { fa_scorefxn_ = newval; }

private:
	// parsible options
	utility::vector1 < core::Size > starting_templates_;
	core::Real stage1_probability_, stage1_increase_cycles_, stage2_increase_cycles_;
	core::Size stage1_1_cycles_, stage1_2_cycles_, stage1_3_cycles_, stage1_4_cycles_;
	// 1mer fragment insertion weight where fragments are not allowed (across anchors) , vs. chunk insertion + big and small frags
	core::Real frag_1mer_insertion_weight_;
	// small fragment insertion weight where big fragments are not allowed (across anchors) , vs. chunk insertion + big frags
	core::Real small_frag_insertion_weight_;
	core::Real big_frag_insertion_weight_; // fragment insertion weight, vs. chunk insertion + small gap frags
	core::Real frag_weight_aligned_; // fragment insertion to the aligned region, vs. unaligned region
	bool auto_frag_insertion_weight_; // automatically set fragment insertion weights
	core::Size max_registry_shift_;
	bool domain_assembly_, add_hetatm_, realign_domains_, realign_domains_stage2_, add_non_init_chunks_, no_global_frame_, linmin_only_;
	core::Real hetatm_cst_weight_;
	core::scoring::ScoreFunctionOP stage1_scorefxn_, stage2_scorefxn_, fa_scorefxn_;
	std::string fa_cst_fn_;
	std::string disulf_file_;
	core::Size cartfrag_overlap_;

	// ddomain options
	core::Real pcut_,hcut_;
	core::Size length_;

	// relax
	core::Size batch_relax_, relax_repeats_;

	// abinitio frag9,frag3 flags
	utility::vector1 <core::fragment::FragSetOP> fragments_big_;  // 9mers/fragA equivalent in AbrelaxApplication
	utility::vector1 <core::fragment::FragSetOP> fragments_small_; // 3mers/fragB equivalent in AbrelaxApplication

	// native pose, aln
	core::pose::PoseOP native_;
	core::sequence::SequenceAlignmentOP aln_;

	// template information (all from hybrid.config)
	utility::vector1 < core::pose::PoseOP > templates_;
	utility::vector1 < std::string > template_fn_;
	utility::vector1 < std::string > template_cst_fn_;
	utility::vector1 < std::string > symmdef_files_;
	utility::vector1 < core::Real > template_weights_;
	utility::vector1 < core::Real > domain_assembly_weights_;
	utility::vector1 < core::Size > template_clusterID_;
	utility::vector1 < protocols::loops::Loops > template_chunks_;
	utility::vector1 < protocols::loops::Loops > template_contigs_;
	utility::vector1 < utility::vector1<core::Size> > template_cst_reses_;
	core::Real template_weights_sum_;
	std::map< Size, utility::vector1 < Size > > clusterID_map_;

	utility::vector1< protocols::loops::Loops > domains_;

	// strand pairings
	std::string pairings_file_;
	utility::vector1<core::Size> sheets_;
	utility::vector1<core::Size> random_sheets_;
	bool filter_templates_;
	utility::vector1< std::pair< core::Size, core::Size > > strand_pairs_;

};

} // hybridization
} // protocols

#endif
