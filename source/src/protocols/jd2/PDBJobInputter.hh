// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/PDBJobInputter.hh
/// @brief  header file for PDBJobInputter class, part of August 2008 job distributor as planned at RosettaCon08
/// @author Steven Lewis smlewi@gmail.com


#ifndef INCLUDED_protocols_jd2_PDBJobInputter_hh
#define INCLUDED_protocols_jd2_PDBJobInputter_hh

//unit headers
#include <protocols/jd2/JobInputter.hh>
#include <protocols/jd2/PDBJobInputter.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>

//utility headers

namespace protocols {
namespace jd2 {

/// @details This is the simplest implementation of JobInputter, which reads from -s/-l and PDB files.
class PDBJobInputter : public protocols::jd2::JobInputter
{
public:

	PDBJobInputter();

	virtual ~PDBJobInputter();

	/// @brief this function is responsible for filling the pose reference with the pose indicated by the job.  The Job object (within its InnerJob) contains a PoseCOP.  This function needs to either fill the pose reference from the InnerJob or, on first demand of a pose from that InnerJob, instantiate the pose, hand off a COP to the InnerJob, and fill the reference.  This implementation uses pose_from_pdb
 	virtual void pose_from_job( core::pose::Pose & pose, JobOP job );

	/// @brief this function determines what jobs exist from -s/-l
	virtual void fill_jobs( Jobs & jobs );

	/// @brief Return the type of input source that the PDBJobInputter is currently
	///  using.
	/// @return Always <em>PDB_FILE</em>.
	virtual JobInputterInputSource::Enum input_source() const;

}; // PDBJobInputter

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_PDBJobInputter_HH
