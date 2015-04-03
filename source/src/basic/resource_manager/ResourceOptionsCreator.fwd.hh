// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/resource_manager/ResourceOptionsCreator.fwd.hh
/// @brief
/// @author

#ifndef INCLUDED_basic_resource_managerResourceOptionsCreator_FWD_HH
#define INCLUDED_basic_resource_managerResourceOptionsCreator_FWD_HH

//utility headers
#include <utility/pointer/owning_ptr.hh>

namespace basic {
namespace resource_manager {

class ResourceOptionsCreator;
typedef utility::pointer::shared_ptr< ResourceOptionsCreator > ResourceOptionsCreatorOP;
typedef utility::pointer::shared_ptr< ResourceOptionsCreator const > ResourceOptionsCreatorCOP;


} // namespace resource_manager
} // namespace basic


#endif //INCLUDED_basic_resource_managerResourceOptionsCreator_hh
