// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


/// @file protocols/antibody/Ab_ModelCDRH3.fwd.hh
/// @brief Build a homology model of an antibody
/// @detailed
///
///// @author Jianqing Xu (xubest@gmail.com)
//

#ifndef INCLUDED_protocols_antibody_ModelCDRH3_fwd_hh
#define INCLUDED_protocols_antibody_ModelCDRH3_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace antibody {

class ModelCDRH3;
typedef utility::pointer::owning_ptr< ModelCDRH3 > ModelCDRH3OP;
typedef utility::pointer::owning_ptr< ModelCDRH3 const > ModelCDRH3COP;

} // antibody
} // protocols

#endif
