// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody_design/LeGrandSasa.fwd.hh
/// @brief 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_scoring_sasa_LEGRANDSASA_FWD_HH
#define INCLUDED_core_scoring_sasa_LEGRANDSASA_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace sasa {

class LeGrandSasa;

	typedef utility::pointer::owning_ptr< LeGrandSasa> LeGrandSasaOP;
	typedef utility::pointer::owning_ptr< LeGrandSasa const> LeGrandSasaCOP; 
}
}
}

#endif	//#ifndef INCLUDED_protocols/antibody_design_LEGRANDSASA_FWD_HH

