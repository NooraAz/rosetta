// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Yuan Liu (wendao@uw.edu)

#ifndef INCLUDED_core_pack_dunbrack_cenrot_SingleResidueCenrotLibrary_fwd_hh
#define INCLUDED_core_pack_dunbrack_cenrot_SingleResidueCenrotLibrary_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace pack {
namespace dunbrack {
namespace cenrot {

class SingleResidueCenrotLibrary;

typedef utility::pointer::owning_ptr< SingleResidueCenrotLibrary > SingleResidueCenrotLibraryOP;
typedef utility::pointer::access_ptr< SingleResidueCenrotLibrary > SingleResidueCenrotLibraryAP;
typedef utility::pointer::owning_ptr< SingleResidueCenrotLibrary const > SingleResidueCenrotLibraryCOP;
typedef utility::pointer::access_ptr< SingleResidueCenrotLibrary const > SingleResidueCenrotLibraryCAP;

} // cenrot
} // dunbrack
} // pack
} // core

#endif // INCLUDED_core_pack_dunbrack_cenrot_SingleResidueCenrotLibrary_HH

