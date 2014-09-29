// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file numeric/kdtree/Wrappedchar.hh
/// @brief silly little class to wrap a char in an OP.
/// @author James Thompson

#ifndef INCLUDED_numeric_kdtree_WrappedChar_hh
#define INCLUDED_numeric_kdtree_WrappedChar_hh

#include <numeric/types.hh>
#include <numeric/kdtree/WrappedPrimitive.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>


namespace numeric {
namespace kdtree {

typedef WrappedPrimitive< char > WrappedChar;
	//typedef utility::pointer::owning_ptr< WrappedChar > WrappedCharOP;

} // kdtree
} // numeric

#endif
