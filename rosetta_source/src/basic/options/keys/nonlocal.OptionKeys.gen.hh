// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   basic/options/keys/nonlocal.OptionKeys.gen.hh
/// @brief  basic::options::OptionKeys collection
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author James M. Thompson (tex@u.washington.edu)

#ifndef INCLUDED_basic_options_keys_nonlocal_OptionKeys_gen_HH
#define INCLUDED_basic_options_keys_nonlocal_OptionKeys_gen_HH

// Unit headers
#include <basic/options/keys/OptionKeys.hh>

namespace basic {
namespace options {
namespace OptionKeys {

namespace nonlocal { extern BooleanOptionKey const nonlocal; }
namespace nonlocal { extern FileOptionKey const moves; }
namespace nonlocal { extern StringOptionKey const mode; }
namespace nonlocal { extern StringOptionKey const builder; }
namespace nonlocal { extern BooleanOptionKey const randomize_missing; }
namespace nonlocal { extern IntegerOptionKey const gap_sampling_extension; }
namespace nonlocal { extern IntegerOptionKey const min_chunk_size; }
namespace nonlocal { extern IntegerOptionKey const max_chunk_size; }
namespace nonlocal { extern RealOptionKey const rdc_weight; }
namespace nonlocal { extern IntegerOptionKey const ramp_constraints_cycles; }
namespace nonlocal { extern BooleanOptionKey const ramp_chainbreaks; }
namespace nonlocal { extern BooleanOptionKey const ramp_constraints; }

} // namespace OptionKeys
} // namespace options
} // namespace basic

#endif
