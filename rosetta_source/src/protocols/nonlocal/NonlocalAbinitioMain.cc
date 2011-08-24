// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/NonlocalAbinitioMain.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Project headers
#include <protocols/nonlocal/NonlocalAbinitio.hh>
#include <protocols/nonlocal/NonlocalAbinitioMain.hh>
#include <protocols/nonlocal/util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>

// Utility headers
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/nonlocal.OptionKeys.gen.hh>
#include <utility/exit.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/vector1.hh>

// External headers
#include <boost/algorithm/string/case_conv.hpp>

// C/C++ headers
#include <iostream>
#include <string>

namespace protocols  {
namespace nonlocal {

static basic::Tracer TR("protocols.nonlocal.NonlocalAbinitioMain");

using namespace std;

void* NonlocalAbinitio_main(void*) {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using protocols::jd2::JobDistributor;
  using protocols::moves::MoverOP;
  using utility::vector1;
  using utility::excn::EXCN_Base;

  // optional configuration
  NonlocalAbinitio::KinematicPolicy policy = get_policy_or_die();

  MoverOP mover;
  if (option[in::file::alignment].user() &&
      option[in::file::template_pdb].user() &&
      option[cm::aln_format].user()) {
    // randomize the coordinates of missing loops to avoid clashes
    if (!option[OptionKeys::nonlocal::randomize_missing].user()) {
      option[OptionKeys::nonlocal::randomize_missing](true);
      TR.Warning << "Enabling -nonlocal:randomize_missing" << std::endl;
    }

    vector1<NLGrouping> groupings;
    protocols::nonlocal::nonlocal_groupings_from_alignment(&groupings);
    mover = new NonlocalAbinitio(groupings, policy);
  } else {
    // read groupings from options
    mover = new NonlocalAbinitio(policy);
  }

  // hand the configured mover to the job distributor for execution
  try {
    JobDistributor::get_instance()->go(mover);
  } catch (EXCN_Base& e) {
    cerr << "Exception: " << endl;
    e.show(cerr);
  }

  // prevent "control may reach end of non-void function" warning
  return 0;
}

NonlocalAbinitio::KinematicPolicy get_policy_or_die() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  std::string mode = option[OptionKeys::nonlocal::mode]();
  boost::to_lower(mode);
  if (!(mode == "rigid" || mode == "semirigid"))
    utility_exit_with_message("Invalid setting for option nonlocal:mode");

  return (mode == "rigid")
      ? NonlocalAbinitio::RIGID : NonlocalAbinitio::SEMI_RIGID;
}

}  // namespace nonlocal
}  // namespace protocols
