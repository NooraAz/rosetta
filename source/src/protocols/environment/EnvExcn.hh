// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/EnvExcn.hh
/// @brief A class that is used to express a mover-specific DoF-unlock on a ProtectedPose. Its destruction expresses a re-locking.
///
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_EnvExcn_hh
#define INCLUDED_protocols_environment_EnvExcn_hh

// Unit Headers

// Package headers
#include <core/environment/EnvExcnBase.hh>
#include <core/environment/EnvCore.hh>
#include <protocols/environment/Environment.hh>
#include <protocols/environment/ProtectedConformation.fwd.hh>

// Project headers

// C++ Headers
#include <string>

// ObjexxFCL Headers

namespace protocols {
namespace environment {

class EXCN_Env_Security_Exception : public core::environment::EXCN_Env_Exception {
	typedef core::environment::EXCN_Env_Exception Parent;
public:
	EXCN_Env_Security_Exception( std::string const& dofname,
                               ProtectedConformation const&,
                               std::string const& mover_name,
                               EnvironmentCAP env ) {
		std::ostringstream msg;
		msg << "ProtectedConformation reported illegal access to '" << dofname
        << "' by mover '" << mover_name << "'.";

    Parent( msg.str(), env );
	}
};

class EXCN_Env_Passport : public core::environment::EXCN_Env_Exception {
  typedef core::environment::EXCN_Env_Exception Parent;
public:
  EXCN_Env_Passport( std::string const& message,
                     std::string const& mover_name,
                     EnvironmentCAP env ) {
    std::ostringstream msg;
    msg << message << " Mover: '" << mover_name << "'.";

    Parent( msg.str(), env );
  }
};

} // environment
} // protocols

#endif //INCLUDED_protocols_environment_EnvExcn_hh
