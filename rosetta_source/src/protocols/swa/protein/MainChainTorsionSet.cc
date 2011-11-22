// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MainChainTorsionClasses
/// @brief a few functions used by several StepWiseProteinAnsatz classes
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/swa/protein/MainChainTorsionSet.hh>

//////////////////////////////////
#include <core/types.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
// AUTO-REMOVED #include <core/id/NamedAtomID.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>

// AUTO-REMOVED #include <core/pose/util.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.tmpl.hh>

// AUTO-REMOVED #include <numeric/angle.functions.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>
// AUTO-REMOVED #include <numeric/conversions.hh>

#include <string>

//Auto Headers
#include <utility/vector1.hh>


using core::Real;
using core::Size;
using core::pose::Pose;

namespace protocols {
namespace swa {
namespace protein {

	MainChainTorsionSet::MainChainTorsionSet( core::Real const & phi, core::Real const & psi, core::Real const & omega ):
		phi_( phi ),
		psi_( psi ),
		omega_( omega )
	{}

	MainChainTorsionSet::MainChainTorsionSet( core::Real const & phi, core::Real const & psi ):
		phi_( phi ),
		psi_( psi ),
		omega_( 180.0 )
	{}

	MainChainTorsionSet::~MainChainTorsionSet(){}

	core::Real const MainChainTorsionSet::phi() const{ return phi_; }
	core::Real const MainChainTorsionSet::psi() const{ return psi_; }
	core::Real const MainChainTorsionSet::omega() const{ return omega_; }

	MainChainTorsionSet
	MainChainTorsionSet::operator=( MainChainTorsionSet const & src )
	{
		phi_ =  src.phi();
		psi_ =  src.psi();
		omega_ =  src.omega();
		return (*this );
	}

}
}
}
