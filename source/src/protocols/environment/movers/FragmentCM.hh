// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/FragmentCM.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_FragmentCM_hh
#define INCLUDED_protocols_environment_FragmentCM_hh

// Unit Headers
#include <protocols/environment/movers/FragmentCM.fwd.hh>

// Package headers
#include <protocols/environment/ClaimingMover.hh>
#include <protocols/simple_moves/FragmentMover.fwd.hh>

// Project headers

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace environment {

class FragmentCM : public protocols::environment::ClaimingMover {

public:
  FragmentCM();

  FragmentCM( simple_moves::FragmentMoverOP, std::string const& );

  virtual void set_label( std::string const& label );

  virtual void set_mover( simple_moves::FragmentMoverOP mover );

  virtual ~FragmentCM();

  virtual claims::EnvClaims yield_claims( core::pose::Pose& );

  virtual void initialize( Pose& pose );

  virtual void apply( Pose& pose );

  virtual std::string get_name() const;

  std::string const& label() const { return label_; }

protected:
  virtual void passport_updated();

  simple_moves::FragmentMoverOP mover() const { return mover_; };

private:
  simple_moves::FragmentMoverOP mover_;
  std::string label_;

}; // end FragmentCM base class

} // environment
} // protocols

#endif //INCLUDED_protocols_environment_FragmentCM_hh
