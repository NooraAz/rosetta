// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:


#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <devel/init.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>

using namespace core::pose;
using namespace core::kinematics;
using namespace core::id;
using namespace utility;
using namespace numeric;


typedef numeric::xyzVector<core::Real> Vec;
typedef numeric::xyzMatrix<core::Real> Mat;

int main (int argc, char *argv[]) {

	try {

	devel::init(argc,argv);
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	utility::vector1<std::string> files = option[in::file::s]();
	for(int ifile = 1; ifile <= (int)files.size(); ++ifile) {
		core::pose::Pose pose;
		core::import_pose::pose_from_pdb(pose,files[ifile]);
		Size nres = pose.n_residue();
		Pose init(pose);
		core::pose::symmetry::make_symmetric_pose(pose);

	}
	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}



