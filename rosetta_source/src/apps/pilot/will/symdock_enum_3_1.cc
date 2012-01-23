#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmData.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/SymmetryInfo.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>
#include <core/id/AtomID.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/sasa.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>
// AUTO-REMOVED #include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/string_util.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <utility/vector1.hh>
static basic::Tracer TR("symdock_enum");

// #ifdef USE_OPENMP
// #pragma omp parallel for schedule(dynamic,1)
// #endif
// #ifdef USE_OPENMP
// #pragma omp critical
// #endif



OPT_1GRP_KEY( Real	  , tcdock, clash_dis		)
OPT_1GRP_KEY( Real	  , tcdock, contact_dis	)
OPT_1GRP_KEY( Boolean , tcdock, intra	)
void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::s );
	NEW_OPT( tcdock::clash_dis	 ,"max acceptable clash dis"	,	3.5 );
	NEW_OPT( tcdock::contact_dis ,"max acceptable contact dis", 10.0 );	
	NEW_OPT( tcdock::intra       ,"include intra 3-3 5-5 contacts", true );	
}


using core::Size;
using core::Real;
using core::pose::Pose;
using utility::vector1;
using ObjexxFCL::fmt::I;
using ObjexxFCL::fmt::F;
using numeric::min;
using numeric::max;

typedef numeric::xyzVector<core::Real> Vec;
typedef numeric::xyzMatrix<core::Real> Mat;
typedef numeric::xyzVector<double> Vecf;
typedef numeric::xyzMatrix<double> Matf;

inline double sqr(double x) { return x*x; }
inline float	sqr( float x) { return x*x; }

inline Real sigmoid( Real const & sqdist, Real const & start, Real const & stop ) {
	if( sqdist > stop*stop ) {
		return 0.0;
	} else if( sqdist < start*start ) {
		return 1.0;
	} else {
		Real dist = sqrt( sqdist );
		return sqr(1.0	- sqr( (dist - start) / (stop - start) ) );
	}
}

void trans_pose( Pose & pose, Vec const & trans ) {
	for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, pose.xyz(aid) + trans );
		}
	}
}

void rot_pose( Pose & pose, Mat const & rot ) {
	for(Size ir = 1; ir <= pose.n_residue(); ++ir) {
		for(Size ia = 1; ia <= pose.residue_type(ir).natoms(); ++ia) {
			core::id::AtomID const aid(core::id::AtomID(ia,ir));
			pose.set_xyz( aid, rot * pose.xyz(aid) );
		}
	}
}

void rot_pose( Pose & pose, Mat const & rot, Vec const & cen ) {
	trans_pose(pose,-cen);
	rot_pose(pose,rot);
	trans_pose(pose,cen);
}

void rot_pose( Pose & pose, Vec const & axis, Real const & ang ) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang));
}

void rot_pose( Pose & pose, Vec const & axis, Real const & ang, Vec const & cen ) {
	rot_pose(pose,rotation_matrix_degrees(axis,ang),cen);
}


void alignaxis(Pose & pose, Vec newaxis, Vec oldaxis, Vec cen = Vec(0,0,0) ) {
	newaxis.normalize();
	oldaxis.normalize();
	Vec axis = newaxis.cross(oldaxis).normalized();
	Real ang = -acos(numeric::max(-1.0,numeric::min(1.0,newaxis.dot(oldaxis))))*180/numeric::constants::d::pi;
	rot_pose(pose,axis,ang,cen);
}

int cbcount_vec(vector1<Vecf> & cba, vector1<Vecf> & cbb) {
	int cbcount = 0;
	for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia)
		for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib)
			if( ib->distance_squared(*ia) < 100.0 ) cbcount++;
	return cbcount;
}

void prune_cb_pairs_dis10(vector1<Vecf> & cba, vector1<Vecf> & cbb) {
	vector1<Vecf> a,b;
	for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
		for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
			if( ib->distance_squared(*ia) < 100.0 ) {
				a.push_back(*ia);
				b.push_back(*ib);
			}
		}
	}
	cba = a;
	cbb = b;
}

int pose_cbcount(Pose const & a, Pose const & b) {
	int count = 0;
	for(Size i = 1; i <= a.n_residue(); ++i) {
		for(Size j = 1; j <= b.n_residue(); ++j) {
			if(a.residue(i).xyz(2).distance_squared(b.residue(j).xyz(2)) < 100.0) {
				count++;
			}
		}
	}
	return count;
}



double
sicsafe(
				vector1<Vecf> const & pa,
				vector1<Vecf> const & pb,
				vector1<Vecf> const & cba,
				vector1<Vecf> const & cbb,
				Vecf ori,
				int & cbcount,
				bool debug = false
){
	Real const CONTACT_D	= basic::options::option[basic::options::OptionKeys::tcdock::contact_dis]();
	Real const CLASH_D		= basic::options::option[basic::options::OptionKeys::tcdock::	clash_dis]();
	Real const CONTACT_D2 = sqr(CONTACT_D);
	Real const CLASH_D2	 = sqr(CLASH_D);
	ori.normalize();
	Vec axis = Vecf(0,0,1).cross(ori);
	Matf R = rotation_matrix( axis, -acos(numeric::max(-1.0,numeric::min(1.0,Vecf(0,0,1).dot(ori)))) );
	Real mindis = 9e9;
	int mni,mnj;
	// for(vector1<Vec>::const_iterator i = pa.begin(); i != pa.end(); ++i) {
	// for(vector1<Vec>::const_iterator j = pb.begin(); j != pb.end(); ++j) {
	for(int ii = 1; ii <= pa.size(); ++ii) {
		Vec const i = R * pa[ii];
		for(int jj = 1; jj <= pb.size(); ++jj) {		
			Vec const j = R * pb[jj];
			Real const dxy2 = (i.x()-j.x())*(i.x()-j.x()) + (i.y()-j.y())*(i.y()-j.y());
			if( dxy2 >= CLASH_D2 ) continue;
			Real const dz = j.z() - i.z() - sqrt(CLASH_D2-dxy2);
			if( dz < mindis) {
				mindis = dz;
				mni = ii;
				mnj = jj;
			}
		}		
	}
	
	cbcount = 0.0;
	for(vector1<Vec>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
		Vec const & va = *ia;
		for(vector1<Vec>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
			Vec const & vb = *ib;
			Real d2 = vb.distance_squared( (va + mindis*ori ) );
			if( d2 < CONTACT_D2 ) {
				cbcount += 1.0;//sigmoid(d2, CLASH_D, CONTACT_D );
			}
		}
	}
	return mindis;
	
}

double
sicfast(
				vector1<Vecf>	pa,
				vector1<Vecf>	pb,
				vector1<Vecf> const & cba,
				vector1<Vecf> const & cbb,
				Vecf ori,
				int & cbcount,
				bool debug = false
){
	double BIN = 2.0;

	// get points, rotated ro ori is 0,0,1, might already be done
	Matf rot = Matf::identity();
	if		 ( ori.dot(Vec(0,0,1)) < -0.99999 ) rot = rotation_matrix( Vec(1,0,0).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
	else if( ori.dot(Vec(0,0,1)) <	0.99999 ) rot = rotation_matrix( Vec(0,0,1).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
	if( rot != Matf::identity() ) {
		for(vector1<Vecf>::iterator ia = pa.begin(); ia != pa.end(); ++ia) *ia = rot*(*ia);
		for(vector1<Vecf>::iterator ib = pb.begin(); ib != pb.end(); ++ib) *ib = rot*(*ib);
	}

	// get bounds for plane hashes
	double xmx1=-9e9,xmn1=9e9,ymx1=-9e9,ymn1=9e9,xmx=-9e9,xmn=9e9,ymx=-9e9,ymn=9e9;
	for(vector1<Vecf>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia) {
		xmx1 = max(xmx1,ia->x()); xmn1 = min(xmn1,ia->x());
		ymx1 = max(ymx1,ia->y()); ymn1 = min(ymn1,ia->y());
	}
	for(vector1<Vecf>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib) {
		xmx = max(xmx,ib->x()); xmn = min(xmn,ib->x());
		ymx = max(ymx,ib->y()); ymn = min(ymn,ib->y());
	}
	xmx = min(xmx,xmx1); xmn = max(xmn,xmn1);
	ymx = min(ymx,ymx1); ymn = max(ymn,ymn1);


	int xlb = (int)floor(xmn/BIN)-2; int xub = (int)ceil(xmx/BIN)+2; // one extra on each side for correctness,
	int ylb = (int)floor(ymn/BIN)-2; int yub = (int)ceil(ymx/BIN)+2; // and one extra for outside atoms

	// TR << "BOUNDS " << xmn << " " << xmx << " " << ymn << " " << ymx << std::endl;
	// TR << "BOUNDS " << xlb << " " << xub << " " << ylb << " " << yub << std::endl;

	// insert points into hashes
	int const xsize = xub-xlb+1;
	int const ysize = yub-ylb+1;
	ObjexxFCL::FArray2D<Vecf> ha(xsize,ysize,Vecf(0,0,-9e9)),hb(xsize,ysize,Vecf(0,0,9e9));
	for(vector1<Vecf>::const_iterator ia = pa.begin(); ia != pa.end(); ++ia) {
		// int const ix = min(xsize,max(1,(int)ceil(ia->x()/BIN)-xlb));
		// int const iy = min(ysize,max(1,(int)ceil(ia->y()/BIN)-ylb));
		int const ix = (int)ceil(ia->x()/BIN)-xlb;
		int const iy = (int)ceil(ia->y()/BIN)-ylb;
		if( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) continue;
		if( ha(ix,iy).z() < ia->z() ) ha(ix,iy) = *ia;
	}
	for(vector1<Vecf>::const_iterator ib = pb.begin(); ib != pb.end(); ++ib) {
		// int const ix = min(xsize,max(1,(int)ceil(ib->x()/BIN)-xlb));
		// int const iy = min(ysize,max(1,(int)ceil(ib->y()/BIN)-ylb));
		int const ix = (int)ceil(ib->x()/BIN)-xlb;
		int const iy = (int)ceil(ib->y()/BIN)-ylb;
		if( ix < 1 || ix > xsize || iy < 1 || iy > ysize ) continue;
		if( hb(ix,iy).z() > ib->z() ) hb(ix,iy) = *ib;
	}

	// check hashes for min dis
	int imna=0,jmna=0,imnb=0,jmnb=0;
	double mindis = 9e9;
	for(int i = 1; i <= xsize; ++i) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
		for(int j = 1; j <= ysize; ++j) {
			for(int k = -2; k <= 2; ++k) {
				if(i+k < 1 || i+k > xsize) continue;
				for(int l = -2; l <= 2; ++l) {
					if(j+l < 1 || j+l > ysize) continue;
					double const xa = ha(i	,j	).x();
					double const ya = ha(i	,j	).y();
					double const xb = hb(i+k,j+l).x();
					double const yb = hb(i+k,j+l).y();
					double const d2 = (xa-xb)*(xa-xb) + (ya-yb)*(ya-yb);

					if( d2 < 16.0 ) {
						double dz = hb(i+k,j+l).z() - ha(i,j).z() - sqrt(16.0-d2);
						if( dz < mindis ) {
							mindis = dz;
							imna = i;
							jmna = j;
							imnb = i+k;
							jmnb = j+l;
						}
					}
				}
			}
		}
	}

	// {
	//	utility::io::ozstream out("cba.pdb");
	//	for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
	//		Vec viz = (*ia) + (mindis*ori);
	//		out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"A"<<I(4,100)<<"		"<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	//	}
	//	out.close();
	// }
	// {
	//	utility::io::ozstream out("cbb.pdb");
	//	for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
	//		Vec viz = (*ib);
	//		out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"B"<<I(4,100)<<"		"<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
	//	}
	//	out.close();
	// }

	cbcount = 0;
	// utility::io::ozstream out("cb8.pdb");
	// TR << "CB0 " << cbcount << std::endl;
	for(vector1<Vecf>::const_iterator ia = cba.begin(); ia != cba.end(); ++ia) {
		for(vector1<Vecf>::const_iterator ib = cbb.begin(); ib != cbb.end(); ++ib) {
			if( ib->distance_squared( (*ia) + (mindis*ori) ) < 100.0 ) {
				cbcount++;
				// Vec viz = (*ia) + (mindis*ori);
				// out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' <<	"VIZ"<<' '<<"A"<<I(4,100)<<"		"<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
				// viz = *ib;
				// out<<"HETATM"<<I(5,1000)<<' '<<"VIZ "<<' ' <<	"VIZ"<<' '<<"B"<<I(4,100)<<"		"<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
			}
		}
	}
	// out.close();
	// TR << "CB1 " << cbcount << std::endl;

	// // rotate points back -- needed iff pa/pb come by reference
	rot = rot.transposed();
	// if( rot != Matf::identity() ) {
	//	for(vector1<Vecf>::iterator ia = pa.begin(); ia != pa.end(); ++ia) *ia = rot*(*ia);
	//	for(vector1<Vecf>::iterator ib = pb.begin(); ib != pb.end(); ++ib) *ib = rot*(*ib);
	// }

	// uncomment this to get hashes in local space
	// rot = Matf::identity();
	// ori = Vec(0,0,1);

	if(debug){
		{
			utility::io::ozstream out("hasha.pdb");
			for(int i = 2; i <= xsize-1; ++i) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
				for(int j = 2; j <= ysize-1; ++j) {
					Vecf viz = rot*ha(i,j) + mindis*ori;
					if(viz.z() < -9e8 || 9e8 < viz.z()) continue;
					out<<"HETATM"<<I(5,1000+i)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"B"<<I(4,100+j)<<"		"<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
				}
			}
			Vecf viz = rot*ha(imna,jmna) + mindis*ori;
			out<<"HETATM"<<I(5,1000+imna)<<' '<<"MIN "<<' ' <<	"MIN"<<' '<<"B"<<I(4,100+jmna)<<"		"<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
			out.close();
		}
		{
			utility::io::ozstream out("hashb.pdb");
			for(int i = 2; i <= xsize-1; ++i) { // skip 1 and N because they contain outside atoms (faster than clashcheck?)
				for(int j = 2; j <= ysize-1; ++j) {
					Vecf viz = rot*hb(i,j);
					if(viz.z() < -9e8 || 9e8 < viz.z()) continue;
					out<<"HETATM"<<I(5,1000+i)<<' '<<"VIZ "<<' ' << "VIZ"<<' '<<"C"<<I(4,100+j)<<"		"<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
				}
			}
			Vecf viz = rot*hb(imnb,jmnb);
			out<<"HETATM"<<I(5,1000+imnb)<<' '<<"MIN "<<' ' <<	"MIN"<<' '<<"C"<<I(4,100+jmnb)<<"		"<<F(8,3,viz.x())<<F(8,3,viz.y())<<F(8,3,viz.z())<<F(6,2,1.0)<<F(6,2,1.0)<<'\n';
			out.close();
		}
	}

	return mindis;
}

double sicfast(
							 Pose const & a,
							 Pose const & b,
							 Vecf ori_in,
							 int & cbcount
){
	// get points, rotated ro ori is 0,0,1
	vector1<Vecf> pa,pb;
	vector1<Vecf> cba,cbb;
	Vecf ori = ori_in.normalized();
	Matf rot = Matf::identity();
	if		 ( ori.dot(Vec(0,0,1)) < -0.999 ) rot = rotation_matrix( Vec(1,0,0).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
	else if( ori.dot(Vec(0,0,1)) <	0.999 ) rot = rotation_matrix( Vec(0,0,1).cross(ori), -acos(Vec(0,0,1).dot(ori)) );
	for(int i = 1; i <= (int)a.n_residue(); ++i) {
		cba.push_back(rot*Vecf(a.residue(i).xyz(2)));
		int const natom = (a.residue(i).name3()=="GLY") ? 4 : 5;
		for(int j = 1; j <= natom; ++j) pa.push_back(rot*Vecf(a.residue(i).xyz(j)));
	}
	for(int i = 1; i <= (int)b.n_residue(); ++i) {
		cbb.push_back(rot*Vecf(b.residue(i).xyz(2)));
		int const natom = (b.residue(i).name3()=="GLY") ? 4 : 5;
		for(int j = 1; j <= natom; ++j) pb.push_back(rot*Vecf(b.residue(i).xyz(j)));
	}
	return sicfast( pa, pb, cba, cbb, Vec(0,0,1), cbcount );
}



void make_dimer(core::pose::Pose & pose) {
	core::pose::Pose t2(pose);
	rot_pose(t2,Vec(0,0,1),180.0);
	for(Size i = 1; i <= t2.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t2.residue(i),1); else pose.append_residue_by_bond(t2.residue(i));
}

void make_trimer(core::pose::Pose & pose) {
	core::pose::Pose t2(pose),t3(pose);
	rot_pose(t2,Vec(0,0,1),120.0);
	rot_pose(t3,Vec(0,0,1),240.0);
	for(Size i = 1; i <= t2.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t2.residue(i),1); else pose.append_residue_by_bond(t2.residue(i));
	for(Size i = 1; i <= t3.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t3.residue(i),1); else pose.append_residue_by_bond(t3.residue(i));
}

void make_pentamer(core::pose::Pose & pose) {
	core::pose::Pose t2(pose),t3(pose),t4(pose),t5(pose);
	rot_pose(t2,Vec(0,0,1), 72.0);
	rot_pose(t3,Vec(0,0,1),144.0);
	rot_pose(t4,Vec(0,0,1),216.0);
	rot_pose(t5,Vec(0,0,1),288.0);
	for(Size i = 1; i <= t2.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t2.residue(i),1); else pose.append_residue_by_bond(t2.residue(i));
	for(Size i = 1; i <= t3.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t3.residue(i),1); else pose.append_residue_by_bond(t3.residue(i));
	for(Size i = 1; i <= t4.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t4.residue(i),1); else pose.append_residue_by_bond(t4.residue(i));
	for(Size i = 1; i <= t5.n_residue(); ++i) if(pose.residue(i).is_lower_terminus()) pose.append_residue_by_jump(t5.residue(i),1); else pose.append_residue_by_bond(t5.residue(i));
}

void tri1_to_tri2(core::pose::Pose & pose, Vecf triaxs, Vecf triax2) {
	Matf r = rotation_matrix_degrees( triaxs.cross(triax2), angle_degrees(triaxs,Vec(0,0,0),triax2) );
	r = r * rotation_matrix_degrees(triaxs,60.0);
	rot_pose(pose,r);
}
void tri1_to_tri2(vector1<Vecf> & pts, Vecf triaxs, Vecf triax2) {
	Matf r = rotation_matrix_degrees( triaxs.cross(triax2), angle_degrees(triaxs,Vec(0,0,0),triax2) );
	r = r * rotation_matrix_degrees(triaxs,60.0);
	for(vector1<Vecf>::iterator i = pts.begin(); i != pts.end(); ++i) *i = r*(*i);
}

void pnt1_to_pnt2(core::pose::Pose & pose, Vecf pntaxs, Vecf pntax2) {
	Matf r = rotation_matrix_degrees( pntaxs.cross(pntax2), angle_degrees(pntaxs,Vec(0,0,0),pntax2) );
	r = r * rotation_matrix_degrees(pntaxs,36.0);
	rot_pose(pose,r);
}
void pnt1_to_pnt2(vector1<Vecf> & pts, Vecf pntaxs, Vecf pntax2) {
	Matf r = rotation_matrix_degrees( pntaxs.cross(pntax2), angle_degrees(pntaxs,Vec(0,0,0),pntax2) );
	r = r * rotation_matrix_degrees(pntaxs,36.0);
	for(vector1<Vecf>::iterator i = pts.begin(); i != pts.end(); ++i) *i = r*(*i);
}

void run( Size itrifile, Size ipntfile ) {
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using namespace core::id;
	using numeric::conversions::radians;

	core::chemical::ResidueTypeSetCAP crs=core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::CENTROID);

	Pose tri_in,pnt_in;
	core::import_pose::pose_from_pdb(tri_in,*crs,option[in::file::s]()[itrifile]);
	core::import_pose::pose_from_pdb(pnt_in,*crs,option[in::file::s]()[ipntfile]);
	std::string trifile = utility::file_basename(option[in::file::s]()[itrifile]);
	//TR << "make trimer and pentamer" << std::endl;
	make_trimer(tri_in);
	//tri_in.dump_pdb("trimer.pdb");
	make_pentamer(pnt_in);
	//pnt_in.dump_pdb("pentamer.pdb");
	//std::exit(-1);
	//Real sasa_tri = core::scoring::calc_total_sasa( tri_in, 2.0 );
	//Real sasa_pnt = core::scoring::calc_total_sasa( pnt_in, 2.0 );	
	

	// set up geometry
	Vecf triaxs = Vec( 0.000000, 0.000000,1.000000).normalized();
	Vecf triax2 = Vec(-0.333333,-0.577350,0.745356).normalized(); // 33.4458470159, 10.42594
	Vecf pntaxs = Vec(-0.607226, 0.000000,0.794529).normalized();
	Vecf pntax2 = Vec(-0.491123,-0.850651,0.187593).normalized(); // 63.4311873349, 5.706642
	Real alpha = angle_degrees(triaxs,Vec(0,0,0),pntaxs);
	rot_pose(pnt_in,Vec(0,1,0),-alpha,Vec(0,0,0));

	ObjexxFCL::FArray3D<float> gdsp(72,120,360,-1);
	ObjexxFCL::FArray3D<float> gcbc(72,120,360,-1);	

	// Vecf tcom(0,0,0),pcom(0,0,0);
	// for(Size i = 1; i <= tri_in.n_residue()/3; ++i) for(Size j = 1; j <= 5; ++j) tcom += tri_in.residue(i).xyz(j);
	// tcom /= double(5*tri_in.n_residue()/3);
	// for(Size i = 1; i <= pnt_in.n_residue()/5; ++i) for(Size j = 1; j <= 5; ++j) pcom += pnt_in.residue(i).xyz(j);
	// pcom /= double(5*tri_in.n_residue()/5);
	// rot_pose(tri_in,triaxs,dihedral_degrees(pntaxs,Vec(0,0,0),triaxs,tcom));
	// rot_pose(pnt_in,pntaxs,dihedral_degrees(triaxs,Vec(0,0,0),pntaxs,pcom));
	// tri_in.dump_pdb("t0.pdb");
	// pnt_in.dump_pdb("p0.pdb");
	// utility_exit_with_message("debug");

	Pose const triinit(tri_in);
	Pose const pntinit(pnt_in);

	int ANGLE_INCR = 3;
	if(ANGLE_INCR != 3) utility_exit_with_message("first ANGLE_INCR must be 3!!!");
	ObjexxFCL::FArray3D<int>	 cbcount3((Size)floor(72.0/ANGLE_INCR),(Size)floor(120.0/ANGLE_INCR),(Size)floor(360.0/ANGLE_INCR),0);
	// 3 deg. sampling
	{
		// compute high/low min dis for pent and tri here, input to sicfast and don't allow any below
		TR << "precomputing pent-pent and tri-tri interactions every 3°" << std::endl;
		vector1<double>						pntmnpos(72/ANGLE_INCR,0),	 trimnpos(120/ANGLE_INCR,0);
		vector1<double>						pntmnneg(72/ANGLE_INCR,0),	 trimnneg(120/ANGLE_INCR,0);
		ObjexxFCL::FArray2D<int>	pntcbpos(72/ANGLE_INCR,97,0),tricbpos(120/ANGLE_INCR,145,0);
		ObjexxFCL::FArray2D<int>	pntcbneg(72/ANGLE_INCR,97,0),tricbneg(120/ANGLE_INCR,145,0);
		{
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
			for(int ipnt = 0; ipnt < 72; ipnt+=ANGLE_INCR) {
				Pose p;
#ifdef USE_OPENMP
#pragma omp critical
#endif
				p = pntinit;
				rot_pose(p,pntaxs,ipnt);
				vector1<Vecf> ppnt,cbp; // precompute these
				for(int i = 1; i <= (int)p.n_residue(); ++i) {
					cbp.push_back(Vecf(p.residue(i).xyz(2)));
					int const natom = (p.residue(i).name3()=="GLY") ? 4 : 5;
					for(int j = 1; j <= natom; ++j) ppnt.push_back(Vecf(p.residue(i).xyz(j)));
				}
				vector1<Vecf> ppn2(ppnt),cb2(cbp);
				pnt1_to_pnt2(ppn2,pntaxs,pntax2);
				pnt1_to_pnt2( cb2,pntaxs,pntax2);				
				int cbcount = 0;
				double const d = sicsafe(ppnt,ppn2,cbp,cb2,(pntax2-pntaxs).normalized(),cbcount,false);
				if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for pntpos! "+ObjexxFCL::string_of(ipnt));
				for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(pntax2-pntaxs).normalized();
				pntmnpos[ipnt/ANGLE_INCR+1] = -d/2.0/sin( angle_radians(pntax2,Vec(0,0,0),pntaxs)/2.0 );
				pntcbpos(ipnt/ANGLE_INCR+1,1) = cbcount;
				// std::cerr << "compute CB" << std::endl;
				prune_cb_pairs_dis10(cbp,cb2);
				for(int i = 2; i <= 97; ++i) {
					for(vector1<Vecf>::iterator iv = cbp.begin(); iv != cbp.end(); ++iv) *iv = (*iv) + 0.1*pntaxs;
					for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) + 0.1*pntax2;
					int cbc = 0; for(Size j = 1; j <= cbp.size(); ++j) if(cbp[j].distance_squared(cb2[j]) < 100) cbc++;
					pntcbpos(ipnt/ANGLE_INCR+1,i) = cbc;
					if(cbc==0) break;
				}
			}
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
			for(int ipnt = 0; ipnt < 72; ipnt+=ANGLE_INCR) {
				Pose p;
#ifdef USE_OPENMP
#pragma omp critical
#endif
				p = pntinit;
				rot_pose(p,pntaxs,ipnt);
				vector1<Vecf> ppnt,cbp; // precompute these
				for(int i = 1; i <= (int)p.n_residue(); ++i) {
					cbp.push_back(Vecf(p.residue(i).xyz(2)));
					int const natom = (p.residue(i).name3()=="GLY") ? 4 : 5;
					for(int j = 1; j <= natom; ++j) ppnt.push_back(Vecf(p.residue(i).xyz(j)));
				}
				vector1<Vecf> ppn2(ppnt),cb2(cbp);
				pnt1_to_pnt2(ppn2,pntaxs,pntax2);
				pnt1_to_pnt2( cb2,pntaxs,pntax2);
				int cbcount = 0;
				double const d = sicsafe(ppnt,ppn2,cbp,cb2,(pntaxs-pntax2).normalized(),cbcount,false);				
				if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for pntneg! "+ObjexxFCL::string_of(ipnt));
				for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(pntaxs-pntax2).normalized();
				pntmnneg[ipnt/ANGLE_INCR+1] = d/2.0/sin( angle_radians(pntax2,Vec(0,0,0),pntaxs)/2.0 );
				pntcbneg(ipnt/ANGLE_INCR+1,1) = cbcount;
				// std::cerr << "compute CB" << std::endl;
				prune_cb_pairs_dis10(cbp,cb2);
				for(int i = 2; i <= 97; ++i) {
					for(vector1<Vecf>::iterator iv = cbp.begin(); iv != cbp.end(); ++iv) *iv = (*iv) - 0.1*pntaxs;
					for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) - 0.1*pntax2;
					int cbc = 0; for(Size j = 1; j <= cbp.size(); ++j) if(cbp[j].distance_squared(cb2[j]) < 100) cbc++;
					pntcbneg(ipnt/ANGLE_INCR+1,i) = cbc;
					if(cbc==0) break;
				}
			}
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
			for(int itri = 0; itri < 120; itri+=ANGLE_INCR) {
				Pose t;
#ifdef USE_OPENMP
#pragma omp critical
#endif
				t = triinit;
				rot_pose(t,triaxs,itri);
				vector1<Vecf> ptri,cbt; // precompute these
				for(int i = 1; i <= (int)t.n_residue(); ++i) {
					cbt.push_back(Vecf(t.residue(i).xyz(2)));
					int const natom = (t.residue(i).name3()=="GLY") ? 4 : 5;
					for(int j = 1; j <= natom; ++j) ptri.push_back(Vecf(t.residue(i).xyz(j)));
				}
				vector1<Vecf> ptr2(ptri),cb2(cbt);
				tri1_to_tri2(ptr2,triaxs,triax2);
				tri1_to_tri2( cb2,triaxs,triax2);				
				int cbcount = 0;
				double const d = sicsafe(ptri,ptr2,cbt,cb2,(triax2-triaxs).normalized(),cbcount,false);
				if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for tripos! "+ObjexxFCL::string_of(itri));
				for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(triax2-triaxs).normalized();
				trimnpos[itri/ANGLE_INCR+1] = -d/2.0/sin( angle_radians(triax2,Vec(0,0,0),triaxs)/2.0 );
				tricbpos(itri/ANGLE_INCR+1,1) = cbcount;
				// std::cerr << "compute CB" << std::endl;
				prune_cb_pairs_dis10(cbt,cb2);
				for(int i = 2; i <= 145; ++i) {
					for(vector1<Vecf>::iterator iv = cbt.begin(); iv != cbt.end(); ++iv) *iv = (*iv) + 0.1*triaxs;
					for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) + 0.1*triax2;
					int cbc = 0; for(Size j = 1; j <= cbt.size(); ++j) if(cbt[j].distance_squared(cb2[j]) < 100) cbc++;
					tricbpos(itri/ANGLE_INCR+1,i) = cbc;
					if(cbc==0) break;
				}
			}
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
			for(int itri = 0; itri < 120; itri+=ANGLE_INCR) {
				Pose t;
#ifdef USE_OPENMP
#pragma omp critical
#endif
				t = triinit;
				rot_pose(t,triaxs,itri);
				vector1<Vecf> ptri,cbt; // precompute these
				for(int i = 1; i <= (int)t.n_residue(); ++i) {
					cbt.push_back(Vecf(t.residue(i).xyz(2)));
					int const natom = (t.residue(i).name3()=="GLY") ? 4 : 5;
					for(int j = 1; j <= natom; ++j) ptri.push_back(Vecf(t.residue(i).xyz(j)));
				}
				vector1<Vecf> ptr2(ptri),cb2(cbt);
				tri1_to_tri2(ptr2,triaxs,triax2);
				tri1_to_tri2( cb2,triaxs,triax2);				
				int cbcount = 0;
				double const d = sicsafe(ptri,ptr2,cbt,cb2,(triaxs-triax2).normalized(),cbcount,false);
				if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for trineg! "+ObjexxFCL::string_of(itri));
				for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(triaxs-triax2).normalized();
				trimnneg[itri/ANGLE_INCR+1] = d/2.0/sin( angle_radians(triax2,Vec(0,0,0),triaxs)/2.0 );
				tricbneg(itri/ANGLE_INCR+1,1) = cbcount;
				// std::cerr << "compute CB" << std::endl;
				prune_cb_pairs_dis10(cbt,cb2);
				for(int i = 2; i <= 145; ++i) {
					for(vector1<Vecf>::iterator iv = cbt.begin(); iv != cbt.end(); ++iv) *iv = (*iv) - 0.1*triaxs;
					for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) - 0.1*triax2;
					int cbc = 0; for(Size j = 1; j <= cbt.size(); ++j) if(cbt[j].distance_squared(cb2[j]) < 100) cbc++;
					tricbneg(itri/ANGLE_INCR+1,i) = cbc;
					if(cbc==0) break;
				}
			}
		}
		
		
		TR << "main loop 1 over ipnt, itri, iori every 3 degrees" << std::endl;
		ObjexxFCL::FArray3D<int>	 cbcount((Size)floor(72.0/ANGLE_INCR),(Size)floor(120.0/ANGLE_INCR),(Size)floor(360.0/ANGLE_INCR),0);
		ObjexxFCL::FArray3D<float> surfdis3((Size)floor(72.0/ANGLE_INCR),(Size)floor(120.0/ANGLE_INCR),(Size)floor(360.0/ANGLE_INCR),0.0);
		{
			double mxd = 0;
			int cbmax = 0, mxiori = 0, mxipnt = 0, mxitri = 0;
			for(int ipnt = 0; ipnt < 72; ipnt+=ANGLE_INCR) {
				if(ipnt%9==0 && ipnt!=0) TR << "	 loop1 " << trifile << " " << (100*ipnt)/72 << "\% done, cbmax:" << cbmax << std::endl;
				Pose p = pntinit;
				rot_pose(p,pntaxs,(Real)ipnt);
				vector1<Vecf> pb,cbb; // precompute these
				for(int i = 1; i <= (int)p.n_residue(); ++i) {
					cbb.push_back(Vecf(p.residue(i).xyz(2)));
					int const natom = (p.residue(i).name3()=="GLY") ? 4 : 5;
					for(int j = 1; j <= natom; ++j) pb.push_back(Vecf(p.residue(i).xyz(j)));
				}
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
				for(int itri = 0; itri < 120; itri+=ANGLE_INCR) {
					Pose t;
#ifdef USE_OPENMP
#pragma omp critical
#endif
					t = triinit;
					rot_pose(t,triaxs,(Real)itri);
					vector1<Vecf> pa,cba; // precompute these
					for(int i = 1; i <= (int)t.n_residue(); ++i) {
						cba.push_back(Vecf(t.residue(i).xyz(2)));
						int const natom = (t.residue(i).name3()=="GLY") ? 4 : 5;
						for(int j = 1; j <= natom; ++j) pa.push_back(Vecf(t.residue(i).xyz(j)));
					}

					int iori = -1, ori_stage = 1;
					bool newstage = true;
					// for(iori = 0; iori < 360; iori+=ANGLE_INCR)
					while(ori_stage < 5) {
						if(newstage) {
							if( ori_stage == 1 || ori_stage == 2 ) iori = (int)( 90.0+double(ANGLE_INCR)/2.0+angle_degrees(triaxs,Vecf(0,0,0),pntaxs));
							if( ori_stage == 3 || ori_stage == 4 ) iori = (int)(270.0+double(ANGLE_INCR)/2.0+angle_degrees(triaxs,Vecf(0,0,0),pntaxs));
							iori = (iori / ANGLE_INCR) * ANGLE_INCR; // round to closest multiple of angle incr
							if( ori_stage == 2 || ori_stage == 4 ) iori -= ANGLE_INCR;
							newstage = false;
						} else {
							if( ori_stage == 1 || ori_stage == 3 ) iori += ANGLE_INCR;
							if( ori_stage == 2 || ori_stage == 4 ) iori -= ANGLE_INCR;
						}
						// TR << "IORI " << iori << std::endl;
						Vecf sicaxis = (rotation_matrix_degrees(-Vec(0,1,0),(Real)iori) * Vec(0,0,1)).normalized();
						int tmpcbc;
						double d = sicfast(pb,pa,cbb,cba,sicaxis,tmpcbc);

						double theta = iori;
						double gamma = theta-alpha;
						double x = d * sin(numeric::conversions::radians(gamma));
						double y = d * cos(numeric::conversions::radians(gamma));
						double w = x / sin(numeric::conversions::radians(alpha));
						double z = x / tan(numeric::conversions::radians(alpha));
						double dpnt = y+z;
						double dtri = w;
						double pntmn,trimn;
						if( w > 0 ) {
							pntmn = pntmnpos[ipnt/ANGLE_INCR+1];
							trimn = trimnpos[itri/ANGLE_INCR+1];
							if( dtri < trimn ) { ori_stage++; newstage=true; continue; };
							if( dpnt < pntmn ) { ori_stage++; newstage=true; continue; };
							int dp = (int)(dpnt-pntmn)*10+1;
							int dt = (int)(dtri-trimn)*10+1;
							// if(ipnt==18 && itri==72 && iori==276)
							//	TR << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbpos(ipnt/ANGLE_INCR+1,dp)
							//		<< "		" << dt << " " << dtri << " " << trimn << " " << tricbpos(itri/ANGLE_INCR+1,dt) << std::endl;
							if( option[tcdock::intra]() && dp <= 97	) tmpcbc += pntcbpos(ipnt/ANGLE_INCR+1,dp);
							if( option[tcdock::intra]() && dt <= 145 ) tmpcbc += tricbpos(itri/ANGLE_INCR+1,dt);
							// TR << "CHK " << dpnt << " " << pntmn << "		" << dtri << " " << trimn << std::endl;
						} else {
							pntmn = pntmnneg[ipnt/ANGLE_INCR+1];
							trimn = trimnneg[itri/ANGLE_INCR+1];
							if( dtri > trimn ) { ori_stage++; newstage=true; continue; };
							if( dpnt > pntmn ) { ori_stage++; newstage=true; continue; };
							int dp = (int)(-dpnt+pntmn)*10+1;
							int dt = (int)(-dtri+trimn)*10+1;
							// if(ipnt==18 && itri==72 && iori==276)
							//	TR << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbneg(ipnt/ANGLE_INCR+1,dp)
							//		<< "		" << dt << " " << dtri << " " << trimn << " " << pntcbneg(itri/ANGLE_INCR+1,dt) << std::endl;
							if( option[tcdock::intra]() && dp <= 97	) tmpcbc += pntcbneg(ipnt/ANGLE_INCR+1,dp);
							if( option[tcdock::intra]() && dt <= 145 ) tmpcbc += tricbneg(itri/ANGLE_INCR+1,dt);
						}

						surfdis3(ipnt/ANGLE_INCR+1,itri/ANGLE_INCR+1,iori/ANGLE_INCR+1) = d;
						cbcount(ipnt/ANGLE_INCR+1,itri/ANGLE_INCR+1,iori/ANGLE_INCR+1) = tmpcbc;
						gdis(ipnt+1,itri+1,iori+1) = d;
						gcbc(ipnt+1,itri+1,iori+1) = tmpcbc;
						// d = sicfast(t,p,sicaxis,cbcount);
						//std::cerr << "trial " << ipnt << " " << itri << " " << iori << " " << d << " " << tmpcbc << std::endl;
#ifdef USE_OPENMP
#pragma omp critical
#endif
						if(tmpcbc > cbmax) {
							cbmax = tmpcbc;
							mxiori = iori;
							mxipnt = ipnt;
							mxitri = itri;
							mxd = d;
						}
					}
					// utility_exit_with_message("testing");
				}
			}
			if(cbmax==0) utility_exit_with_message("tri or pnt too large, no contacts!");
			TR << "MAX3 " << mxipnt << " " << mxitri << " " << mxiori << " " << cbmax << " " << mxd << std::endl;
			// TR << "MAX " << mxipnt << " " << mxitri << " " << mxiori << " " << cbcount(mxipnt/ANGLE_INCR+1,mxitri/ANGLE_INCR+1,mxiori/ANGLE_INCR+1) << " " << surfdis3(mxipnt/ANGLE_INCR+1,mxitri/ANGLE_INCR+1,mxiori/ANGLE_INCR+1) << std::endl;
		}

		cbcount3 = cbcount;
		// int delta = ceil(3.0/(Real)ANGLE_INCR);
		// TR << "scanning results for local maxima +- " << delta*ANGLE_INCR << "°" << std::endl;
		// for(int ipnt = 1; ipnt <=	72/ANGLE_INCR; ++ipnt) {
		//	for(int itri = 1; itri <= 120/ANGLE_INCR; ++itri) {
		//		for(int iori = 1; iori <= 360/ANGLE_INCR; ++iori) {
		//			// std::cerr << ipnt << " " << itri << " " << iori << " " << cbcount(ipnt,itri,iori) << std::endl;
		//			int lmaxcb = 0;
		//			for(int dipnt = -delta; dipnt <= delta; ++dipnt) {
		//				for(int ditri = -delta; ditri <= delta; ++ditri) {
		//					for(int diori = -delta; diori <= delta; ++diori) {
		//						int i = (ipnt+dipnt-1+( 72/ANGLE_INCR)) % ( 72/ANGLE_INCR) + 1;
		//						int j = (itri+ditri-1+(120/ANGLE_INCR)) % (120/ANGLE_INCR) + 1;
		//						int k = (iori+diori-1+(360/ANGLE_INCR)) % (360/ANGLE_INCR) + 1;
		//						if( cbcount(i,j,k) > lmaxcb ) {
		//							lmaxcb = cbcount(i,j,k);
		//						}
		//					}
		//				}
		//			}
		//			if( cbcount(ipnt,itri,iori) >= lmaxcb ) cbcount3(ipnt,itri,iori) = cbcount(ipnt,itri,iori);
		//			if( ipnt==2 && itri==15 && iori==99 ) TR << "NEARMAX: " << cbcount(ipnt,itri,iori) << " " << lmaxcb << std::endl;
		//		}
		//	}
		// }
	}

	//TR << "main loop 2 over top 1000 3deg hits at 1deg" << std::endl;
	vector1<int> hitpnt,hittri,hitori,hitcbc;
	ANGLE_INCR = 1;
	if(ANGLE_INCR != 1) utility_exit_with_message("second ANGLE_INCR must be 1!!!");
	ObjexxFCL::FArray3D<float> surfdis1((Size)floor(72.0/ANGLE_INCR),(Size)floor(120.0/ANGLE_INCR),(Size)floor(360.0/ANGLE_INCR),0.0);
	ObjexxFCL::FArray3D<int>	 cbcount1((Size)floor(72.0/ANGLE_INCR),(Size)floor(120.0/ANGLE_INCR),(Size)floor(360.0/ANGLE_INCR),0);
	// 1 deg. sampling
	{
		// compute high/low min dis for pent and tri here, input to sicfast and don't allow any below
		TR << "precomputing pent-pent and tri-tri interactions every 1°" << std::endl;
		vector1<double>						pntmnpos(72/ANGLE_INCR,0),	 trimnpos(120/ANGLE_INCR,0);
		vector1<double>						pntmnneg(72/ANGLE_INCR,0),	 trimnneg(120/ANGLE_INCR,0);
		vector1<double>						pntdspos(72/ANGLE_INCR,0),	 tridspos(120/ANGLE_INCR,0);
		vector1<double>						pntdsneg(72/ANGLE_INCR,0),	 tridsneg(120/ANGLE_INCR,0);
		ObjexxFCL::FArray2D<int>	pntcbpos(72/ANGLE_INCR,97,0),tricbpos(120/ANGLE_INCR,145,0); // axes moving out
		ObjexxFCL::FArray2D<int>	pntcbneg(72/ANGLE_INCR,97,0),tricbneg(120/ANGLE_INCR,145,0);
		{
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
			for(int ipnt = 0; ipnt < 72; ipnt+=ANGLE_INCR) {
				Pose p;
#ifdef USE_OPENMP
#pragma omp critical
#endif
				p = pntinit;
				rot_pose(p,pntaxs,ipnt);
				vector1<Vecf> ppnt,cbp; // precompute these
				for(int i = 1; i <= (int)p.n_residue(); ++i) {
					cbp.push_back(Vecf(p.residue(i).xyz(2)));
					int const natom = (p.residue(i).name3()=="GLY") ? 4 : 5;
					for(int j = 1; j <= natom; ++j) ppnt.push_back(Vecf(p.residue(i).xyz(j)));
				}
				vector1<Vecf> ppn2(ppnt),cb2(cbp);
				pnt1_to_pnt2(ppn2,pntaxs,pntax2);
				pnt1_to_pnt2( cb2,pntaxs,pntax2);				
				int cbcount = 0;
				double const d = sicsafe(ppnt,ppn2,cbp,cb2,(pntax2-pntaxs).normalized(),cbcount,false);
				if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for pntpos! "+ObjexxFCL::string_of(ipnt));
				for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(pntax2-pntaxs).normalized();
				pntmnpos[ipnt/ANGLE_INCR+1] = -d/2.0/sin( angle_radians(pntax2,Vec(0,0,0),pntaxs)/2.0 );
				pntdspos[ipnt/ANGLE_INCR+1] = -d;
				pntcbpos(ipnt/ANGLE_INCR+1,1) = cbcount;
				// std::cerr << "compute CB" << std::endl;
				prune_cb_pairs_dis10(cbp,cb2);
				for(int i = 2; i <= 97; ++i) {
					for(vector1<Vecf>::iterator iv = cbp.begin(); iv != cbp.end(); ++iv) *iv = (*iv) + 0.1*pntaxs;
					for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) + 0.1*pntax2;
					int cbc = 0; for(Size j = 1; j <= cbp.size(); ++j) if(cbp[j].distance_squared(cb2[j]) < 100) cbc++;
					pntcbpos(ipnt/ANGLE_INCR+1,i) = cbc;
					if(cbc==0) break;
				}
			}
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
			for(int ipnt = 0; ipnt < 72; ipnt+=ANGLE_INCR) {
				Pose p;
#ifdef USE_OPENMP
#pragma omp critical
#endif
				p = pntinit;
				rot_pose(p,pntaxs,ipnt);
				vector1<Vecf> ppnt,cbp; // precompute these
				for(int i = 1; i <= (int)p.n_residue(); ++i) {
					cbp.push_back(Vecf(p.residue(i).xyz(2)));
					int const natom = (p.residue(i).name3()=="GLY") ? 4 : 5;
					for(int j = 1; j <= natom; ++j) ppnt.push_back(Vecf(p.residue(i).xyz(j)));
				}
				vector1<Vecf> ppn2(ppnt),cb2(cbp);
				pnt1_to_pnt2(ppn2,pntaxs,pntax2);
				pnt1_to_pnt2( cb2,pntaxs,pntax2);				
				int cbcount = 0;
				double const d = sicsafe(ppnt,ppn2,cbp,cb2,(pntaxs-pntax2).normalized(),cbcount,false);
				if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for pntneg! "+ObjexxFCL::string_of(ipnt));
				for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(pntaxs-pntax2).normalized();
				pntmnneg[ipnt/ANGLE_INCR+1] = d/2.0/sin( angle_radians(pntax2,Vec(0,0,0),pntaxs)/2.0 );
				pntdsneg[ipnt/ANGLE_INCR+1] = d;
				pntcbneg(ipnt/ANGLE_INCR+1,1) = cbcount;
				// std::cerr << "compute CB" << std::endl;
				prune_cb_pairs_dis10(cbp,cb2);
				for(int i = 2; i <= 97; ++i) {
					for(vector1<Vecf>::iterator iv = cbp.begin(); iv != cbp.end(); ++iv) *iv = (*iv) - 0.1*pntaxs;
					for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) - 0.1*pntax2;
					int cbc = 0; for(Size j = 1; j <= cbp.size(); ++j) if(cbp[j].distance_squared(cb2[j]) < 100) cbc++;
					pntcbneg(ipnt/ANGLE_INCR+1,i) = cbc;
					if(cbc==0) break;
				}
			}
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif			
			for(int itri = 0; itri < 120; itri+=ANGLE_INCR) {
				Pose t;
#ifdef USE_OPENMP
#pragma omp critical
#endif
				t = triinit;
				rot_pose(t,triaxs,itri);
				vector1<Vecf> ptri,cbt; // precompute these
				for(int i = 1; i <= (int)t.n_residue(); ++i) {
					cbt.push_back(Vecf(t.residue(i).xyz(2)));
					int const natom = (t.residue(i).name3()=="GLY") ? 4 : 5;
					for(int j = 1; j <= natom; ++j) ptri.push_back(Vecf(t.residue(i).xyz(j)));
				}
				vector1<Vecf> ptr2(ptri),cb2(cbt);
				tri1_to_tri2(ptr2,triaxs,triax2);
				tri1_to_tri2( cb2,triaxs,triax2);				
				int cbcount = 0;
				double const d = sicsafe(ptri,ptr2,cbt,cb2,(triax2-triaxs).normalized(),cbcount,false);
				if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for tripos! "+ObjexxFCL::string_of(itri));
				for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(triax2-triaxs).normalized();
				trimnpos[itri/ANGLE_INCR+1] = -d/2.0/sin( angle_radians(triax2,Vec(0,0,0),triaxs)/2.0 );
				tridspos[itri/ANGLE_INCR+1] = -d;
				tricbpos(itri/ANGLE_INCR+1,1) = cbcount;
				// std::cerr << "compute CB" << std::endl;
				prune_cb_pairs_dis10(cbt,cb2);
				for(int i = 2; i <= 145; ++i) {
					for(vector1<Vecf>::iterator iv = cbt.begin(); iv != cbt.end(); ++iv) *iv = (*iv) + 0.1*triaxs;
					for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) + 0.1*triax2;
					int cbc = 0; for(Size j = 1; j <= cbt.size(); ++j) if(cbt[j].distance_squared(cb2[j]) < 100) cbc++;
					tricbpos(itri/ANGLE_INCR+1,i) = cbc;
					if(cbc==0) break;
				}
			}
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif			
			for(int itri = 0; itri < 120; itri+=ANGLE_INCR) {
				Pose t;
#ifdef USE_OPENMP
#pragma omp critical
#endif
				t = triinit;
				rot_pose(t,triaxs,itri);
				vector1<Vecf> ptri,cbt; // precompute these
				for(int i = 1; i <= (int)t.n_residue(); ++i) {
					cbt.push_back(Vecf(t.residue(i).xyz(2)));
					int const natom = (t.residue(i).name3()=="GLY") ? 4 : 5;
					for(int j = 1; j <= natom; ++j) ptri.push_back(Vecf(t.residue(i).xyz(j)));
				}
				vector1<Vecf> ptr2(ptri),cb2(cbt);
				tri1_to_tri2(ptr2,triaxs,triax2);
				tri1_to_tri2( cb2,triaxs,triax2);				
				int cbcount = 0;
				double const d = sicsafe(ptri,ptr2,cbt,cb2,(triaxs-triax2).normalized(),cbcount,false);
				if( d > 0 ) utility_exit_with_message("d shouldn't be > 0 for trineg! "+ObjexxFCL::string_of(itri));
				for(vector1<Vecf>::iterator i = cb2.begin(); i != cb2.end(); ++i) *i = (*i) - d*(triaxs-triax2).normalized();
				trimnneg[itri/ANGLE_INCR+1] = d/2.0/sin( angle_radians(triax2,Vec(0,0,0),triaxs)/2.0 );
				tridsneg[itri/ANGLE_INCR+1] = d;
				tricbneg(itri/ANGLE_INCR+1,1) = cbcount;
				// std::cerr << "compute CB" << std::endl;
				prune_cb_pairs_dis10(cbt,cb2);
				for(int i = 2; i <= 145; ++i) {
					for(vector1<Vecf>::iterator iv = cbt.begin(); iv != cbt.end(); ++iv) *iv = (*iv) - 0.1*triaxs;
					for(vector1<Vecf>::iterator iv = cb2.begin(); iv != cb2.end(); ++iv) *iv = (*iv) - 0.1*triax2;
					int cbc = 0; for(Size j = 1; j <= cbt.size(); ++j) if(cbt[j].distance_squared(cb2[j]) < 100) cbc++;
					tricbneg(itri/ANGLE_INCR+1,i) = cbc;
					if(cbc==0) break;
				}
			}
			// {
			// 	using utility::file_basename;
			// 	TR << "dumping 1D stats: " << option[out::file::o]()+"/"+file_basename(option[in::file::s]()[ipntfile])+"_POS_1D.dat" << std::endl;
			// 	{ utility::io::ozstream out(option[out::file::o]()+"/"+file_basename(option[in::file::s]()[ipntfile])+"_POS_1D.dat");
			// 		for(int i = 1; i <=	72; ++i) out << i << " " << pntdspos[i] << " " << pntcbpos(i,1) << std::endl;
			// 		out.close(); }
			// 	{ utility::io::ozstream out(option[out::file::o]()+"/"+file_basename(option[in::file::s]()[ipntfile])+"_NEG_1D.dat");
			// 		for(int i = 1; i <=	72; ++i) out << i << " " << pntdsneg[i] << " " << pntcbneg(i,1) << std::endl;
			// 		out.close(); }
			// 	{ utility::io::ozstream out(option[out::file::o]()+"/"+file_basename(option[in::file::s]()[itrifile])+"_POS_1D.dat");
			// 		for(int i = 1; i <= 120; ++i) out << i << " " << tridspos[i] << " " << tricbpos(i,1) << std::endl;
			// 		out.close(); }
			// 	{ utility::io::ozstream out(option[out::file::o]()+"/"+file_basename(option[in::file::s]()[itrifile])+"_NEG_1D.dat");
			// 		for(int i = 1; i <= 120; ++i) out << i << " " << tridsneg[i] << " " << tricbneg(i,1) << std::endl;
			// 		out.close(); }
			// }
			
		}

		int top1000_3 = 0;
		{
			vector1<int> cbtmp;
			for(Size i = 0; i < cbcount3.size(); ++i) {
				if(cbcount3[i] > 0) cbtmp.push_back(cbcount3[i]);
			}
			std::sort(cbtmp.begin(),cbtmp.end());
			top1000_3 = cbtmp[max(1,(int)cbtmp.size()-499)];
			TR << "scanning top1000 with cbcount3 >= " << top1000_3 << std::endl;
		}
		// assuming ANGLE_INCR = 1!!!!!!!!!!!
		utility::vector1<vector1<int> > pntlmx,trilmx,orilmx;
		for(int ipnt = 0; ipnt < 72; ipnt+=3) {
			for(int itri = 0; itri < 120; itri+=3) {
				for(int iori = 0; iori < 360; iori+=3) {
					if( cbcount3(ipnt/3+1,itri/3+1,iori/3+1) >= top1000_3) {
						// TR << "checking around " << ipnt << " " << itri << " " << iori << " " << cbcount3(ipnt/3+1,itri/3+1,iori/3+1) << std::endl;
						vector1<int> pnt,tri,ori;
						for(int i = -1; i <= 1; ++i) pnt.push_back( (ipnt+i+ 72)% 72 );
						for(int j = -1; j <= 1; ++j) tri.push_back( (itri+j+120)%120 );
						for(int k = -1; k <= 1; ++k) ori.push_back( (iori+k+360)%360 );
						pntlmx.push_back(pnt);
						trilmx.push_back(tri);
						orilmx.push_back(ori);
					}
				}
			}
		}

		{
			int max1 = 0;
#ifdef USE_OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
			for(Size ilmx = 1; ilmx <= pntlmx.size(); ++ilmx) {
				if( (ilmx-1)%50==0 && ilmx!=1) TR << "	loop2 " << trifile << " " << (double(ilmx-1)/5) << "%% done, cbmax: " << max1 << std::endl;
				// TR << "checking around local max # " << ilmx << std::endl;
				int cbmax = 0, mxiori = 0, mxipnt = 0, mxitri = 0;
				for(vector1<int>::const_iterator pipnt = pntlmx[ilmx].begin(); pipnt != pntlmx[ilmx].end(); ++pipnt) {
					int ipnt = *pipnt;
					Pose p;
#ifdef USE_OPENMP
#pragma omp critical
#endif
					p = pntinit;
					rot_pose(p,pntaxs,(Real)ipnt);
					vector1<Vecf> pb,cbb; // precompute these
					for(int i = 1; i <= (int)p.n_residue(); ++i) {
						cbb.push_back(Vecf(p.residue(i).xyz(2)));
						int const natom = (p.residue(i).name3()=="GLY") ? 4 : 5;
						for(int j = 1; j <= natom; ++j) pb.push_back(Vecf(p.residue(i).xyz(j)));
					}
					for(vector1<int>::const_iterator pitri = trilmx[ilmx].begin(); pitri != trilmx[ilmx].end(); ++pitri) {
						int itri = *pitri;
						Pose t;
#ifdef USE_OPENMP
#pragma omp critical
#endif
						t = triinit;
						rot_pose(t,triaxs,(Real)itri);
						vector1<Vecf> pa,cba; // precompute these
						for(int i = 1; i <= (int)t.n_residue(); ++i) {
							cba.push_back(Vecf(t.residue(i).xyz(2)));
							int const natom = (t.residue(i).name3()=="GLY") ? 4 : 5;
							for(int j = 1; j <= natom; ++j) pa.push_back(Vecf(t.residue(i).xyz(j)));
						}
						for(vector1<int>::const_iterator piori = orilmx[ilmx].begin(); piori != orilmx[ilmx].end(); ++piori) {
							int iori = *piori;
							Vecf sicaxis = (rotation_matrix_degrees(-Vec(0,1,0),(Real)iori) * Vec(0,0,1)).normalized();
							int tmpcbc;
							double d = sicfast(pb,pa,cbb,cba,sicaxis,tmpcbc);
							double theta = iori;
							double gamma = theta-alpha;
							double x = d * sin(numeric::conversions::radians(gamma));
							double y = d * cos(numeric::conversions::radians(gamma));
							double w = x / sin(numeric::conversions::radians(alpha));
							double z = x / tan(numeric::conversions::radians(alpha));
							double dpnt = y+z;
							double dtri = w;
							double pntmn,trimn;
							if( w > 0 ) {
								pntmn = pntmnpos[ipnt/ANGLE_INCR+1];
								trimn = trimnpos[itri/ANGLE_INCR+1];
								if( dpnt < pntmn ) continue;
								if( dtri < trimn ) continue;								
								int dp = (int)(dpnt-pntmn)*10+1;
								int dt = (int)(dtri-trimn)*10+1;
								// if(ipnt==18 && itri==72 && iori==276)
								//	TR << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbpos(ipnt/ANGLE_INCR+1,dp)
								//		<< "		" << dt << " " << dtri << " " << trimn << " " << tricbpos(itri/ANGLE_INCR+1,dt) << std::endl;
								if( option[tcdock::intra]() && dp <= 97	) tmpcbc += pntcbpos(ipnt/ANGLE_INCR+1,dp);
								if( option[tcdock::intra]() && dt <= 145 ) tmpcbc += tricbpos(itri/ANGLE_INCR+1,dt);
								// TR << "CHK " << dpnt << " " << pntmn << "		" << dtri << " " << trimn << std::endl;
							} else {
								pntmn = pntmnneg[ipnt/ANGLE_INCR+1];
								trimn = trimnneg[itri/ANGLE_INCR+1];
								if( dpnt > pntmn ) continue;
								if( dtri > trimn ) continue;								
								int dp = (int)(-dpnt+pntmn)*10+1;
								int dt = (int)(-dtri+trimn)*10+1;
								// if(ipnt==18 && itri==72 && iori==276)
								//	TR << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbneg(ipnt/ANGLE_INCR+1,dp)
								//		<< "		" << dt << " " << dtri << " " << trimn << " " << pntcbneg(itri/ANGLE_INCR+1,dt) << std::endl;
								if( option[tcdock::intra]() && dp <= 97	) tmpcbc += pntcbneg(ipnt/ANGLE_INCR+1,dp);
								if( option[tcdock::intra]() && dt <= 145 ) tmpcbc += tricbneg(itri/ANGLE_INCR+1,dt);
							}
#ifdef USE_OPENMP
#pragma omp critical
#endif
							if(tmpcbc > cbmax) {
								cbmax = tmpcbc;
								mxiori = iori;
								mxipnt = ipnt;
								mxitri = itri;
							}
						}
					}
				}
#ifdef USE_OPENMP
#pragma omp critical
#endif
				{
					hitpnt.push_back(mxipnt);
					hittri.push_back(mxitri);
					hitori.push_back(mxiori);
					hitcbc.push_back(cbmax);
					if(cbmax > max1) max1 = cbmax;
				}
				// TR << "HIT " << ilmx << " " << mxipnt << " " << mxitri << " " << mxiori << " " << cbmax << std::endl;
			}
		}

		int top10,max1;
		{
			vector1<int> hittmp = hitcbc;
			std::sort(hittmp.begin(),hittmp.end());
			max1	= hittmp[hittmp.size()  ];
			top10 = hittmp[hittmp.size()-9];
			TR << "top10 " << top10 << std::endl;
		}
		for(Size ihit = 1; ihit <= hitcbc.size(); ++ihit) {
			if(hitcbc[ihit]==max1)	TR << "MAX1 " << hitpnt[ihit] << " " << hittri[ihit] << " " << hitori[ihit] << " " << hitcbc[ihit] << std::endl;
		}

		core::chemical::ResidueTypeSetCAP rs = core::chemical::ChemicalManager::get_instance()->residue_type_set("fa_standard");
// #ifdef USE_OPENMP
// #pragma omp parallel for schedule(dynamic,1)
// #endif		
		for(Size ihit = 1; ihit <= hitcbc.size(); ++ihit) {
			if(hitcbc[ihit] >= top10 ) {
			//if(hitcbc[ihit] >= max1 ) {				
				int ipnt = hitpnt[ihit];
				int itri = hittri[ihit];
				int iori = hitori[ihit];

				Pose p,t;
// #ifdef USE_OPENMP
// #pragma omp critical
// #endif
				{
					p = pntinit;
					t = triinit;					
				}
				rot_pose(p,pntaxs,(Real)ipnt);
				rot_pose(t,triaxs,(Real)itri);

				std::cerr << ihit << " geom" << std::endl;

				vector1<Vecf> pb,cbb; // precompute these
				for(int i = 1; i <= (int)p.n_residue(); ++i) {
					cbb.push_back(Vecf(p.residue(i).xyz(2)));
					int const natom = (p.residue(i).name3()=="GLY") ? 4 : 5;
					for(int j = 1; j <= natom; ++j) pb.push_back(Vecf(p.residue(i).xyz(j)));
				}
				vector1<Vecf> pa,cba; // precompute these
				for(int i = 1; i <= (int)t.n_residue(); ++i) {
					cba.push_back(Vecf(t.residue(i).xyz(2)));
					int const natom = (t.residue(i).name3()=="GLY") ? 4 : 5;
					for(int j = 1; j <= natom; ++j) pa.push_back(Vecf(t.residue(i).xyz(j)));
				}
				Vecf sicaxis = (rotation_matrix_degrees(-Vec(0,1,0),(Real)iori) * Vec(0,0,1)).normalized();
				int tmpcbc;
				double d = sicfast(pb,pa,cbb,cba,sicaxis,tmpcbc);

				double theta = iori;
				double gamma = theta-alpha;
				double x = d * sin(numeric::conversions::radians(gamma));
				double y = d * cos(numeric::conversions::radians(gamma));
				double w = x / sin(numeric::conversions::radians(alpha));
				double z = x / tan(numeric::conversions::radians(alpha));
				double dpnt = y+z;
				double dtri = w;
				double pntmn,trimn;
				trans_pose(p,dpnt*pntaxs);
				trans_pose(t,dtri*triaxs);
				if( w > 0 ) {
					pntmn = pntmnpos[ipnt/ANGLE_INCR+1];
					trimn = trimnpos[itri/ANGLE_INCR+1];
					if( dpnt < pntmn ) { TR << "WARNING: BAD HIT!!!"; continue; }
					if( dtri < trimn ) { TR << "WARNING: BAD HIT!!!"; continue; }
					int dp = (int)(dpnt-pntmn)*10+1;
					int dt = (int)(dtri-trimn)*10+1;
					// if(ipnt==18 && itri==72 && iori==276)
					//	TR << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbpos(ipnt/ANGLE_INCR+1,dp)
					//		<< "		" << dt << " " << dtri << " " << trimn << " " << tricbpos(itri/ANGLE_INCR+1,dt) << std::endl;
					if( option[tcdock::intra]() && dp <= 97	) tmpcbc += pntcbpos(ipnt/ANGLE_INCR+1,dp);
					if( option[tcdock::intra]() && dt <= 145 ) tmpcbc += tricbpos(itri/ANGLE_INCR+1,dt);
					// TR << "CHK " << dpnt << " " << pntmn << "		" << dtri << " " << trimn << std::endl;
				} else {
					pntmn = pntmnneg[ipnt/ANGLE_INCR+1];
					trimn = trimnneg[itri/ANGLE_INCR+1];
					if( dpnt > pntmn ) { TR << "WARNING: BAD HIT!!!"; continue; }
					if( dtri > trimn ) { TR << "WARNING: BAD HIT!!!"; continue; }
					int dp = (int)(-dpnt+pntmn)*10+1;
					int dt = (int)(-dtri+trimn)*10+1;
					// if(ipnt==18 && itri==72 && iori==276)
					//	TR << "DP " << dp << " " << dpnt << " " << pntmn << " " << pntcbneg(ipnt/ANGLE_INCR+1,dp)
					//		<< "		" << dt << " " << dtri << " " << trimn << " " << pntcbneg(itri/ANGLE_INCR+1,dt) << std::endl;
					if( option[tcdock::intra]() && dp <= 97	) tmpcbc += pntcbneg(ipnt/ANGLE_INCR+1,dp);
					if( option[tcdock::intra]() && dt <= 145 ) tmpcbc += tricbneg(itri/ANGLE_INCR+1,dt);
				}
				//Real sizefac = sasa_pnt/5.0 + sasa_tri/3.0;
				Real sizefac = tri_in.n_residue() + pnt_in.n_residue();

				std::string fname;
				fname	= ObjexxFCL::lead_zero_string_of(Size(Real(tmpcbc)/sizefac*1000000),10);
				fname += "_" + utility::file_basename(option[in::file::s]()[1]);
				fname += "_" + utility::file_basename(option[in::file::s]()[2]);
				fname += "_" + ObjexxFCL::lead_zero_string_of(tmpcbc,4);
				fname += "_" + ObjexxFCL::lead_zero_string_of(ipnt,3);
				fname += "_" + ObjexxFCL::lead_zero_string_of(itri,3);
				fname += "_" + ObjexxFCL::lead_zero_string_of(iori,3);
				fname += ".pdb.gz";
				// continue;
				Pose symm;
// #ifdef USE_OPENMP
// #pragma omp critical
// #endif
				{
					std::cerr << ihit << " sym" << std::endl;
					symm.append_residue_by_jump(t.residue(1),1);
					for(Size i = 2; i <= t.n_residue()/3; ++i) {
						if(symm.residue(i-1).is_upper_terminus()) core::pose::remove_upper_terminus_type_from_pose_residue(symm,i-1);
						if(	 t.residue(i	).is_lower_terminus()) core::pose::remove_lower_terminus_type_from_pose_residue(t,i);
						symm.append_residue_by_bond(t.residue(i));
					}
					symm.append_residue_by_jump(p.residue(1),1);
					for(Size i = 2; i <= p.n_residue()/5; ++i) symm.append_residue_by_bond(p.residue(i));
					// core::util::switch_to_residue_type_set(symm,"fa_standard");
					// for(Size i = 1; i <= symm.n_residue(); ++i) core::pose::replace_pose_residue_copying_existing_coordinates(symm,i,rs->name_map("ALA"));
					//TR << "making symm" << std::endl;					
				}

				core::pose::symmetry::make_symmetric_pose(symm);				
// #ifdef USE_OPENMP
// #pragma omp critical
// #endif
				{
					std::cerr << ihit << " dump " << fname << std::endl;				
					core::io::pdb::dump_pdb(symm,option[out::file::o]()+"/"+fname);
					TR << "HIT " << ihit << " p" << ipnt << " t" << itri << " o" << iori << " c" << hitcbc[ihit] << " c" << tmpcbc << " " << Real(tmpcbc)/sizefac << std::endl;
					TR << "    " << w << " "<< dtri << " " <<  trimnpos[itri/ANGLE_INCR+1] << " " << trimnneg[itri/ANGLE_INCR+1] << std::endl;
					TR << "    " << w << " "<< dpnt << " " <<  pntmnpos[ipnt/ANGLE_INCR+1] << " " << pntmnneg[ipnt/ANGLE_INCR+1] << std::endl;				
				}
					// {
					// 	std::cerr << trimnneg[itri/ANGLE_INCR+1] << " " << dtri << " " <<  triaxs.length() << std::endl;
					// 	Pose p(triinit),q(triinit);
					// 	rot_pose(p,triaxs,itri);
					// 	rot_pose(q,triaxs,itri);				
					// 	tri1_to_tri2(q,triaxs,triax2);
					// 	trans_pose(p,triaxs*trimnneg[itri/ANGLE_INCR+1]);
					// 	trans_pose(q,triax2*trimnneg[itri/ANGLE_INCR+1]);				
					// 	p.dump_pdb("test1.pdb");
					// 	q.dump_pdb("test2.pdb");				
					// }

			}
		}
	}
}


int main (int argc, char *argv[]) {
	register_options();
	devel::init(argc,argv);
	using namespace basic::options::OptionKeys;
	for(Size i = 2; i <= basic::options::option[in::file::s]().size(); ++i) {
		run(i,1);
	}
}




//
//








