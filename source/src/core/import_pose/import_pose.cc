// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/import_pose/import_pose.cc
/// @brief  various functions to construct Pose object(s) from PDB(s)
/// @details A temporary copy of the pose_from_pdb code from the demo directory.
/// Will be phased out in favor of file_data routines soon.
/// @author Sergey Lyskov

// Unit headers
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/import_pose_options.hh>

// Project headers
#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/carbohydrates/util.hh>

#include <core/util/metalloproteins_util.hh>

#include <core/pack/pack_missing_sidechains.hh>
#include <core/pack/optimizeH.hh>

#include <core/io/pdb/pdb_dynamic_reader.hh>
#include <core/io/pdb/pdb_dynamic_reader_options.hh>
#include <core/io/pdb/file_data.hh>
#include <core/io/pdb/file_data_options.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>

// External headers
#include <ObjexxFCL/string.functions.hh>
#include <boost/foreach.hpp>


namespace core {
namespace import_pose {

using core::Size;
using core::SSize;

using basic::T;
using basic::Error;
using basic::Warning;

using namespace core::io::pdb;
using namespace ObjexxFCL;

static thread_local basic::Tracer TR( "core.import_pose.import_pose" );

using utility::vector1;


void
read_additional_pdb_data(
	std::string const & s,
	pose::Pose & pose,
	io::pdb::FileData const &,
	bool read_fold_tree
)
{
	ImportPoseOptions options;
	read_additional_pdb_data( s, pose, options, read_fold_tree );
}

void
read_additional_pdb_data(
	std::string const & s,
	pose::Pose & pose,
	ImportPoseOptions const & options,
	bool read_fold_tree
)
{
	// split additional_pdb_data s into newlines
	utility::vector1< std::string > lines;
	Size start=0, i=0;

	while ( start < s.size() ) {
		if ( s[i] == '\n' || s[i] == '\r' /* || i==s.size()-1 */ ) {
			lines.push_back( std::string(s.begin()+start, s.begin()+i) );
			start = i+1;
		}
		i++;
		if ( i == s.size() ) {
			lines.push_back( std::string(s.begin()+start, s.begin()+i) );
			break;
		}
	}

	//Added by Daniel-Adriano Silva, used to read PDBinfo-LABEL
	//TR.Debug << "Setting PDBinfo-labels from PDB file." << std::endl;
	for ( i=1; i<= lines.size(); ++i ) {
		std::string const & line( lines[i] );
		if ( line.size() > 21 && line.substr(0,21) == "REMARK PDBinfo-LABEL:" ) {
			//Parse and split string
			utility::vector1 < std::string > remark_values;
			utility::vector1 < std::string > tmp_remark_values = utility::string_split(line, ' ');
			//Copy non-empty (i.e. !' ') elements to remark_values
			if ( tmp_remark_values .size() > 3 ) {
				for ( Size j=3; j<= tmp_remark_values.size(); ++j ) {
					if ( tmp_remark_values[j] != "" ) {
						remark_values.push_back(tmp_remark_values[j]);
					}
				}
			}
			//Check that we have at least two elements left ([1]=index, [2-n]=PDBinfo-labels)
			if ( remark_values.size() > 1 ) {
				core::Size tmp_ndx=atoi(remark_values[1].c_str());
				if ( tmp_ndx <= pose.total_residue() ) {
					for ( Size j=2; j<= remark_values.size(); ++j ) {
						pose.pdb_info()->add_reslabel(tmp_ndx,remark_values[j]);
					}
				} else {
					TR.Fatal << "pose_io:: PDBinfo-LABEL io failure: " << line << ' ' << pose.total_residue()  << std::endl;
				}
			} else {
				TR.Fatal << "pose_io:: PDBinfo-LABEL io failure: " << line << ' ' << pose.total_residue() << std::endl;
			}
		}
	}

	if ( (!read_fold_tree) && (!options.fold_tree_io()) ) return;

	for ( i=1; i<= lines.size(); ++i ) {
		std::string const & line( lines[i] );
		// Look for fold_tree info
		if ( line.size() >= 16 && line.substr(0,16) == "REMARK FOLD_TREE" ) {
			std::istringstream l( line );
			std::string tag;
			kinematics::FoldTree f;
			l >> tag >> f;
			if ( !l.fail() && Size(f.nres()) == pose.total_residue() ) {
				TR << "setting foldtree from pdb file: " << f << std::endl;
				pose.fold_tree( f );
			} else {
				TR.Fatal << "pose_io:: foldtree io failure: " << line << ' ' << pose.total_residue()
					<< ' ' << f << std::endl;
				utility_exit();
			}
		}
	}
}


pose::PoseOP pose_from_pdb(std::string const & filename, bool read_fold_tree)
{
	pose::PoseOP pose( new pose::Pose() );
	pose_from_pdb( *pose, filename, read_fold_tree);
	return pose;
}


pose::PoseOP pose_from_pdb(chemical::ResidueTypeSet const & residue_set, std::string const & filename,  bool read_fold_tree)
{
	pose::PoseOP pose( new pose::Pose() );
	pose_from_pdb( *pose, residue_set, filename, read_fold_tree);
	return pose;
}

void
pose_from_pdb(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filenames_string,
	bool read_fold_tree
)
{
	ImportPoseOptions options;
	pose_from_pdb(pose, residue_set, filenames_string, options, read_fold_tree);
}

void
pose_from_pdb(
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filenames_string,
	ImportPoseOptions const & options,
	bool read_fold_tree
)
{
	utility::vector1<std::string> filenames = utility::split(filenames_string);

	std::string res;

	BOOST_FOREACH ( std::string filename, filenames ) {
		utility::io::izstream file( filename );
		if ( !file ) {
			TR.Error << "PDB File:" << filename << " not found!" << std::endl;
			utility_exit_with_message( "Cannot open PDB file \"" + filename + "\"" );
		} else {
			TR.Debug << "read file: " << filename << std::endl;
		}
		utility::slurp( file, res );
	}

	//fpd If the conformation is not of type core::Conformation, reset it
	conformation::ConformationOP conformation_op( new conformation::Conformation() );
	if ( !pose.conformation().same_type_as_me( *conformation_op, true ) ) {
		pose.set_new_conformation( conformation_op );
	}

	io::pdb::FileData fd = core::import_pose::PDB_DReader::createFileData(res, options);
	if ( fd.filename == "" ) {
		fd.filename = utility::join(filenames, "_");
	}
	build_pose(fd, pose, residue_set, options);

	// set secondary structure for centroid PDBs
	if ( residue_set.name() == core::chemical::CENTROID ) {
		core::pose::set_ss_from_phipsi( pose );
	}

	// check for foldtree info
	read_additional_pdb_data( res, pose, options, read_fold_tree );
}


void
pose_from_pdb(
	core::pose::Pose & pose,
	std::string const & filename,
	bool read_fold_tree
) {
	ImportPoseOptions options;
	pose_from_pdb( pose, filename, options, read_fold_tree );
}

void
pose_from_pdb(
	pose::Pose & pose,
	std::string const & filename,
	ImportPoseOptions const & options,
	bool read_fold_tree
) {
	using namespace chemical;

	ResidueTypeSetCOP residue_set( options.centroid() ?
		ChemicalManager::get_instance()->residue_type_set( CENTROID ) :
		ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )
	);

	core::import_pose::pose_from_pdb( pose, *residue_set, filename, options, read_fold_tree );
}

utility::vector1< core::pose::PoseOP >
poseOPs_from_pdbs(
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree
) {
	ImportPoseOptions options;
	return poseOPs_from_pdbs( filenames, options, read_fold_tree );
}

utility::vector1< core::pose::PoseOP >
poseOPs_from_pdbs(
	utility::vector1< std::string > const & filenames,
	ImportPoseOptions const & options,
	bool read_fold_tree
) {
	using namespace chemical;

	ResidueTypeSetCOP residue_set( options.centroid() ?
		ChemicalManager::get_instance()->residue_type_set( CENTROID ) :
		ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )
	);

	return core::import_pose::poseOPs_from_pdbs( *residue_set, filenames, options, read_fold_tree );
}

/// @details Only returns full-atom poses
utility::vector1< core::pose::Pose >
poses_from_pdbs(
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree
) {
	using namespace chemical;
	ResidueTypeSetCOP residue_set
		( ChemicalManager::get_instance()->residue_type_set( FA_STANDARD ) );

	return core::import_pose::poses_from_pdbs( *residue_set, filenames, read_fold_tree );
}

utility::vector1< core::pose::Pose >
poses_from_pdbs(
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< std::string > const & filenames,
	bool read_fold_tree
) {
	using namespace chemical;

	using std::string;
	using utility::vector1;
	using core::pose::Pose;

	ImportPoseOptions options;

	vector1< Pose > poses;
	typedef vector1< string >::const_iterator vec_it;
	for ( vec_it it = filenames.begin(), end = filenames.end(); it != end; ++it ) {
		Pose pose;
		core::import_pose::pose_from_pdb( pose, residue_set, *it, options, read_fold_tree );
		poses.push_back( pose );
	}

	return poses;
}

void
read_all_poses(
	const utility::vector1<std::string>& filenames,
	utility::vector1<core::pose::Pose>* poses
)
{
	using core::pose::Pose;
	using std::string;
	using utility::vector1;
	debug_assert(poses);

	for ( vector1<string>::const_iterator i = filenames.begin(); i != filenames.end(); ++i ) {
		Pose pose;
		pose_from_pdb(pose, *i);
		poses->push_back(pose);
	}
}

utility::vector1< core::pose::PoseOP >
poseOPs_from_pdbs(
	chemical::ResidueTypeSet const & residue_set,
	utility::vector1< std::string > const & filenames,
	ImportPoseOptions const & options,
	bool read_fold_tree
) {
	using namespace chemical;

	using std::string;
	using utility::vector1;
	using core::pose::Pose;

	vector1< pose::PoseOP > poses;
	typedef vector1< string >::const_iterator vec_it;
	for ( vec_it it = filenames.begin(), end = filenames.end(); it != end; ++it ) {
		pose::PoseOP pose( new pose::Pose );
		core::import_pose::pose_from_pdb( *pose, residue_set, *it, options, read_fold_tree );
		poses.push_back( pose );
	}

	return poses;
}


void
pose_from_pdb(
	utility::vector1< pose::Pose > & poses,
	std::string const & filename,
	bool read_fold_tree
)
{
	using namespace chemical;
	ResidueTypeSetCOP residue_set(
		ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )
	);
	core::import_pose::pose_from_pdb( poses, *residue_set, filename, read_fold_tree );
}


void
pose_from_pdb(
	utility::vector1< pose::Pose > & poses,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	bool read_fold_tree
)
{
	ImportPoseOptions options;
	pose_from_pdb( poses, residue_set, filename, options, read_fold_tree );
}

void
pose_from_pdb(
	utility::vector1< pose::Pose > & poses,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename,
	ImportPoseOptions const & options,
	bool read_fold_tree
)
{
	// Size fsize;
	pose::Pose pose;
	std::string all_lines, sub_lines;

	utility::io::izstream file( filename );
	if ( !file ) {
		TR.Error << "File:" << filename << " not found!" << std::endl;
		utility_exit_with_message( "Cannot open file " + filename );
	} else {
		TR.Debug << "read file: " << filename << std::endl;
	}

	utility::slurp( file, all_lines );
	// Hacky way to make sure that we find all models. I blame society.
	all_lines = "\n" + all_lines;

	Size pos1 = 0;
	Size pos2 = 0;
	Size n_models = 0;

	// count the number of poses will be reading
	while ( pos1 != std::string::npos )
			{
		pos1 = all_lines.find( "\nMODEL ", pos1 );
		if ( pos1 != std::string::npos ) {
			++n_models;
			++pos1;
		}
	}

	TR.Debug << "Reading " << n_models << " poses." << std::endl;
	// make space for all of our poses
	if ( n_models == 0 ) n_models = 1;
	poses.reserve( poses.size() + n_models );

	pos1 = 0;

	pos1 = all_lines.find( "\nMODEL ", pos1 );

	if ( pos1 != std::string::npos ) {
		pos2 = 0;
		while ( pos2 != std::string::npos ) {
			// set pos1 == "M" in MODEL
			++pos1;

			// pos2 = position of newline character, start somewhere after pos1
			pos2 = all_lines.find( "\nMODEL ", pos1);
			sub_lines = all_lines.substr( pos1, pos2-pos1 ) ;
			pos1 = pos2;

			io::pdb::FileData fd = core::import_pose::PDB_DReader::createFileData( sub_lines, options );
			fd.filename = filename;
			build_pose( fd, pose, residue_set, options );

			// check for foldtree info
			core::import_pose::read_additional_pdb_data( sub_lines, pose, options, read_fold_tree);
			poses.push_back( pose );
		}
	} else {
		FileData fd = core::import_pose::PDB_DReader::createFileData( all_lines, options );
		if ( fd.filename == "" ) {
			fd.filename = filename;
		}
		core::import_pose::build_pose( fd, pose, residue_set, options );
		// check for foldtree info
		core::import_pose::read_additional_pdb_data( all_lines, pose, options, read_fold_tree);
		poses.push_back( pose );
	}
}

void
pose_from_pdbstring(
	pose::Pose & pose,
	std::string const & pdbcontents,
	ImportPoseOptions const & options,
	std::string const & filename
)
{
	io::pdb::FileData fd = import_pose::PDB_DReader::createFileData( pdbcontents, options );
	fd.filename = filename;
	chemical::ResidueTypeSetCOP residue_set
		( chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD ) );
	core::import_pose::build_pose( fd, pose, *residue_set, options );

}

void
pose_from_pdbstring(
	pose::Pose & pose,
	std::string const & pdbcontents,
	std::string const & filename
)
{
	ImportPoseOptions options;
	pose_from_pdbstring( pose, pdbcontents, options, filename );
}


void
pose_from_pdbstring(
	pose::Pose & pose,
	std::string const & pdbcontents,
	chemical::ResidueTypeSet const & residue_set,
	std::string const & filename
){
	ImportPoseOptions options;
	pose_from_pdbstring( pose, pdbcontents, residue_set, options, filename );
}

void
pose_from_pdbstring(
	pose::Pose & pose,
	std::string const & pdbcontents,
	chemical::ResidueTypeSet const & residue_set,
	ImportPoseOptions const & options,
	std::string const & filename
){
	io::pdb::FileData fd = import_pose::PDB_DReader::createFileData( pdbcontents, options );
	fd.filename = filename;
	core::import_pose::build_pose( fd, pose, residue_set, options);
}

void pose_from_pdb_stream(
	pose::Pose & pose,
	std::istream & pdb_stream,
	std::string const & filename,
	ImportPoseOptions const & options
){
	std::string pdb_file_contents;
	utility::slurp( pdb_stream, pdb_file_contents );

	chemical::ResidueTypeSetCOP residue_set( chemical::ChemicalManager::get_instance()->residue_type_set( options.residue_type_set()) );
	pose_from_pdbstring( pose, pdb_file_contents, *residue_set, options, filename);
	read_additional_pdb_data( pdb_file_contents, pose, options, options.read_fold_tree() );
}

void
centroid_pose_from_pdb(
	pose::Pose & pose,
	std::string const & filename,
	bool read_fold_tree
)
{
	using namespace chemical;
	ResidueTypeSetCOP residue_set
		( ChemicalManager::get_instance()->residue_type_set( CENTROID ) );

	core::import_pose::pose_from_pdb( pose, *residue_set, filename, read_fold_tree);
}

void build_pose(
	io::pdb::FileData & fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set
)
{
	ImportPoseOptions options; // read from the command line
	build_pose( fd, pose, residue_set, options);
}

/// @brief Build Rosetta 3 Pose object from FileData.
void build_pose(
	io::pdb::FileData & fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	ImportPoseOptions const & options
)
{
	TR.Debug << "build_pose..." << std::endl;
	build_pose_as_is( fd, pose, residue_set, options);
	TR.Debug << "build_pose... Ok." << std::endl;
}

void build_pose_as_is2(
	io::pdb::FileData & fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	id::AtomID_Mask & missing,
	ImportPoseOptions const & options );

// "super-simple" (C) by Phil
/// @brief Try to Build pose object from pdb 'as-is'. PDB file must be _really_ clean.
void build_pose_as_is(
	io::pdb::FileData & fd,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	ImportPoseOptions const & options
)
{
	id::AtomID_Mask missing( false );

	io::pdb::build_pose_as_is1( fd, pose, residue_set, missing, options );
	build_pose_as_is2( fd, pose, residue_set, missing, options );
}

void build_pose_as_is2(
	io::pdb::FileData & /*fd*/,
	pose::Pose & pose,
	chemical::ResidueTypeSet const & residue_set,
	id::AtomID_Mask & missing,
	ImportPoseOptions const & options
)
{
	using namespace chemical;
	using namespace conformation;

	core::pose::PDBInfoOP pdb_info( new core::pose::PDBInfo(*pose.pdb_info()) );

	if ( !options.skip_set_reasonable_fold_tree() ) {
		set_reasonable_fold_tree( pose );
	}

	// optimize H if using a full-atom residue type set, and no_optH is not specified
	if ( residue_set.name() == FA_STANDARD ) {
		//if pack_missing_density specified, repack residues w/ missing density
		if ( options.pack_missing_sidechains()  && ! options.membrane() ) {
			pack::pack_missing_sidechains( pose, missing );
		}
		// optimize H if using a fullatom residue type set, and no_optH is not specified
		if ( !options.no_optH() ) {
			pack::optimize_H_and_notify( pose, missing );
		}
	}


	// If pose contains carbohydrate residues, assure that their virtual atoms have the correct coordinates.
	if ( pose.conformation().contains_carbohydrate_residues() ) {
		for ( uint i = 1, n_residues = pose.total_residue(); i <= n_residues; ++i ) {
			ResidueType const & res_type = pose.residue_type(i);
			if ( res_type.is_carbohydrate() ) {
				pose::carbohydrates::align_virtual_atoms_in_carbohydrate_residue(pose, i);
			}
		}
	}

	// If the user has set appropriate flags, check whether the pose contains metal ions,
	// and automatically set up covalent bonds and constraints to them.
	if ( options.set_up_metal_bonds() ) {
		core::util::auto_setup_all_metal_bonds(pose, options.metal_bond_LJ_multiplier(), true);
		if ( options.set_up_metal_constraints() ) {
			core::util::auto_setup_all_metal_constraints( pose, options.metal_bond_dist_constraint_multiplier(), options.metal_bond_angle_constraint_multiplier() );
		}
	}

	pose.pdb_info( pdb_info );
}


/// @brief Convert a jump to the given Residue in a FoldTree to a chemical Edge.
/// @param <end_resnum>: the Residue index of the end of the Jump\n
/// <ft>: the FoldTree being modified
/// @return true if a chemical connection was found and the passed FoldTree was changed
/// @note This function is used by set_reasonable_foldtree() for the primary purpose of rebuilding default FoldTrees in
/// cases of branching or chemical conjugation.
/// @author Labonte <JWLabonte@jhu.edu>
/// @author Joseph Harrison
bool
change_jump_to_this_residue_into_chemical_edge(
	core::uint end_resnum,
	core::pose::Pose const & pose,
	core::kinematics::FoldTree & ft )
{
	using namespace std;
	using namespace core::chemical;
	using namespace conformation;

	Residue const & end_residue( pose.residue( end_resnum ) );
	bool found_chemical_connection( false );

	if ( TR.Trace.visible() ) {
		TR.Trace << "Searching for ResidueConnection from residue " << end_resnum << endl;
	}
	// Loop through the ResidueConnections of the end Residue of the old Jump to find and attach the other end of the
	// new chemical Edge.
	Size const n_connections( end_residue.type().n_residue_connections() );
	for ( uint conn_index( 1 ); conn_index <= n_connections; ++conn_index ) {
		if ( end_residue.connection_incomplete( conn_index ) ) {
			continue;  // Allow incomplete connections for design.
		}
		uint const start_resnum( end_residue.connect_map( conn_index ).resid() );
		if ( start_resnum < end_resnum ) {
			Residue const & start_residue( pose.residue( start_resnum ) );
			// Ensure that the connection is either a polymer branch or a ligand of the same chain.
			if ( ( start_residue.is_branch_point() && end_residue.is_branch_lower_terminus() ) ||
					( ( ! end_residue.is_polymer() ) && start_residue.chain() == end_residue.chain() ) ) {
				uint start_atm_index = start_residue.connect_atom( end_residue );
				uint end_atm_index = end_residue.connect_atom( start_residue );

				string start_atm_name = start_residue.atom_name( start_atm_index );
				string end_atm_name = end_residue.atom_name( end_atm_index );

				if ( TR.Trace.visible() ) {
					TR.Trace << "Adding  EDGE " <<
						start_resnum << ' ' << end_resnum << ' ' << start_atm_name << ' ' << end_atm_name <<
						" to the new FoldTree." << endl;
				}
				ft.add_edge( start_resnum, end_resnum, start_atm_name, end_atm_name );

				found_chemical_connection = true;
				break;
			}
		}
	}
	return found_chemical_connection;
}


/// @details All ligand residues and polymer branches have been appended by a jump.
/// This method creates a new fold tree without jumps through ligands, using CHEMICAL edges instead.
void
set_reasonable_fold_tree( pose::Pose & pose )
{
	// An empty pose doesn't have jumps through ligands.
	// (Will encounter a SegFault otherwise)
	if ( pose.total_residue() == 0 ) {
		return;
	}

	using namespace std;
	using namespace core::chemical;
	using namespace core::kinematics;

	if ( TR.Debug.visible() ) {
		TR.Debug << "original fold tree: " << pose.fold_tree() << endl;
	}

	FoldTree const & origft( pose.fold_tree() );
	FoldTree newft;

	// Copy the original fold tree edges to the new fold tree, possibly replacing
	// some jumps to ligand residues with ligand-ligand chemical bonds.
	// As a result, jumps must be renumbered.
	uint last_jump_id( 0 );
	bool prevent_forward_edge( false );  // Use to prevent a custom reverse edge from being overwritten.
	for ( FoldTree::const_iterator i( origft.begin() ), i_end( origft.end() ); i != i_end; ++i ) {
		Edge e( *i );
		if ( TR.Trace.visible() ) {
			TR.Trace << "checking if if " << e << " is reasonable for this Pose..." << endl;
		}
		if ( e.is_jump() ) {
			uint const ii( e.stop() );  // the residue position at the end of the jump
			// Is this jump to a covalent ligand residue or forward-direction polymer branch?
			if ( ( pose.residue( ii ).is_branch_lower_terminus() ) /* polymer branch in forward direction */ ||
					( ! pose.residue( ii ).is_polymer() ) /* covalent ligand residue */ ) {
				if ( TR.Trace.visible() ) {
					TR.Trace << e << " was initially set as a Jump to a branch lower terminus or ligand." << endl;
				}
				bool jump_was_changed( change_jump_to_this_residue_into_chemical_edge( ii, pose, newft ) );

				if ( jump_was_changed ) {
					continue;  // The Edge has already been added to newft by the function above.
				} else /* We couldn't find a chemical connection. */ {
					if ( pose.residue_type( ii ).n_residue_connections() > 0 ) {
						TR.Warning << "Can't find a chemical connection for residue " << ii << " " <<
							pose.residue_type( ii ).name() << endl;
					}
				}
			} else {
				if ( TR.Trace.visible() ) {
					TR.Trace << e << " was initially set as a Jump to a lower terminus." << endl;
				}
				// Figure out if the jump should connect to the LAST residue of this Edge/chain and not the 1st,
				// i.e., is this a C-terminal connection?

				Edge const next_e( *(i + 1) );
				if ( ii == uint( next_e.start() ) && pose.residue( ii ).is_lower_terminus() ) {
					uint const jj( next_e.stop() );  // the residue position at the end of the next Edge
					if ( pose.residue_type( jj ).is_branch_lower_terminus() ) {
						bool jump_was_changed( change_jump_to_this_residue_into_chemical_edge( jj, pose, newft ) );

						if ( jump_was_changed ) {
							// The CHEMICAL Edge has already been added to newft by the function above.
							// Now, we must add a reverse-direction Edge and skip making the usual forward-direction
							// Edge.
							// The Edge has already been added to newft by the function above.
							if ( TR.Trace.visible() ) {
								TR.Trace << "Adding  EDGE " << jj << ' ' << ii << " -1  to the new FoldTree." << endl;
							}
							newft.add_edge( jj, ii, Edge::PEPTIDE);
							prevent_forward_edge = true;
							continue;  // Skip the add_edge() method below; we are done here.
						} else /* We couldn't find a chemical connection. */ {
							if ( pose.residue_type( jj ).n_residue_connections() > 0 ) {
								TR.Warning << "Can't find a chemical connection for residue " << jj << " " <<
									pose.residue_type( jj ).name() << endl;
							}
						}
					}
				}
			}
			// It's just a normal inter-chain jump or we couldn't find a chemical bond to make. Increment the label.
			e.label() = ++last_jump_id;
		} else {
			if ( TR.Trace.visible() ) {
				TR.Trace << e << " was not a Jump." << endl;
			}
		}
		if ( ! prevent_forward_edge ) {
			if ( TR.Trace.visible() ) {
				TR.Trace << "Adding " << e << " to the new FoldTree." << endl;
			}
			newft.add_edge( e );
		} else {
			prevent_forward_edge = false;  // Adding forward Edge was skipped; reset for next Edge.
		}
	}  // next Edge

	runtime_assert( newft.size() > 0 || pose.total_residue() == 0 );  // A valid fold tree must have >= 1 edges.

	if ( TR.Debug.visible() ) {
		TR.Debug << "new fold tree " << newft << endl;
	}

	pose.fold_tree( newft );
}


} // namespace import_pose
} // namespace core
