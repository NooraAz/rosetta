#!/usr/bin/python

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   /GUIs/pyrosetta_toolkit/window_modules/clean_pdb/fix_pdb_window.py
## @brief  Simple GUI for cleaning pdbs for rosetta.
## @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#Rosetta imports
from rosetta.basic.options import get_boolean_option
#Python Imports
import os
import re

#Tkinter Imports
from Tkinter import *
import tkFileDialog
import tkMessageBox

#Toolkit Imports
from window_main import global_variables
from modules.PythonPDB import PythonPDB
from modules.SQLPDB import SQLPDB
from window_modules.ligand_ncaa_ptm_manager.ligand_ncaa_ptm_manager import ligand_ncaa_ptm_manager
from window_modules.options_system.OptionSystemManager import OptionSystemManager

class FixPDBWindow:
     """
     Quick and dirty cleaning of PDB.
     Requires all these main toolkit classes due to it's use of the ligand_ncaa_ptm_manager
     """
     def __init__(self, input_class, score_class, pose, starting_file_path=""):
          self.input_class = input_class
          
          self.ncaa_manager = ligand_ncaa_ptm_manager(input_class, score_class, pose)
          self.PDBlistPath = StringVar()
          if starting_file_path:
               self.PDBlistPath.set(starting_file_path)

          #These control fixonepdb protocol.  (Go button)
          self.remove_hetatm_var= IntVar(); self.remove_alternates_var = IntVar(); self.change_occupancies_var=IntVar()
          self.remove_hetatm_var.set(False); self.remove_alternates_var.set(False); self.change_occupancies_var.set(True)
          self.remove_waters_var = IntVar(); self.remove_waters_var.set(True)
          self.check_rosetta_var=IntVar(); self.check_rosetta_var.set(True)
          self.replace_res_and_atoms_var =IntVar(); self.replace_res_and_atoms_var.set(True)
          #Arrays determined by various functions.
          self.recognized_aa=[];      #[string three_letter_code] All Codes recognized by rosetta.
          self.on_by_default_aa=[];   #[string three_letter_code] All Codes ON by DEFAULT in rosetta.
          self.get_on_by_default_aa()
          self.get_recognized_aa()
          
          self.unrecognized_aa=[];    #[string three_letter_code] All residues in a PDB that are not recognized by any three letter code specified in params/paths
          self.off_by_default_aa = [];#[string three_letter_code] All residues found that are off by default.
          
          self.ignore_unrecognized_option_set = False; #Is the ignore_unrecognized_residue option already set?  If not, we warn the user before proceeding with PDB load.
          
          
          self.cleaned_pdb_path = StringVar(); #Path where the SINGLE pdb was written to.
          
          
     def runfixPDBWindow(self, m, r, c):
          """
          Shows the Window.  Very simple window, but works well enough for what it does. Can use a PDB or PDBList. 
          m = Main Window
          r = Row
          c = Column
          """
          self.defaultdir = os.getcwd(); #Needs to get this from toolkit...
          fixWindow = Toplevel(m)
          fixWindow.title("Setup PDBs for Rosetta")

          #remove_hetatm_var, remove_alternates_var, change_occupancies_var = 1
          #fixDictionary = dict()
          self.pathEntry = Entry(fixWindow, textvariable=self.PDBlistPath)
          self.pathbutton_ = Button(fixWindow, text = "PDBLIST or PDB", command = lambda: self.PDBlistPath.set(self.getfile()))
          
          self.remH20 = Checkbutton(fixWindow, variable = self.remove_waters_var, text = "Remove Water")
          self.remHet = Checkbutton(fixWindow, variable=self.remove_hetatm_var, text = "Remove HETATM")
          
          self.remAlt = Checkbutton(fixWindow, variable = self.remove_alternates_var, text = "Renumber all chains + incorporate insertions")
          self.chaOcc = Checkbutton(fixWindow, variable = self.change_occupancies_var, text = "Change occupancies to 1")
          self.check_unrecognized_button=Checkbutton(fixWindow, variable = self.check_rosetta_var, text="Check for Rosetta unrecognized/off residue types")
          self.replace_res_and_atoms_button=Checkbutton(fixWindow, variable = self.replace_res_and_atoms_var, text = "Replace known unrecognized residues and atoms")
          self.gobutton_ = Button(fixWindow, text = "GO", command = lambda: self.runFixPDB())
          self.loadbutton = Button(fixWindow, text = "Load Cleaned PDB", command = lambda: self.load_cleaned_pdb())
          self.ignore_unrecognized = Button(fixWindow, text = "Set -ignore_unrecognized_res option", command = lambda: self.set_ignore_unrecognized())
          self.pathEntry.grid(row = r, column=c)
          self.pathbutton_.grid(row=r, column = c+1)
          self.remH20.grid(row=r, column=c+2, sticky=W)
          self.remHet.grid(row=r+1, column=c+2, sticky=W)
          self.chaOcc.grid(row=r+2, column=c+2, sticky=W)
          self.remAlt.grid(row=r+3, column=c+2, sticky=W)
          self.check_unrecognized_button.grid(row=r+4, column=c+2, sticky=W)
          self.replace_res_and_atoms_button.grid(row=r+5, column=c+2, sticky=W)
          
          self.gobutton_.grid(row=r+6, column=c,columnspan=2, sticky=W+E);self.loadbutton.grid(row=r+6, column=c+2, sticky=W+E, pady=7)
          
          self.loadbutton.config(state=DISABLED)
          self.ignore_unrecognized.grid(row=r+7, column=c, columnspan=3, sticky=W+E)
          
     def get_recognized_aa(self):
          """
          Gets recognized 3 letter codes from param/patches
          """
          for name in self.ncaa_manager.param_name_map:
               param = self.ncaa_manager.param_name_map[name]
               self.recognized_aa.append(param.three_letter_name.get())

          for name in self.ncaa_manager.patch_name_map:
               param = self.ncaa_manager.patch_name_map[name]
               self.recognized_aa.append(param.three_letter_name.get())
               
          #Residues Rosetta CAN recognize
          for aa in ["DA", "DT", "DG", "DC", "WAT", "TP3"]: self.recognized_aa.append(aa)
          
     def get_on_by_default_aa(self):
          """
          Gets 3 letter codes that are on by default from param/patches
          Sets a dictionary so that these can be auto turned on?
          """
          for name in self.ncaa_manager.param_name_map:
               param = self.ncaa_manager.param_name_map[name]
               if param.rosetta_read_state:
                    self.on_by_default_aa.append(param.three_letter_name.get())

          for name in self.ncaa_manager.patch_name_map:
               param = self.ncaa_manager.patch_name_map[name]
               if param.rosetta_read_state:
                    self.on_by_default_aa.append(param.three_letter_name.get())
          
          #Residues that ARE on by default
          for aa in ["DA", "DT", "DG", "DC", "WAT", "TP3"]: self.on_by_default_aa.append(aa)
          
          
     def load_cleaned_pdb(self):
         
          
          #Any unrecognized - Warn - Rosetta may be unable to recognize these
          self.ignore_unrecognized_option_set = get_boolean_option('in:ignore_unrecognized_res')
          keep_going = True
          if not self.ignore_unrecognized_option_set:
               unique = dict()
               for aa in self.unrecognized_aa: unique[aa]=0
               if unique:
                    aa_string = ""
                    for aa in sorted(unique):
                         aa_string = aa_string+":"+aa
                    keep_going =  tkMessageBox.askyesno(title="Continue?", message = "Rosetta may be unable to read these residues: "+aa_string+". Continue?")
          if not keep_going: return
          
          #Use off-by-default - Enable NCAA If there are any found
          unique = dict()
          for aa in self.off_by_default_aa: unique[aa]=0
          if unique:
               print "Enabling NCAA found in PDB.  These are normally off by default.  Please see scorefunction options in ligand/ptm/ncaa manager"
               three_letter_map = self.ncaa_manager.three_letter_to_prop_map
               for three_letter_code in self.off_by_default_aa:
                    try:
                         self.ncaa_manager.set_prop(three_letter_map[three_letter_code])
                         self.ncaa_manager.enable()
                    except KeyError:
                         pass
               
          self.input_class.load_pose(self.cleaned_pdb_path.get())
     
     def set_ignore_unrecognized(self):
          """
          Uses the option system manager to extend options in order to ignore unrecognized before load.
          """
          manager = OptionSystemManager()
          manager.common_option.set('-ignore_unrecognized_res')
          manager.extend()
     
     def runFixPDB(self):
          """
          Runs the fixpdbprotocol
          """
          print self.PDBlistPath.get()
          filename = self.PDBlistPath.get()
          """
          outdirectory = tkFileDialog.askdirectory(title = "Output DIR", initialdir = global_variables.current_directory)
          if not outdirectory: return
          global_variables.current_directory = outdirectory
          """
          if not filename:return
          outdirectory = os.path.dirname(filename)
          outdirectory = outdirectory + "/CLEANED"
          if not os.path.exists(outdirectory):os.mkdir(outdirectory)
          print "Saving cleaned PDB(s) to "+outdirectory
          
          if re.search(".pdb", filename):
               self.fixonepdb(filename, outdirectory)
               self.cleaned_pdb_path.set(outdirectory+'/'+os.path.basename(filename))
               self.loadbutton.config(state=NORMAL)
               
          elif re.search(".txt", filename):
               print "Fixing All PDBs in list.."
               LIST = open(filename, 'r')
               for pdbpath in LIST:
                    pdbpath = pdbpath.rstrip()
                    print pdbpath
                    self.fixonepdb(pdbpath, outdirectory, True)
               self.loadbutton.config(state=DISABLED)
               LIST.close()
               print "PDBs Fixed..."
          else:
               print "Please choose a .pdb or .txt file..."
               return
          
     def fixonepdb(self, filename, outdir, silence_dialogs=False):
          """
          Goes through each option, runs the specified command on a file name after loading it into PythonPDB. Saves PDB at the end.
          """
          self.clean_pdb = PythonPDB(filename)
          pdb_map = self.clean_pdb.get_pdb_map()
          
          
          if self.replace_res_and_atoms_var.get():
               print "Attempting to change residue names, atom names, and water"
               self.clean_pdb.clean_PDB()
               
          if self.remove_hetatm_var.get():
               print "Removing HETATM Residues"
               self.clean_pdb.remove_hetatm_atoms()
               
          if self.remove_waters_var.get():
               print "Removing Waters"
               self.clean_pdb.remove_waters()
               
          if self.remove_alternates_var.get():
               print "Fixing alternate residue insertions. Any chains with insertions are renumbered from one."
               self.clean_pdb.remove_alternate_residues()
               
          if self.change_occupancies_var.get():
               print "Changing all occupancies to 1"
               self.clean_pdb.change_occupancy()
          

          
          pdbname = os.path.basename(filename)
          
          self.clean_pdb.save_PDB(outdir+"/"+pdbname)
          #self.cleaned_pdb_name.set(outdir+"/"+pdbname)
          
          print "File Saved..."
          
          if not self.unrecognized_aa and not self.off_by_default_aa:
               print "PDB has no unrecognized aa, and all ncaa found are on by default!"
               return
          
          if self.check_rosetta_var.get():
               print "Checking for unrecognized residues"
               self.check_for_unrecognized_aa(self.clean_pdb.get_pdb_map())
               
               print "Checking for off-by-default residues"
               self.check_for_rosetta_off_aa(self.clean_pdb.get_pdb_map())
               
               FILE = open(outdir+"/"+pdbname.split(".")[0]+"_rosetta_check_log.txt", 'w')
               FILE.write("#unrecognized\n")
               unique = dict()
               for aa in self.unrecognized_aa: unique[aa]=0
               
               if unique:
                    print "unrecognized: "+repr(sorted(unique))
               for aa in sorted(unique):
                    FILE.write(aa+"\n")
               
               FILE.write("#off\n")  
               unique = dict()
               for aa in self.off_by_default_aa: unique[aa]=0
               
               if unique:
                    print "off-by-default: "+repr(sorted(unique))
               for aa in sorted(unique):
                    FILE.write(aa+"\n")
               FILE.close
               
               print "Please see "+outdir+"/"+pdbname.split(".")[0]+"_rosetta_check_log.txt"
          
          #Print Unrecognized and Rosetta Off Info.  if not 'silence' variable, pop up.
          #If not silence, ask to Load PDB?
     
     def check_for_unrecognized_aa(self, pdb_map):
          """
          Checks pdb_map for any residues that are unrecognized by any param/patch file.
          Appends self.unrecognized_aa
          """
          
          for num in pdb_map:
               found=False
               for aa in self.recognized_aa:
                    if aa==pdb_map[num]["three_letter_code"]:
                         found=True
                         break
               if not found:        
                    self.unrecognized_aa.append(pdb_map[num]["three_letter_code"])
                    
     
     def check_for_rosetta_off_aa(self, pdb_map):
          """
          Checks pdb_map for any residues that are off by default in Rosetta.
          It COULD autoload them?
          Appends self.off_by_default_aa
          """
          for num in pdb_map:
               found=False
               for aa in self.on_by_default_aa:
                    if aa==pdb_map[num]["three_letter_code"]:
                         found=True
                         break
               if not found:
                    self.off_by_default_aa.append(pdb_map[num]["three_letter_code"])
     
     def getfile(self):
          """
          simply gets a filename...Nessessary for Tkinter Unfortunately.
          """
          filepath = tkFileDialog.askopenfilename(initialdir = self.defaultdir)
          global_variables.current_directory = os.path.dirname(filepath)
          return filepath
