#! /usr/bin/env python

####
#### Executable script to build MSLD ready ligand files
#### JV,LC 01/2022
####

import msld_chk
import msld_mcs
import msld_crn
import msld_prm
import msld_wrt
import msld_aln
# import glob
import sys
import os
import optparse

###
### This script is executed in 2 steps to (1st) build the MCSS
### and (2nd) perform charge renormalization. This allows the user
### to check that the identified MCSS is correct with vis_check.py.
### Thus, the user should manually call msld_py_prep.py twice.
###
### All ligand structure files (mol2) and toppar files must be 
### available prior to running this script
###
### The inFrag list of lists allows you to move core atoms into 
### alchemical fragments at specific sites. Each nested list 
### corresponds to a single site attached to the ligand core.
###
### The AnCore list of lists allows you to move (non-DUM) atoms
### listed as "anchor atoms" (atoms connecting core and fragment
### components) into the core upon charge renormalization.
###


#####################################################################
## (1) Define System and File Variables
##
## (1) Read/gather inputs
inp = optparse.OptionParser(prog = "msld_py_prep", usage = "Create files to run MSLD and ALF.")
inp.add_option('-n', dest = 'sysname', default = 'name', help = 'name of future output files.', type = str)
inp.add_option('-m', dest = 'mcsout', default = 'MCS_for_MSLD.txt', help = 'Maximum common substructure output filename', type = str)
inp.add_option('--overwrite', dest = 'overwrite_mcsout', action = 'store_true', default = False, help = 'Overwrite the MCS file.')
inp.add_option('-r','--refpdb', dest = 'refstruct', default = None, help = '3D structure of the reference ligand (in PDB only for now)', type = str)
inp.add_option('-c','--ligcsv', dest = 'ligcsv', default = None, help = 'A CSV file consisting of two columns. (Molecule name), (SMILES)', type = str)
inp.add_option('-V','--prnlev', dest = 'prnlev', default = 5, help = 'Level of verbosity.', type = int)
inp.add_option('--wdir', dest = 'wdir', default = '.', help = 'Directory where all the ligands are located (may be required if -o is provided.)')
inp.add_option('--odir', dest = 'outdir', default = None, help = 'MSLD output directory', type = str)
inp.add_option('--cgenff', dest = 'qcgenff', action = 'store_true', default = True, help = 'Are CGenFF/ParamChem parameters being used?')
inp.add_option('--cutoff', dest = 'cutoff', default = 0.8, help = 'RMSD & distance cutoff to differentiate different atoms', type = float)
inp.add_option('--align', dest = 'align', action = 'store_true', default = False, help = 'Align ligands to a reference ligand structure (also writes MOL2 files).')
inp.add_option('--sdf', dest = 'writeSDF', action = 'store_true', default = False, help = 'Write a SDF file for all the ligands from their MOL2 files.')
inp.add_option('--dqcheck', dest = 'ChkQChange', action = 'store_true', default = True, help = 'Option for msld_CRN.')

if len(sys.argv[1:]) == 0:
    print('Incomplete inputs.')
    sys.exit()

opts, args = inp.parse_args()

def die(arg_msg):
    sys.stderr.write('\n\n' + arg_msg +'\n')
    quit()

#####################################################################
# input curation
sysname = opts.sysname
wdir = opts.wdir

if opts.outdir is None:
    outdir = f'build.{sysname}'
else:
    outdir = opts.outdir

cgenff=opts.qcgenff
prnlev = opts.prnlev



if opts.align:
    ## (1.a) ALign molecules
    msld_aln.MsldALN(refpdb = opts.refstruct, \
        ligcsv = opts.ligcsv, \
        savedir = opts.wdir, \
        writeSDF = True, \
        pattern = None, \
        verbosity = prnlev)
    
    die('<msld_py_prep> (REMOVE) Quitting after aligning!!!.')


### MAIN ###
inFrag=[[]]  # reflig core atoms to include in each fragment at each site (list of nsub lists)
AnCore=[[]]  # anchor atoms at each site to include in the core (list of nsub lists)


#####################################################################
if (not os.path.exists(opts.mcsout)) or (opts.overwrite_mcsout):

    
    ## (2) Check molfile and toppar files before getting started
    print('<msld_py_prep> Running some checks before identifying the Maximum Common Substructure (MCS)...')
    msld_chk.MsldCHK(ligcsv = opts.ligcsv, \
        odir = opts.wdir, \
        verbosity = prnlev, **{'rename_atoms' : False})
    print("\n<msld_py_prep> chk finished\n\n")

    
    #####################################################################
    ## (3) Maximum Common SubStruct Search with bonded-environments
    ## "mcsout" = results of the search; edit this file to manual edit ligand splicing
    ## cutoff = RMSD & distance cutoff to differentiate different atoms
    ## change debug to True to get more stdout printed 
    
    # reflig = msld_mcs.MsldMCS(molfile,mcsout,cutoff= opts.cutoff,debug=False)
    reflig = msld_mcs.MsldMCS(ligcsv = opts.ligcsv, \
        odir = opts.wdir, \
        mcsout = opts.mcsout, \
        names_col = 0, \
        cutoff= opts.cutoff, \
        debug=True, \
        verbosity = prnlev)

    print("\n\n")
    print(f"<msld_py_prep> MCS results printed to : {opts.mcsout}")
    print(f"<msld_py_prep> Reference Ligand is : {reflig}")

    die("<msld_py_prep> (REMOVE) Quitting after rendering MCS.")


#####################################################################
## (4) Perform Charge-Renormalization 
## To manually move atoms from the core into alchemical fragments, 
## update the "inFrag" variable above. "Anchor atoms" (and connected Hs)
## are automatically included in each alchemical fragment unless 
## specifically stated to be "AnCore"

msld_crn.MsldCRN('MCS_for_MSLD.yml',\
    outdir,\
    inFrag,\
    AnCore,\
    ChkQChange=opts.ChkQChange,\
    verbosity=prnlev,\
    debug=True)


#####################################################################
## (5) Write Ligand Parameters & the Charmm ALF input scripts
# msld_prm.MsldPRM(outdir,cgenff,verbose=False,debug=False)
msld_prm.MsldPRM(mcsout = 'MCS_for_MSLD.yml', \
    wdir = wdir,\
    outdir = outdir, \
    cgenff = cgenff, \
    verbosity = prnlev, \
    debug=True)


# msld_wrt.writeALF_Files(sysname,outdir,cgenff, yml = True)
msld_wrt.writeALF_Files(sysname, mcsout = 'MCS_for_MSLD.yml', \
    outdir = outdir, \
    cgenff = cgenff, \
    wd = wdir, \
    charmm=True, 
    yml = True)

"""
## Final Notes to the user
print(f"default TOPPAR parameters copied into {outdir}. Check to make sure these work for your system!")


## FINISHED
"""
