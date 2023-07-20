#! /usr/bin/env python

##
## Check molfile and toppar files before getting started
##

import numpy as np
import glob
import os
import re
from pandas import read_csv

class CHK_Error(Exception):
    import sys
    sys.exit

def MsldCHK(ligcsv, odir, verbosity = 5, names_col = 0, smiles_col = 1, **kwargs):
    """
    Check the following points in each mol2file:
     - each atom is uniquely named (within a single mol2 file)
     - each mol2 file has an associated rtf file
       (parse ParamChem stream files into rtf/prm files)

    Parameters
    ----------
    ligcsv: []
        csv file with ligand names and their SMILES strings.

    odir: str
        directory where ligand MOL2, RTF, PRM and STR files are expected to be.

    verbosity: int, default=5
        level of verbosity

    """

    ligDF = read_csv(ligcsv)
    mol2flist = [os.path.join(odir,f'{name}.mol2') for name in ligDF.iloc[:,names_col]]

    qrename = False
    print(kwargs)
    try:
        if kwargs['rename_atoms']:
            import sanitize_mol2
            qrename = True
    except KeyError:
        print('Not renaming atoms in the MOL2 file.')
        qrename = False
        pass


    # Some quick file checks:
    for mdx,fmol2 in enumerate(mol2flist):
        
        molbase = re.sub(r'\.mol2$', '', os.path.basename(fmol2))

        rtfname = os.path.join(os.path.dirname(fmol2), f'{molbase}.rtf')
        prmname = os.path.join(os.path.dirname(fmol2), f'{molbase}.prm')
        strname = os.path.join(os.path.dirname(fmol2), f'{molbase}.str')

        # Check to make sure the mol2 files exist
        # if len(glob.glob(mols[mol]+'.mol2')) == 0:
        if not os.path.isfile(fmol2):
            raise CHK_Error(f'{fmol2} does not exist - check and resubmit')


        if qrename:   
            # Uniquely rename atoms (required for PyPrep Scripts)
            sanitize_mol2.rename_atoms_in_mol2(arg_infile = fmol2) 


        # Check to see if every atom is uniquely named - exit if not
        chk=0
        fp=open(fmol2,'r')
        line=fp.readline()
        aname=[]
        while line:
            if line[0:13] == '@<TRIPOS>ATOM':
                line=fp.readline()
                while line[0:13] != '@<TRIPOS>BOND':
                    tmp=line.split()[1]
                    if not (tmp in aname):
                        aname.append(tmp)
                    else:
                        chk=1
                    line=fp.readline()
                break
            else:
                line=fp.readline()
        fp.close()
        if chk:
            raise CHK_Error(f'{molbase}.mol2 atom names are NOT unique - fix this in both structure and toppar files - then resubmit')

        ## If atom names are not unique - you must rename them before

        # Check to make sure rtf/prm files exist (parse str files)
        # if len(glob.glob(mols[mol]+'.rtf')) == 0:
        if not os.path.exists(rtfname):
            if verbosity >= 3:
                print(f'<msld_chk> Expected RTF file for molecule {mdx} [{rtfname}] not found. Looking for STR file as a last ditched effort.')
            # check to see if a ParamChem stream file exists:
            # if len(glob.glob(mols[mol]+'.str')) == 0:
            if not os.path.exists(strname):
                raise CHK_Error(f'<msld_chk> Expected STR file for molecule {mdx} [{strname}] not found - check and resubmit')
            
            else:
                print(f'<msld_chk> Found STR file for molecule {mdx} [{strname}].')
                print(f'<msld_chk> Writing RTF and PRM files from STR file for molecule {mdx}...')
                # assume the stream file is from paramchem
                fp = open(strname, 'r') #fp=open(mols[mol]+'.str','r')
                rp = open(rtfname, 'w') #rp=open(mols[mol]+'.rtf','w')
                pp = open(prmname, 'w') #pp=open(mols[mol]+'.prm','w')

                line=fp.readline()
                while line:
                    # rtf first
                    if line[0:8] == 'read rtf':
                        line=fp.readline()
                        while line[0:3] != 'END':
                            rp.write(line)
                            line=fp.readline()
                        rp.write(line)
                        line=fp.readline()
                    # prm next
                    elif line[0:10] == 'read param':
                        line=fp.readline()
                        while line[0:3] != 'END':
                            pp.write(line)
                            line=fp.readline()
                        pp.write(line)
                        line=fp.readline()
                    else:
                        line=fp.readline()

                fp.close()
                rp.close()
                pp.close()

        # Check that atom names in the mol2 file match the atom names in the rtf
        print(f'<msld_chk> Checking for unmatched atom names in the MOL2 and RTF of molecule {mdx}...')
        a2name=[]
        fp= open(rtfname, 'r')  #open(mols[mol]+'.rtf','r')
        for line in fp:
            if line[0:4] == 'ATOM':
                a2name.append(line.split()[1])
        fp.close()

        cnt=0
        for at1 in aname:
            for at2 in a2name:
                if at1 == at2:
                    cnt+=1
        if cnt != len(aname):
            raise CHK_Error(f'Atom names in {fmol2} do not match those found in {rtfname}! Check and resubmit.')

        print(f'<msld_chk> Check for molecule {mdx} passed.\n\n')
        # mols[mol] ready to go

    return

