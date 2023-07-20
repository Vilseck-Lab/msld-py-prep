
import os
import re
import msld_wrt

# callback for the "Align" button
def MsldALN(refpdb, ligcsv, savedir = '.', writeSDF = False, pattern = None, verbosity = 5)->list():

    """
    
    """
    
    from rdkit import Chem
    target = Chem.MolFromPDBFile(refpdb, removeHs = False)



    ## TODO
    # Verify alignment target is specified 
    """
    if not (target or pattern):
        show_warning('No Molecule Selection!','Please specify the\
                    molecule you wish to align to.')
        return None

    if pattern and not (dataset or cmd.get_names('objects')):
        show_warning('No Molecules to Align!', 'Please load or specify\
                    the molecule(s) you wish to align.')
        return None
    """

    # Generate 3D conformers and align using RDKit if SMILES csv is
    # specified
    if target is not None:
        refpdb_bname = re.sub(r'.pdb$','',os.path.basename(refpdb))
        
        if not os.path.exists(os.path.join(savedir,f"{refpdb_bname}.sdf")):
            # cmd.save(os.path.join(savedir,f"{target}.sdf"),selection=target,format='sdf')
            fsdf = Chem.SDWriter(os.path.join(savedir, f"{refpdb_bname}.sdf"))
            fsdf.write(target)
            fsdf.close()

        generate3DMols(fn = ligcsv,wd=savedir,\
            target=os.path.join(savedir,f"{refpdb_bname}.sdf"), \
            writeSDF = writeSDF)


    else:
        pass

    return 
#END


def generate3DMols(fn, wd, smiles_col=1, names_col=0, delim=',',pattern=None,target=None, writeSDF = False):
    """
    3D conformer generation helper functions

    Parameters
    ----------
    fn: 
        CSV file
    wd:
        directory where files will be written to.

    smiles_col: int, default = 1
        0-indexed column ID for the column which has the ligand SMILES string in the csv file `fn`

    names_col: int, default = 0

    delim: str, default=`,`

    pattern: str, default = None

    target: str, default = NOne
        path to the reference molecule's SDF file

    writeSDF: bool, default = False
        True <==> write the SDF file of all the ligands into a SDF file.
    
    """

    import generate_3D
    import rename_atoms

    mols = generate_3D.RDKit_Tools(fn,smiles_col=smiles_col,names_col=names_col,delim=delim)
    print(mols.df)
    print('SMILES list: ', mols.smiles_list)
    
    generated_mols = mols.generate_3D(ref=target)
    print(f"<generate3DMols> MOL files for the molecules have been generated and saved in {os.path.abspath(wd)}")
    print(generated_mols)

    

    print('<generate3DMols> Writing and refining MOL2 files for all the substrates .. ')
    for molname in generated_mols:

        # convert MOL files into MOL2 files using openbabel
        molFile = os.path.join(wd, f"{molname}.mol")
        mol2File = os.path.join(wd, f"{molname}.mol2")
        msld_wrt.openbabel_convert(molFile, 'mol', mol2File, 'mol2')
        
        # Change atom naming to 4 characters
        truncate_at_names(mol2File)


        # Uniquely rename atoms (required for PyPrep Scripts)
        rename_atoms.rename_atoms(mol2File) 

        if writeSDF:
            sdfFile = os.path.join(wd, f"{molname}.sdf")
            msld_wrt.openbabel_convert(mol2File, 'mol2', \
                sdfFile, 'sdf')
    return 

        
            
#END



def truncate_at_names(arg_mol2File)->None:
    """
    PyMOL outputs atom names with greater than 4 characters for
    halogens. Need to truncate for CHARMM.

    Notes
    -----
    Overwrites the input `arg_mol2File`.
    """

    import string 
    from copy import deepcopy
    
    with open(arg_mol2File,'r') as f:
        lines = f.readlines()
    
    newlines = deepcopy(lines)
    idx = [i for i in range(len(lines)) if lines[i].startswith('@<TRIPOS>ATOM') or lines[i].startswith('@<TRIPOS>BOND')]
    
    atnames = [l.rstrip().split()[1] for l in lines[idx[0]+1:idx[1]]]
    atnames = [at for at in atnames]
    alphabet = list(string.ascii_uppercase)
    for i,line in enumerate(lines[idx[0]+1:idx[1]]):
        if any(patt in line for patt in ['Cl','CL','Br','BR']):
            atname = line.rstrip().split()[1]
            if len(atname) <= 4:
                continue
            newatname = atname[0:3]+atname[-1]
            atnames[i] == newatname
            it=0
            while newatname in atnames:
                newatname = atname[0:3]+alphabet[it]
                it+=1
            newline = line.replace(atname, newatname)
            newlines[i+idx[0]+1] = newline
        
    with open(arg_mol2File,'w') as f:
        f.write("".join(newlines))

    return

        