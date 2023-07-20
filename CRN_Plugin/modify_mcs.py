
def removeFromCore(mcsout, tosite, pymolspec, out=None):
    """

    In:
    mcsout    : (str) MCS file name
    tosite    : (int) site to move specified atom to
    pymolspec : (str) pymol selection
    out       : (str) if specified, new MCS file will be
                written. If None, mcsout will be modified
                in place

    Out:
    newmcsout : Modified mcsout
 
    """
    if out:
        newmcsout = out
    else:
        newmcsout = mcsout

    ind = [x.index(y) for y in x if y.startswith('CORE') or y.startswith('ANCHOR']
    mol_names = [line.split()[0] for line in x[indices[0]+1:indices[1]] if line.split()]
    core_atoms = [line.split()[1:] for line in x[indices[0]+1:indices[1]] if line.split()] 
 
    print(mol_names)
    print(core_atoms)

