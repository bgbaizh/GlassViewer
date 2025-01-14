import numpy as np
import gzip
import glassviewer.catom as pca
from ase import Atom, Atoms
import gzip
import io
import os

#new function to wrap over mdtraj objects
def read_snap(mdobject, check_triclinic=False):
    """
    Function to read from an MDTraj atoms objects

    Parameters
    ----------
    mdobject : MDTraj Atoms object
        name of the MDTraj atoms object

    triclinic : bool, optional
        True if the configuration is triclinic

    """
    #We have to process atoms and atomic objects from ase
    #Known issues lammps -dump modified format
    #first get box
    a = np.array(mdobject.unitcell_vectors[0][0])
    b = np.array(mdobject.unitcell_vectors[0][1])
    c = np.array(mdobject.unitcell_vectors[0][2])

    box = np.array([a, b, c])

    #box and box dims are set. Now handle atoms
    chems = np.array([atom.name for atom in mdobject.topology.atoms])
    atomsymbols = np.unique(chems)
    atomtypes = np.array(range(1, len(atomsymbols)+1))
    typedict = dict(zip(atomsymbols, atomtypes))

    #now start parsing atoms
    atoms = []
    positions = mdobject.xyz[0]
    for count, position in enumerate(positions):
        atom = pca.Atom()
        atom.pos = list(position)
        atom.id = count+1
        atom.type = typedict[chems[count]]
        atom.loc = count

        customdict = {'species': chems[count]}
        atom.custom = customdict
        atoms.append(atom)

    return atoms, box

def write_snap(**kwargs):
    raise NotImplementedError("write method for mdtraj is not implemented")

def split_snaps(**kwargs):
    raise NotImplementedError("split method for mdtraj is not implemented")

def convert_snap(**kwargs):
    raise NotImplementedError("convert method for mdtraj is not implemented")
