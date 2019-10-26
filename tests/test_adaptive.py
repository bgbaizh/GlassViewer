import pytest
import os
import numpy as np
import pyscal.core as pc
import pyscal.crystal_structures as pcs

def test_neighbors_sann():
    #create some atoms
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [6, 6, 6])
    sys = pc.System()
    sys.atoms = atoms
    sys.box = boxdims


    #then lets find neighbors
    #SANN algo test
    sys.find_neighbors(method = 'cutoff', cutoff='sann')
    #any atom should have 8 neighbors
    atoms = sys.atoms
    assert atoms[0].coordination == 14

    sys.find_neighbors(method = 'cutoff', cutoff='sann', threshold=1)
    #any atom should have 8 neighbors
    atoms = sys.atoms
    assert atoms[0].coordination == 14

    sys.calculate_q(8, averaged=True)
    q = sys.get_qvals(8, averaged=True)
    assert np.round(np.mean(np.array(q), dtype=np.float64), decimals=2) == 0.45


    #sys.calculate_q(7, averaged=True)
    #q = sys.get_qvals(7, averaged=True)
    #summ = 0
    #counts = len(q)
    #for qq in q:
    #    summ += qq
    #avgg = summ/float(counts)


    #assert avgg == 0.05



def test_neighbors_adaptive():
    #create some atoms
    atoms, boxdims = pcs.make_crystal('bcc', repetitions = [6, 6, 6])
    sys = pc.System()
    sys.atoms = atoms
    sys.box = boxdims

    #then lets find neighbors
    #SANN algo test
    sys.find_neighbors(method = 'cutoff', cutoff='adaptive')
    #any atom should have 8 neighbors
    atoms = sys.atoms
    assert atoms[0].coordination == 14
