import pyscal as pc
import os
import pyscal.traj_process as ptp
import matplotlib.pyplot as plt
import numpy as np
import string
from scipy.io import savemat
import pickle

#qs=np.array(pickle.load(open('qs.pkl','rb')))
#qs=np.array(pickle.load(open('qssd.pkl','rb')))
framenum=180
#atomnum=qs.shape[2]
histnum=100#直方图分割数量
file='postest'
format="poscar"
sys1 = pc.System()
sys1.read_inputfile(file, format=format)
sys1.find_neighbors(method="cutoff", cutoff=0)
sys1.calculate_q([2,4,5,6],averaged=True)
r,hist,m=sys1.calculate_rdf()
