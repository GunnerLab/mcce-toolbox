"""
Script to run MCCE simulation at different charges for water molecules.
"""
import os
import sys
import numpy as np
from scipy import stats
from pymcce.automated_mcce import MCCEParams
from pymcce.mcce_simulation import Simulation
from pymcce.utils import write_watpdb_from_coords, get_last_prot_at_index

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import pylab
import seaborn as sns
from sklearn.neighbors import KernelDensity
from numpy import *
import pandas as pd
sns.set(style="white")




mu = []
ab_indices = []
n_wat = []
dipole_x = []
dipole_y = []
dipole_z = []
data_dir = '/home/yzhang/Dropbox/ProtonHopping/data/gramicidin/simulations/input_struct_dp_groups_positive_t3p'
prefix = "run_restart200_000001_1_update" 
print("Processing %s:" % prefix)
msdat = os.path.join(data_dir, prefix, "ms.dat")
head3lst = os.path.join(data_dir, prefix, "head3.lst")
fort38 = os.path.join(data_dir, prefix, "fort.38")
step2out = os.path.join(data_dir, prefix, "step2_out.pdb")

msa = Simulation(msdat, head3lst, fort38)
msa.parse_trajectory(sample_frequency=10)
msa.parse_struct(step2out)
conf_dipoles = msa.calculate_dipoles()
print(conf_dipoles['HOH01W0233_006'])
numbers = zeros(msa.trajectory.shape[0], dtype="int64")
for i in range(msa.trajectory.shape[0]):
    microstate_conf_ids = msa.trajectory[i, :]
    numbers[i] = msa.state_counts[i]
    dps = []
    curr_wat_ids = []
    for index, c in enumerate(microstate_conf_ids):
        conf_name = msa.conf_id_name_map[c + 1]
        if "DM" not in conf_name:
            dpX = conf_dipoles[conf_name][0]
            dpY = conf_dipoles[conf_name][1]
            dpZ = conf_dipoles[conf_name][2]
            dps.append([dpX, dpY, dpZ])
    dps = np.array(dps)
    x = sum(dps[:, 0])
    dipole_x.append(x)
    y = sum(dps[:, 1])
    dipole_y.append(y)
    z = sum(dps[:, 2])
    dipole_z.append(z)

dipole_ms = pd.DataFrame({'x': dipole_x, 'y': dipole_y, 'z': dipole_z, 'count': numbers})
dipole_ms.to_csv('dipole_microstate.csv')
