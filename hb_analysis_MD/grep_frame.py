import mdtraj as md
import time
import os

t0 = time.clock()
traj_index = 80
os.system('mkdir /archive/yzhang/PS2/quinone/md_traj/traj_%s' % traj_index)
topo = '/archive/yzhang/PS2/quinone/md_traj/step5_charmm2omm_keep.psf'
traj = '/archive/yzhang/PS2/quinone/md_traj/step7_%s.dcd' % traj_index

trajectory = md.load_dcd(traj, top=topo)

for i in range(0, 100, 2):
	trajectory[i].save_pdb('/archive/yzhang/PS2/quinone/md_traj/traj_%s/frame%s_%s.pdb' % (traj_index, i, traj_index), force_overwrite=True, bfactors=None)
t3 = time.clock()
print('total time:', t3)
