import MDAnalysis as mda
import sys
import numpy as np
from scipy.spatial import distance
from matplotlib.image import NonUniformImage
import matplotlib.pyplot as plt
from MDAnalysis.lib.distances import distance_array, calc_angles, calc_bonds

np.set_printoptions(threshold=sys.maxsize)

import math

def read_xyz(fin):
    """ read a xyz file from file handle
    Parameters
    ----------
    fin : file handle
        file to read from
    Returns
    -------
    fin : open file
    xyz : namedtuple
        returns a named tuple with coords, title and list of atomtypes.
    See Also
    --------
    write_xyz
    """
    
    natoms = int(fin.readline())
    title = fin.readline()[:-1]
    coords = np.zeros([natoms, 3], dtype="float64")
    atomtypes = []
    for x in coords:
        line = fin.readline().split()
        atomtypes.append(line[0])
        x[:] = list(map(float, line[1:4]))

    return atomtypes


path="./gyration_traj.xyz"
traj="./simulation.xc.xyz"

u = mda.Universe(path,guess_bonds=True,vdwradii={'O':1.5, 'H':1.5, 'M':1}, dt=0.0005,box=[11.295,11.295,11.295, 90, 90, 90])
u.dimensions=[11.295,11.295,11.295, 90, 90, 90]
mda.core.flags['use_periodic_selections'] = True
a = mda.Universe(traj,guess_bonds=True,vdwradii={'O':1.5, 'H':1.5, 'D':1.5}, dt=0.0005,box=[11.295,11.295,11.295, 90, 90, 90])

d2o=np.zeros(u.trajectory.n_frames)
h2o=np.zeros(u.trajectory.n_frames)
hdo=np.zeros(u.trajectory.n_frames)
h_hdo=np.zeros(u.trajectory.n_frames)
d_hdo=np.zeros(u.trajectory.n_frames)
tmp_h=np.zeros(u.trajectory.n_frames)
tmp_d=np.zeros(u.trajectory.n_frames)
h_h2o=np.zeros(u.trajectory.n_frames)
d_d2o=np.zeros(u.trajectory.n_frames)
count_hd=np.zeros(u.trajectory.n_frames)



fin=open("simulation.xc.xyz")

for ts in u.trajectory:
    h_shell=u.select_atoms("name H and around 2.3 name Mn",updating=True)
    o_shell=u.select_atoms("name O and around 1.3 (name H and around 2.3 name Mn)",updating=True)
    for j in range(20):
        atoms=read_xyz(fin)
        name_atoms=np.asarray(atoms)
        
        H1=name_atoms[o_shell.atoms.ids]
        H2=name_atoms[o_shell.atoms.ids+1]
        
        HD=name_atoms[h_shell.atoms.ids-1]
        tmp_h[ts.frame]=tmp_h[ts.frame]+list(HD).count("H")
        tmp_d[ts.frame]=tmp_d[ts.frame]+list(HD).count("D")
        count_hd[ts.frame]=count_hd[ts.frame]+len(HD)

        #print(j)
        for i in range(len(o_shell.atoms.ids)):
            if ( H1[i] == "H" and H2[i] == "H" ):   #H2O
                h2o[ts.frame]=h2o[ts.frame]+1
                if o_shell.atoms.ids[i]+1 in h_shell.atoms.ids :
                    if H1[i] == "H" :
                        h_h2o[ts.frame] = h_h2o[ts.frame] + 1
                if o_shell.atoms.ids[i]+2 in h_shell.atoms.ids :
                    if H2[i] == "H" :
                        h_h2o[ts.frame] = h_h2o[ts.frame] + 1                
                
            if ( H1[i] == "D" and H2[i] == "D" ):   #D2O
                d2o[ts.frame]=d2o[ts.frame]+1
                if o_shell.atoms.ids[i]+1 in h_shell.atoms.ids :
                    if H1[i] == "D" :
                        d_d2o[ts.frame] = d_d2o[ts.frame] + 1
                if o_shell.atoms.ids[i]+2 in h_shell.atoms.ids :
                    if H2[i] == "D" :
                        d_d2o[ts.frame] = d_d2o[ts.frame] + 1                   
                
                
                
            if ( H1[i] != H2[i] ):                                   #HDO
                hdo[ts.frame]=hdo[ts.frame]+1
                if o_shell.atoms.ids[i]+1 in h_shell.atoms.ids :
                    #print("H1",H1[i],o_shell.atoms.ids[i]+1,h_shell.atoms.ids)
                    if H1[i] == "H" :
                        h_hdo[ts.frame] = h_hdo[ts.frame] + 1
                    if H1[i] == "D":
                        d_hdo[ts.frame] = d_hdo[ts.frame] + 1
                if o_shell.atoms.ids[i]+2 in h_shell.atoms.ids :
                    #print("H2",H2[i],o_shell.atoms.ids[i]+2,h_shell.atoms.ids)
                    if H2[i] == "H" :
                        h_hdo[ts.frame] = h_hdo[ts.frame] + 1
                    if H2[i] == "D":
                        d_hdo[ts.frame] = d_hdo[ts.frame] + 1

print("#H2O", np.sum(h2o[300:]),"#D2O:",np.sum(d2o[300:]),"#DHO",np.sum(hdo[300:]))
print("#H2O/D2O",np.sum(h2o[300:])/np.sum(d2o[300:]),"#H/D", np.sum(tmp_h[300:])/np.sum(tmp_d[300:]),"#H/D(HDO)",np.sum(h_hdo[300:])/np.sum(d_hdo[300:]))
print("#H_total, H(H2O), H(HDO)",np.sum(tmp_h[300:]),np.sum(h_h2o[300:]),np.sum(h_hdo[300:]))
print("#D_total, D(D2O), D(HDO)",np.sum(tmp_d[300:]),np.sum(d_d2o[300:]),np.sum(d_hdo[300:]))
