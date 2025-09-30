import numpy as np
import pyvista as pv
import sys

fname = sys.argv[1]

data= np.loadtxt(fname,skiprows=1)
positions = data[:,0:3]
radius = data[:,-1]

print("\n Summary:")
print(f" N = {radius.size}")
cloud = pv.PolyData(positions)

print("\n Saving vtp file for paraview rendering: particles.vtp")
cloud.save("particles.vtp")   

print("\n Writing lethe list insertion file ")
lethe_data=np.zeros((radius.size,10))
lethe_data[:,0:3]=positions
lethe_data[:,-1]=radius*2
header = "p_x; p_y; p_z; v_x; v_y; v_z; w_x; w_y; w_z; diameters;   T;"
np.savetxt("insertion_file.dat", lethe_data, header=header)
