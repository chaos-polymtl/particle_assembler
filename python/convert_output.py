import numpy as np
import pyvista as pv
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='Convert particle assembler output to VTP format for ParaView visualization')
parser.add_argument('-i','--input', type=str, help='Input file containing particle positions and radii (format: x y z radius)')
parser.add_argument('--output-vtp', type=str, default='particles.vtp', help='Output VTP file name (default: particles.vtp)')
parser.add_argument('--output-lethe', type=str, default='insertion_file.dat', help='Output Lethe insertion file name (default: insertion_file.dat)')

args = parser.parse_args()

# Load data
data = np.loadtxt(args.input, skiprows=1)
positions = data[:, 0:3]
radius = data[:, -1]

print("\n Summary:")
print(f" N = {radius.size}")
cloud = pv.PolyData(positions)

print(f"\n Saving vtp file for paraview rendering: {args.output_vtp}")
cloud.save(args.output_vtp)   

print(f"\n Writing lethe list insertion file: {args.output_lethe}")
lethe_data = np.zeros((radius.size, 10))
lethe_data[:, 0:3] = positions
lethe_data[:, -1] = radius * 2
header = "p_x; p_y; p_z; v_x; v_y; v_z; w_x; w_y; w_z; diameters;   T;"
np.savetxt(args.output_lethe, lethe_data, header=header,delimiter=" ; ")
