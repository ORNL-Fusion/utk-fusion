# This script takes in an STL mesh of the wall geometry and outputs a file of a
# cross section of it at every degree formatted for MAFOT. This file is then
# passed into MAFOT via the -W flag.
#
# Example usage: python3 generate_3d_wall_for_mafot path_to_stl_file [num_coordinates]
#   - path_to_stl_file is the full path to the STL mesh file.
#   - num_coordinates is how many coordinates you would like each cross section
#      to have.
import sys
import numpy as np
import trimesh
from tqdm import tqdm

def main(stl_path, num_coords, angle_step=1.0):
    """
    Main routine to generate coordinates of each toroidal slice of a supplied
    3D wall, formatted for input to MAFOT via the -W flag.

    stl_path (str): Path to the STL file.
    num_coords (int): Number of coordinates for each toroidal slice of the wall.
    angle_step (float): Currently only 1.0 is supported due to it being
      hard-coded into MAFOT, but this may change in the future.
    """

    # First load in the STL file with a soft exit if failed.
    try:
        mesh = trimesh.load_mesh(stl_path)
        print("Loaded mesh file: {}".format(stl_path.split("/")[-1]))
    except ValueError:
        print("Error! STL file not found: {}".format(stl_path))
        return None

    # Require watertightness.
    if not mesh.is_watertight:
        print("Error! Mesh must be watertight.")
        return None

    # We want to take slices of the 3D volume to get a successive sequence of
    # 2D cross-sections, one at each degree. We do this by giving trimesh the
    # surface normal of the intersecting plane.
    sections = []
    degrees = np.arange(0.0, 360.0, angle_step)
    centroid = mesh.centroid
    print("Generating cross sections...")
    for i in tqdm(range(0, len(degrees))):
        deg = degrees[i]

        # We add 90 since we want the surface normal at this angle (which is
        # 90 degrees away).
        x_norm = np.cos(np.radians(deg + 90)) + centroid[0]
        y_norm = np.sin(np.radians(deg + 90)) + centroid[1]

        # Grab the 2D cross-section, which actually will include two cross-
        # sections of the vessel. We want the one on the right. This is to say
        # if phi = [0, 180), then we want to keep the vertices where x < 0,
        # and if phi = [180, 360) where x > 0.
        #slice = mesh.section(plane_origin=[0,0,0],
        #  plane_normal=[x_norm, y_norm, 0])
        slice = mesh.section(plane_origin=centroid,
          plane_normal=[x_norm, y_norm, 0])
        slice_2D, to_3D = slice.to_planar()
        verts = np.array(slice_2D.vertices.tolist())

        if deg >= 0 and deg < 180:
            keep = verts[:, 1] >= 0
            x_cs = verts[:, 0][keep]  # cs = cross-section
            y_cs = verts[:, 1][keep]

        # The y's are actually R, and the x's are actually Z in machine
        # coordinates, and the machine angle is actually 360 - deg since it
        # goes clockwise.
        r = y_cs / 1000
        z = x_cs / 1000
        phi_mach = 360 - deg

        sections.append({"degree":deg, "x_norm":x_norm, "y_norm":y_norm,
          "x_cs":x_cs, "y_cs":y_cs, "r":r, "z":z, "phi_mach":phi_mach,
          "slice_2D":slice_2D})

    # Optional flag to make a gif of the wall as you progress around the vessel.
    if make_gif:
        import imageio
        import os



if __name__ == "__main__":

    print()
    print("*****************************************")
    print("* Generate MAFOT input file for 3D wall *")
    print("*****************************************")

    if len(sys.argv) in [2, 3]:
        stl_path = sys.argv[1]
        if len(sys.argv) == 2:
            num_coords = 100
        else:
            num_coords = int(sys.argv[2])
        main(stl_path, num_coords)

    else:
        print("Error: Only pass in the number of coordinates for each" + \
          " toroidal\nslice of the wall.")
        print("Example: python3 generate_3d_wall_for_mafot path_to_stl.stl 150")
