# This script takes in an STL mesh of the wall geometry and outputs a file of a
# cross section of it at every degree formatted for MAFOT. This file is then
# passed into MAFOT via the -W flag.
#
# Example usage: python3 generate_3d_wall_for_mafot path_to_stl_file [make_gif]
#   - path_to_stl_file is the full path to the STL mesh file.
#   - make_gif is whether or not to produce a gif showing the wall at each angle
import sys
import math
import trimesh
import numpy as np
from tqdm import tqdm
import matplotlib.pyplot as plt
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


# Sort points poloidally.
def sort_polar(center, rs, zs, poly):

    # Pop the entries that fall within the polygon. These are the ones
    # we want to sort.
    points = [Point(rs[i], zs[i]) for i in range(0, len(rs))]
    inside = np.array([poly.contains(p) for p in points])

    # The Points object consists of an __array__ attribute, and numpy
    # does not like creating an array of arrays that may be different
    # shapes (i.e. np.array(points)), so the suggested fix is to
    # created an empty array of type object and then fill it.
    points_np = np.empty(len(points), dtype=object)
    for i in range(0, len(points)):
        points_np[i] = points[i]

    # Create a new array that we will put the sorted points into,
    # starting at first_idx.
    new_points = points_np[~inside]
    sort_points = points_np[inside]
    first_idx = np.argwhere(inside==True)[0]

    # Sort the points according to polar angle. We could in theory use
    # the centroid but for greater flexibility best to just allow the
    # user defined center coordinate. 2*pi - angle to go clockwise since that
    # just feels natural.
    pol_ang = [2 * np.pi - math.atan2(p.y - center[1], p.x - center[0]) for p in sort_points]
    sort_idx = np.argsort(pol_ang)

    # Insert into the new_points array at the first index where the
    # sorted points were pulled from.
    new_points = np.insert(new_points, first_idx, sort_points[sort_idx])

    # Break back down into r and z for output.
    rs = [p.x for p in new_points]
    zs = [p.y for p in new_points]

    return rs, zs

def main(stl_path, angle_step=1.0, make_gif=False):
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

    # Polygons defining regions to sort and re-sort points.
    polys = []
    polys.append(Polygon([(-99, -99),    (99,   -99),   (99,   99),    (-99,  99)]))
    polys.append(Polygon([(2.00, 0.35),  (2.60, 0.35),  (2.60, 1.10),  (2.00, 1.10)]))
    polys.append(Polygon([(1.80, -1.10), (2.60, -1.10), (2.60, -0.35), (1.80, -0.35)]))
    polys.append(Polygon([(1.45, 1.08),  (1.45, 1.30),  (1.55, 1.30),  (1.55, 1.08)]))
    polys.append(Polygon([(1.20, 1.12),  (1.38, 1.14),  (1.47, 1.38),  (1.20, 1.38)]))
    polys.append(Polygon([(1.32, 1.28),  (1.32, 1.38),  (1.44, 1.38),  (1.44, 1.28)]))
    polys.append(Polygon([(1.69, 0.95),  (1.69, 1.33),  (2.21, 1.33),  (2.21, 0.95)]))

    # We want to take slices of the 3D volume to get a successive sequence of
    # 2D cross-sections, one at each degree. We do this by giving trimesh the
    # surface normal of the intersecting plane.
    sections = []
    degrees = np.arange(0.0, 360.0, angle_step)
    print("Generating cross sections...")
    for i in tqdm(range(0, len(degrees))):
        deg = degrees[i]

        # We add 90 since we want the surface normal at this angle (which is
        # 90 degrees away).
        x_norm = np.cos(np.radians(deg + 90))
        y_norm = np.sin(np.radians(deg + 90))

        # Can use this to take slices at each degree, returning the X, Y
        # coordinates.
        slice = trimesh.intersections.mesh_plane(mesh, plane_origin=[0,0,0],
          plane_normal=[x_norm, y_norm, 0])

        # Slice is a sequence of 3D lines, where each entry is the start and
        # and coordinates of the line [(X0, Y0, Z0), (X1, Y1, Z1)]. Each line
        # segment is connected, so each coordinate is repeated for the two
        # line segements it shares, but it's not consistent as to if the 1's
        # connect to the next segment's 0's or 1's, so we will go through and
        # pull out one coordinate at a time, avoiding duplicates.
        points = []
        for pointset in slice:
            for point in pointset:
                point = list(point)
                if point not in points:
                    points.append(point)
        points = np.array(points)

        # At this point can define R and Z, in m.
        r = np.sqrt(np.square(points[:,0]) + np.square(points[:,1])) / 1000
        z = points[:,2] / 1000

        # Keep the correct half of the cross-section.
        if deg > 0 and deg <= 180:
            keep = points[:,0] < 0
        else:
            keep = points[:,0] > 0
        x_cs = points[:,0][keep]
        y_cs = points[:,1][keep]
        r = r[keep]
        z = z[keep]

        # Sort using the middle of the vessel as a center. This doesn't get all
        # of them in the correct order, but it's a start.
        # Open thought: Could a sorting algorithm be made such that you divide
        # a region into subregions, and then sort points in a region clockwise,
        # where you change the number of regions as a knob? Hmmmm...
        r, z = sort_polar((1.1, 0.0), r, z, polys[0])

        # Re-sort the port locations. The middle port seems to get sorted just
        # fine via the first sort. This sorts the ports above and below the
        # midplane.
        r, z = sort_polar((2.00, 0.70), r, z, polys[1])
        r, z = sort_polar((2.00, -0.60), r, z, polys[2])

        # Re-sort the SAS region.
        r, z = sort_polar((1.44, 1.10), r, z, polys[3])

        # Re-sort the upper divertor area.
        r, z = sort_polar((1.30, 1.14), r, z, polys[4])

        # Re-sort the UOB area.
        r, z = sort_polar((1.32, 1.32), r, z, polys[5])

        # Re-sort this top port area.
        r, z = sort_polar((1.95, 1.00), r, z, polys[6])

        # Machine coordinates are clockwise so 360 - phi.
        phi_mach = 360 - deg

        sections.append({"degree":deg, "x_norm":x_norm, "y_norm":y_norm,
          "x_cs":x_cs, "y_cs":y_cs, "r":r, "z":z, "phi_mach":phi_mach})

    # Create the MAFOT input wall file. First see what the maximum number of
    # points are for a wall section.
    max_pts = 0
    for i in range(0, len(sections)):
        pts = len(sections[i]["r"])
        if pts > max_pts:
            max_pts = pts

    # Open file and print the header info.
    with open("mafot_3D_wall.dat", "w") as f:
        f.write("# DIII-D 3D wall: (R,Z) for every degree: phi = 0,...,359;" + \
          " leading integer gives number of points per plane\n")
        f.write("# The first integer gives the maximum number of points in" + \
          " any plane\n")
        f.write(str(max_pts) + "\n")

        # Then write the R, Z coordinates with the number of points first.
        for i in range(0, len(sections)):
            num_pts = len(sections[i]["r"])
            f.write(str(num_pts) + "\n")
            for j in range(0, num_pts):
                f.write("{:7.5f} {:7.5f}\n".format(sections[i]["r"][j], sections[i]["z"][j]))


    # Optional flag to make a gif of the wall as you progress around the vessel.
    if make_gif:
        import imageio
        import os

        fnames = []
        print("Creating gif...")
        for i in tqdm(range(0, len(sections))):
            fig, ax = plt.subplots()
            ax.plot(sections[i]["r"], sections[i]["z"], color="r", zorder=9)
            ax.scatter(sections[i]["r"], sections[i]["z"], color="k", zorder=10, s=5)
            for poly in polys[1:]:
                ax.plot(*poly.exterior.xy, color="b")

            ax.set_aspect("equal")
            ax.text(2.0, 1.3, r"$\phi$={:}$^\circ$".format(sections[i]["phi_mach"]))
            ax.set_xlim([0.95, 2.75])
            ax.set_ylim(-1.5, 1.5)

            fname = "wall_plots/{}.png".format(i)
            fnames.append(fname)
            fig.savefig(fname)
            plt.close(fig)

        # Now build the gif and clean up by removing the plots.
        print("Saving as gif...")
        with imageio.get_writer('wall.gif', mode="I") as writer:
            for fname in fnames:
                image = imageio.imread(fname)
                writer.append_data(image)
        for fname in fnames:
            os.remove(fname)


if __name__ == "__main__":

    print()
    print("*****************************************")
    print("* Generate MAFOT input file for 3D wall *")
    print("*****************************************")

    if len(sys.argv) in [2, 3]:
        stl_path = sys.argv[1]

        # Identify if we want to make the gif or not.
        if len(sys.argv) == 2:
            make_gif = False
        elif len(sys.argv) == 3:
            if sys.argv[2].lower() == "true":
                make_gif = True
            else:
                make_gif = False

        # Run main program.
        main(stl_path, make_gif=make_gif)

    else:
        print("Error! Only pass in the path to the STL file. If there are" + \
          "spaces in your path enclose the whole path in quotes.")
        print("Example: python3 generate_3d_wall_for_mafot path_to_stl.stl")
