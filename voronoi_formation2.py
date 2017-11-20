import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import sys


n=0
points0=[]

#https://stackoverflow.com/questions/6583573/how-to-read-numbers-from-file-in-python
# This opens a handle to your file, in 'r' read mode
file_handle = open(sys.argv[1], 'r')
# Read in all the lines of your file into a list of lines
lines_list = file_handle.readlines()
# Do a double-nested list comprehension to get the rest of the data into your matrix
my_data = [[float(val) for val in line.split()] for line in lines_list[0:]]

for i in xrange(len(my_data)):
    points0.append([])
    points0[n].append(my_data[i][2]/68);
    points0[n].append(my_data[i][1]/105);
    n=n+1
    
    
#adding some points "at infinity" so that the voronoi diagram will cover the whole pitch

points0.append([])
points0[n].append(10)
points0[n].append(10)
n=n+1
points0.append([])
points0[n].append(-10)
points0[n].append(10)
n=n+1
points0.append([])
points0[n].append(-10)
points0[n].append(-10)
n=n+1
points0.append([])
points0[n].append(10)
points0[n].append(-10)


#this is stupid, left-over for historical reasons
points1 = np.array(points0)
Lx=68
Ly=105


for i in points1:
    i[0]=i[0]*Lx
    i[1]=i[1]*Ly


#https://gist.github.com/pv/8036995
def voronoi_finite_polygons_2d(vor, radius=None):
    """
    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = []
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):
        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions.append(vertices)
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions.append(new_region.tolist())

    return new_regions, np.asarray(new_vertices)

vor = Voronoi(points1)

regions, vertices = voronoi_finite_polygons_2d(vor)

dpi=100
Nx0=400
Nx=(Nx0/dpi)
Ny=(Nx0*Ly/Lx/dpi)

plt.figure(figsize=(Nx,Ny),dpi=100)

#https://gist.github.com/pv/8036995
# colorize
for region in regions:
    polygon = vertices[region]
    plt.fill(*zip(*polygon), alpha=0.4)

plt.plot(points1[:,0], points1[:,1], 'ko')
plt.tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='off', right='off', left='off', labelleft='off') 
linie = np.array([[0,Ly/2.],[Lx,Ly/2.]])
plt.plot(linie[:,0],linie[:,1],'k',marker='')

plt.xlim(-0.0*Lx, 1.0*Lx)
plt.ylim(-0.0*Lx, 1.0*Ly)



#plt.show()
plt.savefig(sys.argv[len(sys.argv)-1],dpi=dpi,bbox_inches='tight')
