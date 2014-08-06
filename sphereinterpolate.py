import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import special
from scipy.interpolate import Rbf
from collections import namedtuple
from mayavi import mlab

# Nice aliases
pi = np.pi
cos = np.cos
sin = np.sin

# Sphere Making Function


def coordinates(r, n):
    phi, theta = np.mgrid[0:pi:n, 0:2 * pi:n]
    Coor = namedtuple('Coor', 'r phi theta x y z')
    r = r
    x = r * sin(phi) * cos(theta)
    y = r * sin(phi) * sin(theta)
    z = r * cos(phi)
    return Coor(r, phi, theta, x, y, z)

# Creating the sphere
fun = coordinates(3, 100j)
interpphi, interptheta = np.mgrid[0:pi:100j, 0:2 * pi:100j]


# Defining function to colour sphere
def colorfunction(m, n, theta, phi):
    return special.sph_harm(m, n, theta, phi).real
norm = colors.Normalize()

fun = coordinates(1, 100j)
interp = coordinates(1, 5j)

def sphere_norm(x,y):
    return np.arccos(np.dot(x,y))

harmonic=colorfunction(1, 3, fun.theta, fun.phi)

rbf=Rbf(interp.phi, interp.theta, colorfunction(1, 3, interp.theta, interp.phi))
interp_values=rbf(fun.phi,fun.theta)


# rbf=Rbf(interp.x, interp.y, interp.z, colorfunction(1, 3, interp.theta, interp.phi))
# interp_values=rbf(fun.x,fun.y,fun.z)

mlab.figure()
mlab.mesh(fun.x, fun.y, fun.z,scalars=harmonic)
mlab.points3d(interp.x, interp.y, interp.z,scale_factor=0.1)
mlab.show()



# fig = plt.figure(figsize=(13,5))

# #Function Figure
# ax1 = fig.add_subplot(121, projection='3d')
# ax1.scatter(
#     interp.x, interp.y, interp.z, 
#     color='r', marker='o',
#     zorder=1, 
#     s=40
#     )

# ax1.plot_surface(
#     fun.x, fun.y, fun.z,  
#     rstride=1, 
#     cstride=1, 
#     alpha=0.5,
#     facecolors=cm.jet(norm(harmonic)),
#     zorder=-1
#     )

# ax2 = fig.add_subplot(122, projection='3d')
# ax2.plot_surface(
#     fun.x, fun.y, fun.z,  
#     rstride=1, 
#     cstride=1, 
#     alpha=1,
#     facecolors=cm.jet(norm(interp_values))
#     )
# fig.savefig('sphere.png')
# plt.show()
