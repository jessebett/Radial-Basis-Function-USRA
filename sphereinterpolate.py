import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import special
from collections import namedtuple

# Create a sphere
pi = np.pi
cos = np.cos
sin = np.sin
fun= spheregrid(100j)
interpphi, interptheta = np.mgrid[0:pi:100j, 0:2 * pi:100j]

def spheregrid(n):
    Grid=namedtuple('Grid', 'phi theta')
    phi, theta = np.mgrid[0:pi:n, 0:2*pi:n]
    return Grid(phi,theta)

def coordinates(r,phi,theta):
    Coor=namedtuple('Coor', 'x y z')
    r = r
    x = r * sin(phi) * cos(theta)
    y = r * sin(phi) * sin(theta)
    z = r * cos(phi)
    return Coor(x,y,z)

colorfunction = special.sph_harm(3, 4, theta, phi).real
norm = colors.Normalize()

fun=coordinates(3, funphi, funtheta)
interp=coordinates(3, )

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(fun.x, fun.y, fun.z)
ax.plot_surface(
    fun.x, fun.y, fun.z,  rstride=1, cstride=1, facecolors=cm.jet(norm(colorfunction)))
fig.savefig('sphere.png')
plt.show()


