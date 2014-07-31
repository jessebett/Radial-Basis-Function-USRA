import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import special
from collections import namedtuple

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

fun = coordinates(3, 100j)
interp = coordinates(3, 100j)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
# ax.scatter(fun.x, fun.y, fun.z)
ax.plot_surface(
    fun.x, fun.y, fun.z,  rstride=1, cstride=1, facecolors=cm.jet(norm(colorfunction(1, 3, fun.theta, fun.phi))))
fig.savefig('sphere.png')
plt.show()
