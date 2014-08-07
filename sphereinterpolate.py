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

# Creating a sphere in Cartesian and Sphereical
# Saves coordinates as named tuples


def coordinates(r, n):
    phi, theta = np.mgrid[0:pi:n, 0:2 * pi:n]
    Coor = namedtuple('Coor', 'r phi theta x y z')
    r = r
    x = r * sin(phi) * cos(theta)
    y = r * sin(phi) * sin(theta)
    z = r * cos(phi)
    return Coor(r, phi, theta, x, y, z)

# Creating a sphere
# fine is coordinates on a fine grid
# interp is coordinates on coarse grid for training interpolation
fine = coordinates(1, 100j)
interp = coordinates(1, 6j)


# Defining finection to colour sphere
# Here we are using a spherical harmonic
def harmonic(m, n, theta, phi):
    return special.sph_harm(m, n, theta, phi).real
norm = colors.Normalize()

# One example of the harmonic function, for testing
harmonic13_fine = harmonic(1, 3, fine.theta, fine.phi)
harmonic13_coarse = harmonic(1, 3, interp.theta, interp.phi)


# Train the interpolation using interp coordinates
rbf = Rbf(interp.phi, interp.theta, harmonic13_coarse)
# The result of the interpolation on fine coordinates
interp_values = rbf(fine.phi, fine.theta)

error = harmonic13_fine - interp_values
L_infinity = np.max(np.abs(error))

# rbf=Rbf(interp.x, interp.y, interp.z, harmonic13_coarse)
# interp_values=rbf(fine.x,fine.y,fine.z)


def make_figures(fine_coor, interp_coor, fine_function, coarse_function, interp_results):
    # Figure of harmoinc function on sphere in fine cordinates
    # Points3d showing interpolation training points coloured to their value
    mlab.figure()
    vmax, vmin = np.max(fine_function), np.min(fine_function)
    mlab.mesh(fine_coor.x, fine_coor.y, fine_coor.z,
              scalars=fine_function, vmax=vmax, vmin=vmin)
    mlab.points3d(interp_coor.x, interp_coor.y, interp_coor.z, coarse_function,
                  scale_factor=0.1, scale_mode='none', vmax=vmax, vmin=vmin)
    # mlab.savefig('interppointssphere.png')

    # Figure showing results of rbf interpolation
    mlab.figure()
    mlab.mesh(fine_coor.x, fine_coor.y, fine_coor.z,
              scalars=interp_results, vmax=vmax, vmin=vmin)
    mlab.points3d(interp_coor.x, interp_coor.y, interp_coor.z, coarse_function,
                  scale_factor=0.1, scale_mode='none', vmax=vmax, vmin=vmin)
    # mlab.savefig('interpsphere.png')

    mlab.figure()
    mlab.mesh(fine_coor.x, fine_coor.y, fine_coor.z,
              scalars=error, vmax=L_infinity, vmin=-L_infinity)
    mlab.colorbar(title='Error', orientation='vertical')
    # mlab.points3d(interp_coor.x, interp_coor.y, interp_coor.z, scalars, scale_factor=0.1, scale_mode='none',vmax=vmax, vmin=vmin)
    # mlab.savefig('interpsphere.png')

    mlab.show()

make_figures(fine, interp, harmonic13_fine, harmonic13_coarse, interp_values)


# A different norm potentially
def sphere_norm(x, y):
    return np.arccos(np.dot(x, y))
