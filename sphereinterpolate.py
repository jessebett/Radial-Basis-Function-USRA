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
    Harmonic = namedtuple('Harmonic', 'fine coarse')
    return Harmonic(special.sph_harm(m, n, fine.theta, fine.phi).real, special.sph_harm(m, n, coarse.theta, coarse.phi).real)
norm = colors.Normalize()


def rbf_interpolate(fine_coor, interp_coor, coarse_function):
    # Train the interpolation using interp coordinates
    rbf = Rbf(interp_coor.phi, interp_coor.theta, coarse_function)
    # The result of the interpolation on fine coordinates
    return rbf(fine_coor.phi, fine_coor.theta)


def interp_error(fine_function, interp_results):
    Error = namedtuple('Error', 'errors max')
    errors = fine_function - interp_results
    error_max = np.max(np.abs(errors))
    return Error(errors, error_max)


# rbf=Rbf(interp.x, interp.y, interp.z, harmonic13_coarse)
# interp_values=rbf(fine.x,fine.y,fine.z)


def make_figures(fine_coor, interp_coor, fine_function, coarse_function, interp_results, error):
    # Figure of harmoinc function on sphere in fine cordinates
    # Points3d showing interpolation training points coloured to their value
    mlab.figure()
    vmax, vmin = np.max(fine_function), np.min(fine_function)
    mlab.mesh(fine_coor.x, fine_coor.y, fine_coor.z,
              scalars=fine_function, vmax=vmax, vmin=vmin)
    mlab.points3d(interp_coor.x, interp_coor.y, interp_coor.z, coarse_function,
                  scale_factor=0.1, scale_mode='none', vmax=vmax, vmin=vmin)
    mlab.colorbar(title='Spherical Harmonic', orientation='vertical')
    # mlab.savefig('interppointssphere.png')

    # Figure showing results of rbf interpolation
    mlab.figure()
    mlab.mesh(fine_coor.x, fine_coor.y, fine_coor.z,
              scalars=interp_results, vmax=vmax, vmin=vmin)
    mlab.points3d(interp_coor.x, interp_coor.y, interp_coor.z, coarse_function,
                  scale_factor=0.1, scale_mode='none', vmax=vmax, vmin=vmin)
    mlab.colorbar(title='Interpolation', orientation='vertical')
    # mlab.savefig('interpsphere.png')

    mlab.figure()
    mlab.mesh(fine_coor.x, fine_coor.y, fine_coor.z,
              scalars=error.errors, vmax=error.max, vmin=-error.max)
    mlab.colorbar(title='Error', orientation='vertical')
    # mlab.points3d(interp_coor.x, interp_coor.y, interp_coor.z, scalars, scale_factor=0.1, scale_mode='none',vmax=vmax, vmin=vmin)
    # mlab.savefig('interpsphere.png')

    mlab.show()


# One example of the harmonic function, for testing
harmonic13_fine = harmonic(1, 3, fine.theta, fine.phi)
harmonic13_coarse = harmonic(1, 3, interp.theta, interp.phi)
interp_values = rbf_interpolate(fine, interp, harmonic13_coarse)
error = interp_error(harmonic13_fine, interp_values)

make_figures(
    fine, interp, harmonic13_fine, harmonic13_coarse, interp_values, error)


# A different norm potentially
def sphere_norm(x, y):
    return np.arccos(np.dot(x, y))
