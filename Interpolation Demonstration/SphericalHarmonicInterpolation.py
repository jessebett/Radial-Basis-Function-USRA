import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import pi, cos, sin
from scipy import special
from scipy.interpolate import Rbf
from collections import namedtuple
from mayavi import mlab
from math import asin, sqrt
from dipy.core.sphere import Sphere


# Creating a sphere in Cartesian and Sphereical
# Saves coordinates as named tuples
def coordinates(n_fine, n_coarse):
    def make_coor(n):
        phi, theta = np.mgrid[0:pi:n, 0:2 * pi:n]
        Coor = namedtuple('Coor', 'r phi theta x y z')
        r = 1
        x = r * sin(phi) * cos(theta)
        y = r * sin(phi) * sin(theta)
        z = r * cos(phi)
        return Coor(r, phi, theta, x, y, z)

    def uniform_spherical_distribution(N):
        """n points distributed evenly on the surface of a unit sphere"""
        pts = []
        r = 1
        inc = pi * (3 - sqrt(5))
        off = 2 / float(N)
        for k in range(0, int(N)):
            y = k * off - 1 + (off / 2)
            r = sqrt(1 - y * y)
            phi = k * inc
            pts.append([cos(phi) * r, y, sin(phi) * r])
        return np.array(pts)

    def appendSpherical_np(xyz):
        ptsnew = np.hstack((xyz, np.zeros(xyz.shape)))
        xy = xyz[:, 0] ** 2 + xyz[:, 1] ** 2
        ptsnew[:, 3] = np.sqrt(xy + xyz[:, 2] ** 2)
        # for elevation angle defined from Z-axis down
        ptsnew[:, 4] = np.arctan2(np.sqrt(xy), xyz[:, 2])
        # ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle
        # defined from XY-plane up
        ptsnew[:, 5] = np.arctan2(xyz[:, 1], xyz[:, 0])
        return ptsnew

    def make_uni_coor(n):
        Coor = namedtuple('Coor', 'theta phi x y z')
        pts = uniform_spherical_distribution(n)
        pts = appendSpherical_np(pts)

        return Coor(pts[:, 5], pts[:, 4], pts[:, 0], pts[:, 1], pts[:, 2])

    Coordinates = namedtuple('Coordinates', 'fine coarse ')
    # return Coordinates(make_coor(n_fine),make_coor(n_coarse))
    return Coordinates(make_coor(n_fine), make_uni_coor(n_coarse))

def make_sphere(coordinates):
    return Sphere(coordinates.x, coordinates.y, coordinates.z)
# Defining finection to colour sphere
# Here we are using a spherical harmonic
def harmonic(m, n, coor):
    Harmonic = namedtuple('Harmonic', 'fine coarse')
    return Harmonic(
        special.sph_harm(m, n, coor.fine.theta, coor.fine.phi).real,
        special.sph_harm(m, n, coor.coarse.theta, coor.coarse.phi).real
    )
norm = colors.Normalize()


def angle(x1, x2):
    xx = np.arccos((x1 * x2).sum(axis=0))
    xx[np.isnan(xx)] = 0
    return xx


def spherical_dist(pos1, pos2):
    r = 1
    pos1 = pos1[:, None]
    cos_lat1 = np.cos(pos1[0])
    cos_lat2 = np.cos(pos2[0])
    cos_lat_d = np.cos(pos1[0] - pos2[0])
    cos_lon_d = np.cos(pos1[1] - pos2[1])
    return r * np.arccos(cos_lat_d - cos_lat1 * cos_lat2 * (1 - cos_lon_d))[0]


def rbf_interpolate(coor, coarse_function, epsilon=None):
    # Train the interpolation using interp coordinates
    rbf = Rbf(coor.coarse.x, coor.coarse.y, 
              coor.coarse.z, coarse_function, norm=angle, epsilon=epsilon)
    print rbf.epsilon
    # The result of the interpolation on fine coordinates
    return rbf(coor.fine.x, coor.fine.y, coor.fine.z)


def interp_error(fine_function, interp_results):
    Error = namedtuple('Error', 'errors max')
    errors = fine_function - interp_results
    error_max = np.max(np.abs(errors))
    return Error(errors, error_max)


# rbf=Rbf(interp.x, interp.y, interp.z, harmonic13_coarse)
# interp_values=rbf(fine.x,fine.y,fine.z)


def make_figures(coor, fun, interp_results, error):
    # Figure of harmoinc function on sphere in fine cordinates
    # Points3d showing interpolation training points coloured to their value
    mlab.figure()
    vmax, vmin = np.max(fun.fine), np.min(fun.fine)
    mlab.mesh(coor.fine.x, coor.fine.y, coor.fine.z,
              scalars=fun.fine, vmax=vmax, vmin=vmin)
    mlab.points3d(coor.coarse.x, coor.coarse.y, coor.coarse.z, fun.coarse,
                  scale_factor=0.1, scale_mode='none', vmax=vmax, vmin=vmin)
    mlab.colorbar(title='Spherical Harmonic', orientation='vertical')
    # mlab.savefig('interppointssphere.png')

    # Figure showing results of rbf interpolation
    mlab.figure()
    mlab.mesh(coor.fine.x, coor.fine.y, coor.fine.z,
              scalars=interp_results, vmax=vmax, vmin=vmin)
    mlab.points3d(coor.coarse.x, coor.coarse.y, coor.coarse.z, fun.coarse,
                  scale_factor=0.1, scale_mode='none', vmax=vmax, vmin=vmin)
    mlab.colorbar(title='Interpolation', orientation='vertical')
    # mlab.savefig('interpsphere.png')

    mlab.figure()
    mlab.mesh(coor.fine.x, coor.fine.y, coor.fine.z,
              scalars=error.errors, vmax=error.max, vmin=-error.max)
    mlab.colorbar(title='Error', orientation='vertical')
    # mlab.points3d(coor.coarse.x, coor.coarse.y, coor.coarse.z, scalars, scale_factor=0.1, scale_mode='none',vmax=vmax, vmin=vmin)
    # mlab.savefig('interpsphere.png')

    mlab.show()


def harmonic_combo(*args):
    fine_combo = np.add([arg.fine for arg in args])
    coarse_combo = np.add([arg.coarse for arg in args])
    # for arg in args:
    #     print np.shape(arg.fine)
    #     print np.add(arg.fine)
    Combo = namedtuple('Combo', 'fine coarse')
    return Combo(fine_combo, coarse_combo)


# Creating a sphere
# fine is coordinates on a fine grid
# coarse is coordinates on coarse grid for training interpolation
coordinates = coordinates(300j, 350)

sphere = make_sphere(coordinates.coarse)

# One example of the harmonic function, for testing
function = harmonic(3, 4, coordinates)
# harmonic12= harmonic(2, 3, coordinates)

interp_values = rbf_interpolate(coordinates, function.coarse)
error = interp_error(function.fine, interp_values)


# make_figures(coordinates, function, interp_values, error)

errors=[]
epsilons = np.linspace(1.4,5,20)
for epsilon in epsilons:
    interp_values = rbf_interpolate(coordinates, function.coarse, epsilon=epsilon)
    error = interp_error(function.fine, interp_values)
    maxerror =  error.max
    errors.append(maxerror)

print errors
print np.shape(errors)
print np.shape(epsilons)
ax = plt.subplots()
plt.plot(epsilons,errors, 'bo')
optimum = errors.index(np.min(errors))
plt.plot(epsilons[optimum],errors[optimum],'ro')
plt.savefig('optimizationcurve.png')
plt.show()

print epsilons[optimum]


interp_values = rbf_interpolate(coordinates, function.coarse, epsilon=epsilons[optimum])
error = interp_error(function.fine, interp_values)
make_figures(coordinates, function, interp_values, error)