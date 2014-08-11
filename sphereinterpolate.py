import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import special
from scipy.interpolate import Rbf
from collections import namedtuple
from mayavi import mlab
from math import asin

# Nice aliases
pi = np.pi
cos = np.cos
sin = np.sin

# Creating a sphere in Cartesian and Sphereical
# Saves coordinates as named tuples


def coordinates(n_fine, n_coarse):
    def make_coor(n):
        phi, theta = np.mgrid[0.1:pi-0.1:n, 0.1:2 * pi-0.1:n]
        phi = np.ravel(phi)
        theta = np.ravel(theta)
        Coor = namedtuple('Coor', 'r phi theta x y z')
        r = 1
        x = r * sin(phi) * cos(theta)
        y = r * sin(phi) * sin(theta)
        z = r * cos(phi)
        return Coor(r, phi, theta, x, y, z)
    def rand_sphere(n):
         """n points distributed evenly on the surface of a unit sphere""" 
        z = 2 * random.rand(n) - 1   # uniform in -1, 1
        t = 2 * pi * random.rand(n)   # uniform in 0, 2*pi
         x = sqrt(1 - z**2) * cos(t)
         y = sqrt(1 - z**2) * sin(t)
         return x, y, z
    def appendSpherical_np(xyz):
        ptsnew = np.hstack((xyz, np.zeros(xyz.shape)))
        xy = xyz[:,0]**2 + xyz[:,1]**2
        ptsnew[:,3] = np.sqrt(xy + xyz[:,2]**2)
        ptsnew[:,4] = np.arctan2(np.sqrt(xy), xyz[:,2]) # for elevation angle defined from Z-axis down
        #ptsnew[:,4] = np.arctan2(xyz[:,2], np.sqrt(xy)) # for elevation angle defined from XY-plane up
        ptsnew[:,5] = np.arctan2(xyz[:,1], xyz[:,0])
        return ptsnew
    

    Coordinates= namedtuple('Coordinates', 'fine coarse')
    # return Coordinates(make_coor(n_fine),make_coor(n_coarse))
    return Coordinates(rand(n_fine),make_coor(n_coarse))



# Defining finection to colour sphere
# Here we are using a spherical harmonic
def harmonic(m, n, coor):
    Harmonic = namedtuple('Harmonic', 'fine coarse')
    return Harmonic(
        special.sph_harm(m, n, coor.fine.theta, coor.fine.phi).real, 
        special.sph_harm(m, n, coor.coarse.theta, coor.coarse.phi).real
        )
norm = colors.Normalize()




def spherical_dist(pos1, pos2):
    r=1
    pos1=pos1[:,None]
    cos_lat1 = np.cos(pos1[0])
    cos_lat2 = np.cos(pos2[0])
    cos_lat_d = np.cos(pos1[0] - pos2[0])
    cos_lon_d = np.cos(pos1[1] - pos2[1])
    return r * np.arccos(cos_lat_d - cos_lat1 * cos_lat2 * (1 - cos_lon_d))[0]


def rbf_interpolate(coor, coarse_function):
    # Train the interpolation using interp coordinates
    rbf = Rbf(coor.coarse.phi, coor.coarse.theta, coarse_function,norm=spherical_dist)
    # The result of the interpolation on fine coordinates
    return rbf(coor.fine.phi, coor.fine.theta)


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
    Combo=namedtuple('Combo', 'fine coarse')
    return Combo(fine_combo,coarse_combo)




# Creating a sphere
# fine is coordinates on a fine grid
# coarse is coordinates on coarse grid for training interpolation
coordinates = coordinates(100j, 6j)

# One example of the harmonic function, for testing
function= harmonic(3, 4, coordinates)
# harmonic12= harmonic(2, 3, coordinates)


interp_values = rbf_interpolate(coordinates, function.coarse)
error = interp_error(function.fine, interp_values)

make_figures(coordinates, function, interp_values, error)





# A different norm potentially
def periodic_norm(x1, x2):
    def norm(x1,x2,bound):
       return bound/2-np.abs(bound/2-np.abs(x1-x2))

    t=norm(x1[0],x2[0],pi)
    p=norm(x1[1],x2[1],2*pi)

    return np.linalg.norm(zip(t,p))