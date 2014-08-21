#Import libraries
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import pi, cos, sin
from math import asin, sqrt
#For interpolating
from scipy import special
from scipy.interpolate import Rbf
#For clean variable names
from collections import namedtuple
#For visualizing interpolation output
from mayavi import mlab
#For geometry on surface of sphere
from dipy.core.sphere import Sphere

#Function to handle generation of coordinate systems
#Coarse coordinates refer to data site coordinates for training RBF
#Fine coordinates refer to interpolation sites
def coordinates(n_fine, n_coarse):
    '''Creates fine and coarse coordiante systems'''
    def make_coor(n):
        '''Creates points on the surface of sphere using lat-lon grid points'''
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
        '''Appends spherical coordinates to array of Cartesian coordinates'''
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
        '''Make named tuple of unifromly distrubed points on sphere'''
        Coor = namedtuple('Coor', 'theta phi x y z')
        pts = uniform_spherical_distribution(n)
        pts = appendSpherical_np(pts)

        return Coor(pts[:, 5], pts[:, 4], pts[:, 0], pts[:, 1], pts[:, 2])

    Coordinates = namedtuple('Coordinates', 'fine coarse ')

    return Coordinates(make_coor(n_fine), make_uni_coor(n_coarse))

def make_sphere(coordinates):
    '''Create DiPy sphere object from sphere points in R3'''
    return Sphere(coordinates.x, coordinates.y, coordinates.z)

#Defining the spherical harmonic to interpolate
def harmonic(m, l, coor):
    '''Produce m,l spherical harmonic at coarse and fine coordinates'''
    Harmonic = namedtuple('Harmonic', 'fine coarse')
    return Harmonic(
        special.sph_harm(m, l, coor.fine.theta, coor.fine.phi).real,
        special.sph_harm(m, l, coor.coarse.theta, coor.coarse.phi).real
    )

def angle(x1, x2):
    '''Distance metric on the surface of the unit sphere'''
    xx = np.arccos((x1 * x2).sum(axis=0))
    xx[np.isnan(xx)] = 0
    return xx

def rbf_interpolate(coor, coarse_function, epsilon=None):
    '''Radial Basis Function Interpolation from coarse sites to fine cooridnates'''
    # Train the interpolation using interp coordinates
    rbf = Rbf(coor.coarse.x, coor.coarse.y, 
              coor.coarse.z, coarse_function, norm=angle, epsilon=epsilon)
    # The result of the interpolation on fine coordinates
    return rbf(coor.fine.x, coor.fine.y, coor.fine.z)


def interp_error(fine_function, interp_results):
    '''Error between interpolated function and actual function'''
    Error = namedtuple('Error', 'errors max')
    errors = fine_function - interp_results
    error_max = np.max(np.abs(errors))
    return Error(errors, error_max)




def make_figures(coor, fun, interp_results, error):
    '''Produce MayaVi figures for interpolation results'''
    mlab.figure()
    vmax, vmin = np.max(fun.fine), np.min(fun.fine)
    mlab.mesh(coor.fine.x, coor.fine.y, coor.fine.z,
              scalars=fun.fine, vmax=vmax, vmin=vmin)
    mlab.points3d(coor.coarse.x, coor.coarse.y, coor.coarse.z, fun.coarse,
                  scale_factor=0.1, scale_mode='none', vmax=vmax, vmin=vmin)
    mlab.colorbar(title='Spherical Harmonic', orientation='vertical')
    mlab.savefig('Figures/functionsphere.png')

    # Figure showing results of rbf interpolation
    mlab.figure()
    mlab.mesh(coor.fine.x, coor.fine.y, coor.fine.z,
              scalars=interp_results, vmax=vmax, vmin=vmin)
    mlab.points3d(coor.coarse.x, coor.coarse.y, coor.coarse.z, fun.coarse,
                  scale_factor=0.1, scale_mode='none', vmax=vmax, vmin=vmin)
    mlab.colorbar(title='Interpolation', orientation='vertical')
    mlab.savefig('Figures/interpsphere.png')

    mlab.figure()
    mlab.mesh(coor.fine.x, coor.fine.y, coor.fine.z,
              scalars=error.errors, vmax=error.max, vmin=-error.max)
    mlab.colorbar(title='Error', orientation='vertical')
    # mlab.points3d(coor.coarse.x, coor.coarse.y, coor.coarse.z, scalars, scale_factor=0.1, scale_mode='none',vmax=vmax, vmin=vmin)
    mlab.savefig('Figures/errorsphere.png')

    mlab.show()



# Create coordinate system
#Fine coordinate given as complex number
coordinates = coordinates(300j, 350)

#Create sphere object
sphere = make_sphere(coordinates.coarse)
                                                                                                                                                                    
# Harmonic function we'll interpolate
function = harmonic(3, 4, coordinates)

#Preform the RBF interpolation
interp_values = rbf_interpolate(coordinates, function.coarse)

#Test the error between function and interpolation
error = interp_error(function.fine, interp_values)

#Create the MayaVi Figures
# make_figures(coordinates, function, interp_values, error)


#Optimize choice of epsilon to minimize error
errors=[]
epsilons = np.linspace(1.4,5,20)
for epsilon in epsilons:
    interp_values = rbf_interpolate(coordinates, function.coarse, epsilon=epsilon)
    error = interp_error(function.fine, interp_values)
    maxerror =  error.max
    errors.append(maxerror)

fig, ax = plt.subplots()
fig.suptitle('RBF Interpolation Epsilon Optimization')
ax.set_xlabel('Shape Parameter, $\epsilon$')
ax.set_ylabel('Maximum Error')
plt.plot(epsilons,errors, 'bo')
optimum = errors.index(np.min(errors))
plt.plot(epsilons[optimum],errors[optimum],'ro')
plt.savefig('Figures/optimizationcurve')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
plt.show()



interp_values = rbf_interpolate(coordinates, function.coarse, epsilon=epsilons[optimum])
error = interp_error(function.fine, interp_values)
make_figures(coordinates, function, interp_values, error)