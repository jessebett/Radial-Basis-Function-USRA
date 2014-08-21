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

def savefig(name):
    mlab.savefig('Figures/'+name)

def make_coor(n):
    '''Creates points on the surface of sphere using lat-lon grid points'''
    phi, theta = np.mgrid[0:pi:n, 0:2 * pi:n]
    Coor = namedtuple('Coor', 'r phi theta x y z')
    r = 1
    x = r * sin(phi) * cos(theta)
    y = r * sin(phi) * sin(theta)
    z = r * cos(phi)
    return Coor(r, phi, theta, x, y, z)

def makefinegrid():
    coor=make_coor(20j)

    mlab.figure()
    mlab.mesh(coor.x, coor.y, coor.z)
    mlab.points3d(coor.x, coor.y, coor.z, 
                      scale_factor=0.1)
    savefig('finegrid.png')


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

def make_sphere(coordinates):
    '''Create DiPy sphere object from sphere points in R3'''
    return Sphere(coordinates.x, coordinates.y, coordinates.z)



def makecoarsegrid():
    coor=make_uni_coor(100)
    sphere= make_sphere(coor)


    mlab.figure()
    mlab.points3d(coor.x, coor.y, coor.z, 
                      scale_factor=0.1)
    mlab.triangular_mesh(coor.x,coor.y,coor.z,sphere.faces)
    savefig('coarsegrid.png')
makecoarsegrid()