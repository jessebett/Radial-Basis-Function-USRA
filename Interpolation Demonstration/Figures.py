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

    mlab.figure(bgcolor=None)
    mlab.mesh(coor.x, coor.y, coor.z)
    mlab.points3d(coor.x, coor.y, coor.z, 
                      scale_factor=0.1)
    mlab.savefig('finegrid.png')
makefinegrid()