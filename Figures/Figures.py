
# The SciPy interpolate library has a Radial Basis Function called Rbf
from scipy.interpolate import Rbf
import numpy as np
import matplotlib.pyplot as plt

plt.close("all")


def rgb(colour):
    return tuple(x / 255.0 for x in colour)
buf = 0.4
grey = '#ECE4D8' #
pink = '#cd596b' #
orange = '#a61717' #
blue = '#557398' #
red = '#D25252'
green = '#009926' #
white = '#8b786f' #


def jesseaxis(ax, x, y):
    ax.set_axis_bgcolor(grey)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['bottom'].set_color(white)
    ax.spines['left'].set_color(white)
    ax.set_xlim(min(x) - buf, max(x) + buf)
    ax.set_ylim(min(y) - buf, max(y) + buf)


def saveimg(name):
    plt.savefig(name + '.png', bbox_inches='tight',
                pad_inches=0.1, facecolor=grey, dpi=500)

####FIGURE SHOWING DEFINITION OF INTERPOLATION####
#FIG 1#
def makeinterpdef():
    # Defining the data sites
    fig1 = plt.figure(facecolor=grey)
    ax1 = plt.axes()

    x = np.linspace(-6, 6, 20)
    y = (5 * np.sin(x) + x)
    y = y + max(y)
    ax1.scatter(x, y, color=pink, marker='o')
    ax1.scatter(x, [-buf] * np.size(y), color=orange, marker='o')
    ax1.vlines(x, -buf, y, pink, linestyles='dashed', alpha=0.4)
    ax1.annotate('$x_i$', (-0.3, 0.4), size=30, color=orange)
    ax1.annotate('$f_i$', (-0.30, 9), size=30, color=pink)
    ax1.annotate('$s(x)$', (4, 9), size=30, color=blue)
    jesseaxis(ax1, x, y)
    # Defining our interpolation function from the data sites
    rbf = Rbf(x, y, epsilon=0.2)
    # Applying the interpolation
    xi = np.linspace(-6, 6, 100)
    yi = rbf(xi)
    plt.plot(xi, yi, color=blue)
    saveimg('interpdef')
    # plt.show()


###INTERP VS APPROX###
#FIG 2#
def makeinterpvsapprox():
    plt.close()
    x = np.array([0.0, 1.5, 2.0, 3.0,  4.0,  4.5])
    y = np.array([0.0, 0.8, 0.9, 0.1, -0.8, -1.0])
    z = np.polyfit(x, y, 3)
    p = np.poly1d(z)
    p30 = np.poly1d(np.polyfit(x, y, 5))
    xp = np.linspace(-2, 6, 100)
    ax2 = plt.axes()
    plt.plot(x, y, 'o', color=pink)
    plt.plot(xp, p(xp), '--', color=green)
    plt.plot(xp, p30(xp), '-', color=orange)
    ax2.annotate('Interpolation', (3.3, 0.5), size=20, color=orange)
    ax2.annotate('Approximation', (0, -0.5), size=20, color=green)
    plt.ylim(-2, 2)

    jesseaxis(ax2, x, y)
    saveimg('interpvsapprox')
    # plt.show()

###POLYNOMIAL INTERP###
#FIG 3#
def makepolyinterp():
    plt.close()
    x = np.array([0.0, 1.5, 2.0, 3.0,  4.0,  4.5])
    y = np.array([0.0, 0.8, 0.9, 0.1, -0.8, -1.0])
    z = np.polyfit(x, y, 3)
    p = np.poly1d(z)
    p30 = np.poly1d(np.polyfit(x, y, 5))
    xp = np.linspace(-2, 6, 100)
    ax3 = plt.axes()
    plt.plot(x, y, 'o', color=pink)
    #plt.plot(xp, p(xp), '--', color=green)
    plt.plot(xp, p30(xp), '-', color=orange)
    # ax3.annotate('$s(x)=-0.02988 x^5 + 0.417 x^4 - 2.018 x^3 + 3.694 x^2 - 1.722 x - 5.511e^{-14}$', (-0.2, -1.2), size=13, color=orange)
    ax3.annotate('$s(x)$', (1.5, 0.3), size=30, color=orange)
    plt.ylim(-2, 2)

    jesseaxis(ax3, x, y)
    saveimg('polyinterp')
    # plt.show()

###BASIC FUNCTION###
#FIG 4#
def makebasicfun():
    plt.close()
    x = np.linspace(-4, 4, 1000)
    ax = plt.axes()
    plt.plot(x, np.abs(x), '-', color=green)
    plt.plot(0, 0, 'o', color=orange)

    # ax3.annotate('$s(x)=-0.02988 x^5 + 0.417 x^4 - 2.018 x^3 + 3.694 x^2 - 1.722 x - 5.511e^{-14}$', (-0.2, -1.2), size=13, color=orange)
    ax.annotate('Basic Function', (-1.9, 2.59), size=30, color=green)
    ax.annotate('Center', (-0.9, -0.5), size=30, color=orange)

    jesseaxis(ax, x, np.abs(x))
    plt.ylim(-1, 3)
    saveimg('basicfun')
    #plt.show()

###SHOWING BASIC FUNCTION AS BASIS###
#FIG 5#
def makebasicbasis():
    plt.close()
    # Defining the data sites
    fig1 = plt.figure(facecolor=grey)
    ax1 = plt.axes()

    x = np.linspace(-6, 6, 20)
    y = (5 * np.sin(x) + x)
    y = y + max(y)
    ax1.scatter(x, y, color=pink, marker='o')
    ax1.scatter(x, [-buf] * np.size(y), color=orange, marker='o')
    ax1.vlines(x, -buf, y, pink, linestyles='dashed', alpha=0.4)
    ax1.annotate('$x_i$', (-0.3, 0.4), size=30, color=orange)
    ax1.annotate('$f_i$', (-0.30, 7.5), size=30, color=pink)
    ax1.annotate('$\psi_i$', (4, 1), size=30, color=green)

    plt.plot(x, np.abs(x - x[10]) - buf, '-', color=green)
    # ax1.annotate('$s(x)$', (4, 9), size=30, color=blue)
    jesseaxis(ax1, x, y)
    # Defining our interpolation function from the data sites
    # rbf = Rbf(x, y, epsilon=0.2)
    # Applying the interpolation
    # xi = np.linspace(-6, 6, 100)
    # yi = rbf(xi)
    # plt.plot(xi, yi, color=blue)
    saveimg('basicbasis')
    #plt.show()

###BASIC FUNCTION WITH Xi AT CENTER###
#FIG 6#
def makebasicfunxi():
    plt.close()
    x = np.linspace(-4, 4, 1000)
    ax = plt.axes()
    plt.plot(x, np.abs(x), '-', color=green)
    plt.plot(0, 0, 'o', color=orange)

    ax.annotate('Basic Function', (-1.9, 2.59), size=30, color=green)
    ax.annotate('$x_i$', (-0.2, -0.4), size=30, color=orange)

    jesseaxis(ax, x, np.abs(x))
    plt.ylim(-1, 3)
    saveimg('basicfunxi')
    #plt.show()

###MULTIQUADRIC KERNEL AS BASIC FUNCTION###
#FIG 7#
def makekernelfun():
    plt.close()
    x = np.linspace(-4, 4, 1000)
    ax = plt.axes()
    plt.plot(x, np.sqrt(x ** 2 + 0.4 ** 2), '-', color=green)
    plt.plot(0, 0, 'o', color=orange)
    plt.vlines(0, 0.4, 0, pink, linestyles='dashed', alpha=0.7)

    ax.annotate('Multiquadric Kernel', (-2.45, 2.59), size=30, color=green)
    ax.annotate('Center', (-0.8, -0.4), size=30, color=orange)
    ax.annotate('$c$', (0.2, 0.1), size=30, color=pink)
    jesseaxis(ax, x, x)
    plt.ylim(-1, 3)
    saveimg('kernelfun')
    #plt.show()

###MULTIQUADRIC KERNEL AS BASIS FUNCTION##
#FIG 8#
def makekernelbasis():
    # Defining the data sites
    plt.close()
    fig1 = plt.figure(facecolor=grey)
    ax1 = plt.axes()

    x = np.linspace(-6, 6, 20)
    xi = np.linspace(-6, 6, 200)
    y = (5 * np.sin(x) + x)
    y = y + max(y)
    ax1.scatter(x, y, color=pink, marker='o')
    ax1.scatter(x, [-buf] * np.size(y), color=orange, marker='o')
    ax1.vlines(x, -buf, y, pink, linestyles='dashed', alpha=0.4)
    ax1.annotate('$x_i$', (-0.3, 0.4), size=30, color=orange)
    ax1.annotate('$f_i$', (-0.30, 7.5), size=30, color=pink)
    ax1.annotate('$\psi_i$', (4, 1), size=30, color=green)

    plt.plot(xi, np.sqrt((xi - x[10]) ** 2 + 0.4 ** 2) - buf, '-', color=green)
    # ax1.annotate('$s(x)$', (4, 9), size=30, color=blue)
    jesseaxis(ax1, x, y)
    # Defining our interpolation function from the data sites
    # rbf = Rbf(x, y, epsilon=0.2)
    # Applying the interpolation
    # xi = np.linspace(-6, 6, 100)
    # yi = rbf(xi)
    # plt.plot(xi, yi, color=blue)
    saveimg('kernelbasis')
    #plt.show()

#MultiQuadric Kernel#
def makeMultiQuadric():
    plt.close()
    x = np.linspace(-4, 4, 1000)
    ax = plt.axes()
    plt.plot(x, np.sqrt(1 + (np.abs(x) * (1 / 0.8)) ** 2), '-', color=green, linewidth=3.5)
    plt.plot(0, 0, 'o', color=orange)
    jesseaxis(ax, x, x)
    ax.axis('off')
    plt.ylim(-0.2, 3)
    saveimg('multiquadric')
    # plt.show()

#Gaussian Kernel#
def makeGaussian():
    plt.close()
    x = np.linspace(-4, 4, 1000)
    ax = plt.axes()
    plt.plot(x, np.exp(-(np.abs(x) * 0.4) ** 2), '-', color=green, linewidth=3.5)
    plt.plot(0, 0, 'o', color=orange)
    jesseaxis(ax, x, x)
    ax.axis('off')
    plt.ylim(-0.2, 1.2)
    saveimg('gaussian')
    # plt.show()

#Inverse Quadratic Kernel#
def makeInverseQuadratic():
    plt.close()
    x = np.linspace(-4, 4, 1000)
    ax = plt.axes()
    plt.plot(x, (1 + (np.abs(x) * 0.4) ** 2) ** (-1), '-', color=green, linewidth=3.5)
    plt.plot(0, 0, 'o', color=orange)
    jesseaxis(ax, x, x)
    ax.axis('off')
    plt.ylim(-0.2, 1.2)
    saveimg('inversequadratic')
    # plt.show()

#inverse MultiQuadric Kernel#
def makeInverseMultiQuadric():
    plt.close()
    x = np.linspace(-4, 4, 1000)
    ax = plt.axes()
    plt.plot(x, np.sqrt(1 + (np.abs(x) * 0.4) ** 2) ** (-1), '-', color=green, linewidth=3.5)
    plt.plot(0, 0, 'o', color=orange)
    jesseaxis(ax, x, x)
    ax.axis('off')
    plt.ylim(-0.2, 1.2)
    saveimg('inversemultiquadric')
    # plt.show()

#Make the above kernels#
def makekernels():
    makeGaussian()
    makeMultiQuadric()
    makeInverseQuadratic()
    makeInverseMultiQuadric()

def quickpiece(x):
    return np.piecewise(x,[np.abs(x)<1,np.abs(x)>=1],[1,0])


#Demonstrate conditioning problems in RBF#
def makeconditioned():
    plt.close()
    eps1=1
    eps2=2
    eps3=3

    # Defining the data sites
    fig1 = plt.figure(1,facecolor=grey)
    ax1 = plt.axes()
    x = np.linspace(-2,2,20)
    y = quickpiece(x)
    y = y + max(y)
    ax1.scatter(x, y, color=pink, marker='o')
    ax1.scatter(x, [-buf] * np.size(y), color=orange, marker='o')
    ax1.vlines(x, -buf, y, pink, linestyles='dashed', alpha=0.4)
    rbf = Rbf(x, y, epsilon=eps1)
    # Applying the interpolation
    xi = np.linspace(-2,2, 1000)
    yi = rbf(xi)
    ax1.annotate('$\epsilon=$'+str(eps1), (0,max(yi)), size=30, color=pink)
    plt.plot(xi, yi, color=blue)
    jesseaxis(ax1, x, yi)
    saveimg('conditioned')
    # plt.show()
    # plt.close(1)

    fig2 = plt.figure(2,facecolor=grey)
    ax1 = plt.axes()
    x = np.linspace(-2,2,20)
    y = quickpiece(x)
    y = y + max(y)
    ax1.scatter(x, y, color=pink, marker='o')
    ax1.scatter(x, [-buf] * np.size(y), color=orange, marker='o')
    ax1.vlines(x, -buf, y, pink, linestyles='dashed', alpha=0.4)
    rbf = Rbf(x, y, epsilon=eps2)
    # Applying the interpolation
    xi = np.linspace(-2,2, 1000)
    yi = rbf(xi)
    ax1.annotate('$\epsilon=$'+str(eps2), (0,max(yi)), size=30, color=pink)
    plt.plot(xi, yi, color=blue)
    jesseaxis(ax1, x, yi)
    saveimg('illconditioned')

    fig3 = plt.figure(3,facecolor=grey)
    ax1 = plt.axes()
    x = np.linspace(-2,2,20)
    y = quickpiece(x)
    y = y + max(y)
    ax1.scatter(x, y, color=pink, marker='o')
    ax1.scatter(x, [-buf] * np.size(y), color=orange, marker='o')
    ax1.vlines(x, -buf, y, pink, linestyles='dashed', alpha=0.4)
    rbf = Rbf(x, y, epsilon=eps3)
    # Applying the interpolation
    xi = np.linspace(-2,2, 1000)
    yi = rbf(xi)
    ax1.annotate('$\epsilon=$'+str(eps3), (0,max(yi)), size=30, color=pink)
    plt.plot(xi, yi, color=blue)
    jesseaxis(ax1, x, yi)
    saveimg('veryillconditioned')


    #plt.show()

#Demonstrate how epsillon changes the basis function#
def makebasisgaus():
    # Defining the data sites
    fig1 = plt.figure(facecolor=grey)
    ax1 = plt.axes()

    x = np.linspace(-6, 6, 20)
    xi = np.linspace(-6, 6, 200)
    y = (5 * np.sin(x) + x)
    y = y + max(y)
    ax1.scatter(x, y, color=pink, marker='o')
    ax1.scatter(x, [-buf] * np.size(y), color=orange, marker='o')
    ax1.vlines(x, -buf, y, pink, linestyles='dashed', alpha=0.4)
    ax1.annotate('$x_i$', (-0.3, 0.4), size=30, color=orange)
    ax1.annotate('$f_i$', (-0.30, 7.5), size=30, color=pink)
    ax1.annotate('$\psi_i$', (4, 1), size=30, color=green)

    plt.plot(xi, np.exp(-(0.4*xi)**2) - buf, '-', color=green)
    # ax1.annotate('$s(x)$', (4, 9), size=30, color=blue)
    jesseaxis(ax1, x, y)
    # Defining our interpolation function from the data sites
    # rbf = Rbf(x, y, epsilon=0.2)
    # Applying the interpolation
    # xi = np.linspace(-6, 6, 100)
    # yi = rbf(xi)
    # plt.plot(xi, yi, color=blue)
    saveimg('basisgaus1')

    fig2 = plt.figure(facecolor=grey)
    ax2 = plt.axes()

    x = np.linspace(-6, 6, 20)
    xi = np.linspace(-6, 6, 200)
    y = (5 * np.sin(x) + x)
    y = y + max(y)
    ax2.scatter(x, y, color=pink, marker='o')
    ax2.scatter(x, [-buf] * np.size(y), color=orange, marker='o')
    ax2.vlines(x, -buf, y, pink, linestyles='dashed', alpha=0.4)
    ax2.annotate('$x_i$', (-0.3, 0.4), size=30, color=orange)
    ax2.annotate('$f_i$', (-0.30, 7.5), size=30, color=pink)
    ax2.annotate('$\psi_i$', (4, 1), size=30, color=green)

    plt.plot(xi, np.exp(-(1*xi)**2) - buf, '-', color=green)
    # ax2.annotate('$s(x)$', (4, 9), size=30, color=blue)
    jesseaxis(ax2, x, y)
    # Defining our interpolation function from the data sites
    # rbf = Rbf(x, y, epsilon=0.2)
    # Applying the interpolation
    # xi = np.linspace(-6, 6, 100)
    # yi = rbf(xi)
    # plt.plot(xi, yi, color=blue)
    saveimg('basisgaus2')


    #plt.show()

def makeallfigs():
    makeinterpdef()
    makeinterpvsapprox()
    makepolyinterp()
    makebasicfun()
    makebasicbasis()
    makebasicfunxi()
    makekernelfun()
    makekernelbasis()
    makekernels()
    makeconditioned()
    makebasisgaus()

makeallfigs()
# makeGaussian()