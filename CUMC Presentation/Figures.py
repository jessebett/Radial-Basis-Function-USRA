
# The SciPy interpolate library has a Radial Basis Function called Rbf
from scipy.interpolate import Rbf
import numpy as np
import matplotlib.pyplot as plt

plt.close("all")


def rgb(colour):
    return tuple(x / 255.0 for x in colour)
buf = 0.4
grey = '#181818'
pink = '#FC0964'
orange = '#FD971F'
blue = '#66D9EF'
red = '#D25252'
green = '#7FB347'


def jesseaxis(ax,x,y):
    ax.set_axis_bgcolor(grey)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['bottom'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.set_xlim(min(x) - buf, max(x) + buf)
    ax.set_ylim(min(y) - buf, max(y) + buf)

def saveimg(name):
    plt.savefig('Images/'+name+'.png', bbox_inches='tight',
               pad_inches=0.1, facecolor=grey, dpi=500)

####FIGURE 1####
def makefig1():
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
    jesseaxis(ax1,x,y)
    # Defining our interpolation function from the data sites
    rbf = Rbf(x, y, epsilon=0.2)
    # Applying the interpolation
    xi = np.linspace(-6, 6, 100)
    yi = rbf(xi)
    plt.plot(xi, yi, color=blue)
    saveimg('fig1')
    # plt.show()


###INTERP VS APPROX###
def makefig2():
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

    jesseaxis(ax2,x,y)
    saveimg('fig2')
    # plt.show()

###POLYNOMIAL INTERP###


def makefig3():
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

    jesseaxis(ax3,x,y)
    saveimg('fig3')
    # plt.show()

###POLYNOMIAL INTERP###


def makefig4():
    plt.close()
    x = np.linspace(-4, 4, 1000)
    ax = plt.axes()
    plt.plot(x, np.abs(x), '-', color=green)
    plt.plot(0, 0, 'o', color=pink)

    # ax3.annotate('$s(x)=-0.02988 x^5 + 0.417 x^4 - 2.018 x^3 + 3.694 x^2 - 1.722 x - 5.511e^{-14}$', (-0.2, -1.2), size=13, color=orange)
    ax.annotate('Basic Function', (-1.9, 2.59), size=30, color=green)
    ax.annotate('Center', (-0.9, -0.5), size=30, color=pink)

    jesseaxis(ax,x,np.abs(x))
    plt.ylim(-1, 3)
    saveimg('fig4')
    plt.show()

def makefig5():
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
    ax1.annotate('$\psi_i$', (4,1), size=30, color=green)

    plt.plot(x, np.abs(x-x[10])-buf, '-', color=green)
    # ax1.annotate('$s(x)$', (4, 9), size=30, color=blue)
    jesseaxis(ax1,x,y)
    # Defining our interpolation function from the data sites
    # rbf = Rbf(x, y, epsilon=0.2)
    # Applying the interpolation
    # xi = np.linspace(-6, 6, 100)
    # yi = rbf(xi)
    # plt.plot(xi, yi, color=blue)
    saveimg('fig5')
    plt.show()
makefig5()
