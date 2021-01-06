# Generated with SMOP  0.41-beta
# from smop.libsmop import *
# yang3D_test1.m
"""
Test of the 3D finite spheroidal source
"""
from numpy import arange, zeros, dot
import matplotlib.pyplot as plt
import pandas as pd

from dmodelspy import Spheroid


# ellipsoid parameters
model = Spheroid(0, 0, 3000, a=1000, asrat=0.5, P_G=0.1,
                 mu=9.6e9, nu=0.25, theta=90, phi=0)

x=arange(-30000, 30000, 500)
y = x.copy()
z = zeros(x.shape)
csv_kwargs = {'delim_whitespace': True, 'names': ['X', 'Y'], 'comment': '%'}
xlim = 10

def plot_one(ax, xa, ya, xfem, yfem, ax_title):
    ax.plot(xa, ya, 'r.', label='analytical')
    ax.plot(xfem, yfem, 'b', label='FEM')
    ax.set_title(ax_title)
    
# plot displacement
print('TEST1 - vertical spheroid, free surface, changing aspect ratio')
for asrat, Atxt in zip([0.99, 0.75, 0.5, 0.25],
                   ['100', '075', '050', '025']):
    model.asrat = asrat
    u, v, w = model.calc_displ(x,y,z)
    fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)

    U = pd.read_csv(f'FEM_SPHEROID/Vertical spheroid/u_A{Atxt}.txt', **csv_kwargs)
    V = pd.read_csv(f'FEM_SPHEROID/Vertical spheroid/v_A{Atxt}.txt', **csv_kwargs)
    W = pd.read_csv(f'FEM_SPHEROID/Vertical spheroid/w_A{Atxt}.txt', **csv_kwargs)

    plot_one(axs[0], 0.001*x, u, 0.001*U['X'], U['Y'], 'U')
    plot_one(axs[1], 0.001*x, v, 0.001*V['X'], V['Y'], 'V')
    plot_one(axs[2], 0.001*x, w, 0.001*W['X'], W['Y'], 'W')

    axs[0].set_ylabel('Displacement, in meters')
    axs[0].set_xlim([-xlim, xlim])
    axs[1].set_xlabel('X, in kilometers')
    axs[2].legend(frameon=False)
    
    fig.suptitle('Changing A (aspect ratio): z=0, $\\theta = {:g}^\circ $, $\phi = {:g}^\circ $, A={:g}'
                 .format(model.theta, model.phi, model.asrat))

plt.show()

print('TEST2 - tilted spheroid, free surface, changing tilt')
model.asrat = 0.5
for theta, Ttxt in zip([90, 60, 30, 0.01],
                       ['90', '60', '30', '00']):
    model.theta = theta
    u, v, w = model.calc_displ(x,y,z)
    # plot displacement
    fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)

    U = pd.read_csv(f'FEM_SPHEROID/Tilted spheroid/u_A05_t{Ttxt}_p00.txt', **csv_kwargs)
    V = pd.read_csv(f'FEM_SPHEROID/Tilted spheroid/v_A05_t{Ttxt}_p00.txt', **csv_kwargs)
    W = pd.read_csv(f'FEM_SPHEROID/Tilted spheroid/w_A05_t{Ttxt}_p00.txt', **csv_kwargs)

    plot_one(axs[0], 0.001*x, u, 0.001*U['X'], U['Y'], 'U')
    plot_one(axs[1], 0.001*x, v, 0.001*V['X'], V['Y'], 'V')
    plot_one(axs[2], 0.001*x, w, 0.001*W['X'], W['Y'], 'W')

    axs[0].set_ylabel('Displacement, in meters')
    axs[0].set_xlim([-xlim, xlim])
    axs[1].set_xlabel('X, in kilometers')
    axs[2].legend(frameon=False)
    
    fig.suptitle('Changing $\\theta$ (tilt): z=0, $\\theta = {:g}^\circ $, $\phi = {:g}^\circ $, A={:g}'
                 .format(model.theta, model.phi, model.asrat))

plt.show()

print('TEST3 - horizontal spheroid, free surface, changing rotation')
model.theta=0.01
for phi, Ptxt in zip([90., 60., 30., 0.],
                     ['90', '60', '30', '00']):
    model.phi = phi
    u, v, w = model.calc_displ(x,y,z)
    # plot displacement
    fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)

    U = pd.read_csv(f'FEM_SPHEROID/Rotated spheroid/u_A05_t00_p{Ptxt}.txt', **csv_kwargs)
    V = pd.read_csv(f'FEM_SPHEROID/Rotated spheroid/v_A05_t00_p{Ptxt}.txt', **csv_kwargs)
    W = pd.read_csv(f'FEM_SPHEROID/Rotated spheroid/w_A05_t00_p{Ptxt}.txt', **csv_kwargs)

    plot_one(axs[0], 0.001*x, u, 0.001*U['X'], U['Y'], 'U')
    plot_one(axs[1], 0.001*x, v, 0.001*V['X'], V['Y'], 'V')
    plot_one(axs[2], 0.001*x, w, 0.001*W['X'], W['Y'], 'W')

    axs[0].set_ylabel('Displacement, in meters')
    axs[0].set_xlim([-xlim, xlim])
    axs[1].set_xlabel('X, in kilometers')
    axs[2].legend(frameon=False)
    
    fig.suptitle('Changing $\\phi$ (rotation): z=0, $\\theta = {:g}^\circ $, $\phi = {:g}^\circ $, A={:g}'
                 .format(model.theta, model.phi, model.asrat))
plt.show()

print("That's all folks!")