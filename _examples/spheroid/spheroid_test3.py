"""
Test of the 3D finite spheroidal source

TEST3 - rotated spheroid, free surface
"""
from numpy import arange, zeros, dot
import matplotlib.pyplot as plt
import pandas as pd

from dmodelspy import spheroid


# ellipsoid parameters
x0=0
y0=0
z0=3000

a=1000
A=0.5
P_G=0.1
mu=9600000000.0
nu=0.25
theta=0.01

x=(arange(- 30000,30000,500))
y=x.copy()
z=zeros(x.shape)
csv_kwargs = {'delim_whitespace': True, 'names': ['X', 'Y'], 'comment': '%'}

for phi, Ptxt in zip([90., 60., 30., 0.],
                     ['90', '60', '30', '00']):
    u, v, w, *_ = spheroid(x0, y0, z0, a, A, P_G, mu, nu, theta, phi, x, y, z)
    # plot displacement
    XLim=10
    fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)

    U = pd.read_csv(f'FEM/Rotated spheroid/u_A05_t00_p{Ptxt}.txt', **csv_kwargs)
    ax = axs[0]
    ax.plot(0.001*x,u,'r.', label='analytical')
    ax.plot(0.001*U['X'],U['Y'],'b', label='FEM')
    ax.set_title('U')
    ax.set_ylabel('Displacement, in meters')
    ax.set_xlim([-XLim, XLim])
    ax.legend(frameon=False)

    V = pd.read_csv(f'FEM/Rotated spheroid/v_A05_t00_p{Ptxt}.txt', **csv_kwargs)
    ax = axs[1]
    ax.plot(dot(0.001,x),v,'r.', label='analytical')
    ax.plot(0.001*V['X'],V['Y'],'b', label='FEM')
    ax.set_title('V')
    ax.set_xlabel('X, in kilometers')
    ax.set_xlim([-XLim, XLim])
    ax.legend(frameon=False)

    W = pd.read_csv(f'FEM/Rotated spheroid/w_A05_t00_p{Ptxt}.txt', **csv_kwargs)
    ax = axs[2]
    ax.plot(dot(0.001,x),w,'r.', label='analytical')
    ax.plot(dot(0.001,W['X']),W['Y'],'b', label='FEM')
    ax.set_title('W')
    ax.set_xlim([-XLim,XLim])
    #ax.legend(frameon=False)
    
    fig.suptitle(f'z=0, $\\theta = 0^\circ $, $\phi = {phi:g}^\circ $, A=0.5')

plt.show()
print("That's all folks!")