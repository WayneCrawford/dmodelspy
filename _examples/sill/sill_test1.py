"""
Test of the 3D sill-like source

Test different depths
"""    
from numpy import sqrt, arange, dot, finfo, multiply
from scipy.interpolate import CubicSpline as spline
import matplotlib.pyplot as plt
import pandas as pd

from dmodelspy import sill


eps = finfo('float').eps

# penny-crack parameters
x0 = 0
y0 = 0
z0 = 1000
P_G = 0.01
a = 1000
nu = 0.25

# profile
x = 1000.0 *(arange(eps, 5, 0.05))
y = dot(2, x)
r = sqrt(x**2 + y**2)
XLim = 10
csv_kwargs = {'delim_whitespace': True, 'names': ['X', 'Y'], 'comment': '%'}

for z, Ztxt in zip([0, 500], ['', '_500m']):
    fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)
    u, v, w, dV, *_ = sill(x0, y0, z0, P_G, a, nu, x, y, z)
    U = pd.read_csv(f'FEM/sill_ur_r{Ztxt}.txt', **csv_kwargs)
    W = pd.read_csv(f'FEM/sill_w_r{Ztxt}.txt', **csv_kwargs)
    cs=spline(U['X'],U['Y'])
    UR=cs(r)

    ax = axs[0]
    ax.plot(dot(0.001, r), u, 'r.', label='analytical')
    ax.plot(dot(0.001, r), multiply(UR, x) / r, 'b', label='FEM')
    ax.set_title('U')
    ax.set_ylabel('displacement, in meters')
    ax.set_xlim([0, XLim])

    ax = axs[1]
    ax.plot(dot(0.001, r), v, 'r.', label='analytical')
    ax.plot(dot(0.001, r), multiply(UR, y) / r, 'b', label='FEM')
    ax.set_title('V')
    ax.set_xlabel('radial distance, in kilometers')

    ax = axs[2]
    ax.plot(dot(0.001, r), w, 'r.', label='analytical')
    ax.plot(dot(0.001, W['X']), W['Y'], 'b', label='FEM')
    ax.set_title('W')
    ax.legend(frameon=False)

    fig.suptitle(f'SILL-LIKE SOURCE: $z$={z}m, $a$={a}m, $z_0$={z0}m')
plt.show()
print("That's all folks!")