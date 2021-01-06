"""
Test 3D sill-like source, different depths
"""    
from numpy import sqrt, arange, dot, finfo, multiply
from scipy.interpolate import CubicSpline as spline
import matplotlib.pyplot as plt
import pandas as pd

from dmodelspy import Sill


eps = finfo('float').eps

# penny-crack parameters
sill = Sill(0, 0, 1000, P_G=0.01, a=1000., nu=0.25)
# print(sill)

# profile
x = 1000.0 *(arange(eps, 5, 0.05))
xlim = 10
y = dot(2, x)
r = sqrt(x**2 + y**2)

csv_kwargs = {'delim_whitespace': True, 'names': ['X', 'Y'], 'comment': '%'}

def plot_one(ax, xa, ya, xfem, yfem, ax_title):
    ax.plot(xa, ya, 'r.', label='analytical')
    ax.plot(xfem, yfem, 'b', label='FEM')
    ax.set_title(ax_title)
    
# Plot displacements for two different observer depths
for z, Ztxt in zip([0, 500], ['', '_500m']):
    fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)
    u, v, w = sill.calc_displ(x, y, z)
    U = pd.read_csv(f'FEM_SILL/sill_ur_r{Ztxt}.txt', **csv_kwargs)
    W = pd.read_csv(f'FEM_SILL/sill_w_r{Ztxt}.txt', **csv_kwargs)
    cs=spline(U['X'],U['Y'])
    UR=cs(r)

    plot_one(axs[0], dot(0.001, r), u, dot(0.001, r), multiply(UR, x) / r, 'U')
    plot_one(axs[1], dot(0.001, r), v, dot(0.001, r), multiply(UR, y) / r, 'V')
    plot_one(axs[2], dot(0.001, r), w, dot(0.001, W['X']), W['Y'], 'W')

    axs[0].set_ylabel('displacement, in meters')
    axs[0].set_xlim([0, xlim])
    axs[1].set_xlabel('radial distance, in kilometers')
    axs[2].legend(frameon=False)

    fig.suptitle(f'SILL-LIKE SOURCE: $z$={z}m, $a$={sill.a}m, $z_0$={sill.z0}m')
plt.show()

# Plot ground tilt
dwdx, dwdy, *_ = sill.calc_tilt(x, y, 0)
# read FEM models for radial and vertical deformations
DWDX = pd.read_csv('FEM_SILL/sill_dwdx_x.txt', **csv_kwargs)
DWDY = pd.read_csv('FEM_SILL/sill_dwdy_y.txt', **csv_kwargs)

fig, axs = plt.subplots(1, 2, sharex=True, sharey=True)
plot_one(axs[0], 0.001*x, dwdx, 0.001*DWDX['X'], DWDX['Y'], 'North')
plot_one(axs[1], 0.001*y, dwdy, 0.001*DWDY['X'], DWDY['Y'], 'East')

axs[0].set_xlabel('X, in kilometers')
axs[0].set_ylabel('Tilt')
axs[0].set_xlim([0, xlim])
axs[1].set_xlabel('Y, in kilometers')
axs[1].legend(frameon=False)
fig.suptitle(f'SILL-LIKE SOURCE: $z$=0, $a$={sill.a} m, $z_0$={sill.z0} m')
fig.tight_layout()

# Plot strain
eea, g1, g2 = sill.calc_strain(x, y, 500)
EEA = pd.read_csv('FEM_SILL/sill_eea_r_500m.txt', **csv_kwargs)
G1 = pd.read_csv('FEM_SILL/sill_gamma1_r_500m.txt', **csv_kwargs)
G2 = pd.read_csv('FEM_SILL/sill_gamma2_r_500m.txt', **csv_kwargs)

fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)
plot_one(axs[0], 0.001*r, eea, 0.001*EEA['X'], EEA['Y'], '$e_a$')
plot_one(axs[1], 0.001*r, g1, 0.001*G1['X'], G1['Y'], '$\\gamma_1$')
plot_one(axs[2], 0.001*r, g2, 0.001*G2['X'], G2['Y'], '$\\gamma_2$')

axs[0].set_ylabel('Strain')
axs[0].set_xlim([0, xlim])
axs[1].set_xlabel('Radial distance, in kilometers')
axs[2].legend(frameon=False)
fig.suptitle(f'SILL-LIKE SOURCE: $z$=500 m, $a$={sill.a} m, $z_0$={sill.z0} m')
fig.tight_layout()

plt.show()

print("That's all folks!")