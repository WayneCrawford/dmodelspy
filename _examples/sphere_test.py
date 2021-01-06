"""
Test of the 2D finite spherical source
"""
import matplotlib.pyplot as plt
import pandas as pd
from numpy import arange

from dmodelspy import Sphere2D, Sphere3D


# read input data for flat half-space model (McTigue, 1987)
with open('inputdata_sphere.txt') as f:
    z0 = float(f.readline().split()[-1])
    P_G = float(f.readline().split()[-1])
    a = float(f.readline().split()[-1])
    nu = float(f.readline().split()[-1])

# load FEM model
kwargs = {'delim_whitespace': True, 'header': 0, 'names': ['X', 'Y']}
U = pd.read_csv('FEM_SPHERE/udispl-0m.txt', **kwargs)
V = pd.read_csv('FEM_SPHERE/vdispl-0m.txt', **kwargs)
W = pd.read_csv('FEM_SPHERE/wdispl-0m.txt', **kwargs)

# DISPLACEMENT ************************************************************
# compute deformation at the surface
z = 0
x = (arange(1e-07, 5000, 100))
y = x.copy()
sphere3D = Sphere3D(0, 0, z0, P_G, a, nu)
u3D, v3D, w3D = sphere3D.calc_displ(x, y, z)
dwdx3D, dwdy3D = sphere3D.calc_tilt(x, y, z)
# print(f'{u3.shape=}, {v3.shape=}, {w3.shape=}')
sphere2D = Sphere2D(0, 0, z0, P_G, a, nu)
u2D, v2D, w2D = sphere2D.calc_displ(x, y)
dwdx2D, dwdy2D = sphere2D.calc_tilt(x, y)
# print(f'{u.shape=}, {v.shape=}, {w.shape=}, {dwdx.shape=}, {dwdy.shape=}')

def plot_one(ax, x2D, y2D, x3D, y3D, xfem, yfem, ax_title):
    ax.plot(x2D, y2D, 'r.', label='2D')
    ax.plot(x3D, y3D, 'g--', label='3D')
    ax.plot(xfem, yfem, 'b', label='FEM')
    ax.set_title(ax_title)
    
# compare 3D, 2D and FEM models
fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)
plot_one(axs[0], x, u2D, x, u3D, U['X'], U['Y'], 'X')
plot_one(axs[1], x, v2D, x, v3D, V['X'], V['Y'], 'Y')
plot_one(axs[2], x, w2D, x, w3D, W['X'], W['Y'], 'Z')

axs[0].set_ylabel('displacement (m) at z=0')
axs[1].set_xlabel('x (m)')
axs[2].legend(frameon=False)
plt.tight_layout()

# GROUND TILT *************************************************************
# load FE model
DX = pd.read_csv('FEM_SPHERE/dwdx-0m.txt', delim_whitespace=True, header=0,
                 names=['X', 'Y'])
DY = pd.read_csv('FEM_SPHERE/dwdy-0m.txt', delim_whitespace=True, header=0,
                 names=['X', 'Y'])

# compare ground tilt against FE model
fig, axs = plt.subplots(1, 2, sharex=True, sharey=True)
plot_one(axs[0], x, dwdx2D, x, dwdx3D, DX['X'], DX['Y'], 'E')
plot_one(axs[1], x, dwdy2D, x, dwdy3D, DY['X'], DY['Y'], 'N')

axs[0].set_xlabel('x (m)')
axs[0].set_ylabel('ground tilt')
axs[0].legend(frameon=False)
axs[1].set_xlabel('x (m)')
plt.tight_layout()

plt.show()
# filename = 'comparetiltFEM'; print('-djpeg',filename);

# PRINT TABLE *************************************************************
table = pd.DataFrame({'x_km': x/1.e3,
                      'y_km': y/1.e3,
                      'u_mm': 1.e3*u2D,
                      'v_mm': 1.e3*v2D,
                      'w_mm': 1.e3*w2D,
                      'dwdx_M': 1.e6*dwdx2D,
                      'dwdy_M': 1.e6*dwdy2D,
                      'u3_mm': 1.e3*u3D.flatten(),
                      'v3_mm': 1.e3*v3D.flatten(),
                      'w3_mm': 1.e3*w3D.flatten()})
table.to_csv('Table1.txt', index=False, float_format='%6.4f')