"""
Test of the 2D finite spherical source
"""
import matplotlib.pyplot as plt
import pandas as pd
from numpy import arange

from dmodelspy import sphere2D, sphere3D_displ


# read input data for flat half-space model (McTigue, 1987)
with open('inputdata.txt') as f:
    z0 = float(f.readline().split()[-1])
    P_G = float(f.readline().split()[-1])
    a = float(f.readline().split()[-1])
    nu = float(f.readline().split()[-1])

# load FEM model
kwargs = {'delim_whitespace': True, 'header': 0, 'names': ['X', 'Y']}
U = pd.read_csv('FEM/udispl-0m.txt', **kwargs)
V = pd.read_csv('FEM/vdispl-0m.txt', **kwargs)
W = pd.read_csv('FEM/wdispl-0m.txt', **kwargs)

# DISPLACEMENT ************************************************************
# compute deformation at the surface
z = 0
x = (arange(1e-07, 5000, 100))
y = x.copy()
u3, v3, w3 = sphere3D_displ(0, 0, z0, P_G, a, nu, x, y, z)
# print(f'{u3.shape=}, {v3.shape=}, {w3.shape=}')
u, v, w, dwdx, dwdy = sphere2D(0, 0, z0, P_G, a, nu, x, y)
# print(f'{u.shape=}, {v.shape=}, {w.shape=}, {dwdx.shape=}, {dwdy.shape=}')

# compare 3D, 2D and FEM models
fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)
ax=axs[0]
ax.plot(x, u3)
ax.plot(x, u, 'r.')
ax.plot(U['X'], U['Y'], 'g--')
ax.set_xlabel('x (m)')
ax.set_ylabel('displacement (m) at z=0')
ax.set_title('X')

ax=axs[1]
ax.plot(x, v3)
ax.plot(x, v, 'r.')
ax.plot(V['X'], V['Y'], 'g--')
ax.set_xlabel('x (m)')
ax.set_title('Y')

ax=axs[2]
ax.plot(x, w3, label='3D')
ax.plot(x, w, 'r.', label='2D')
ax.plot(W['X'], W['Y'], 'g--', label='FEM')
ax.set_xlabel('x (m)')
ax.set_title('Z')
ax.legend(frameon=False)
plt.tight_layout()
# legend('boxoff')
# filename = 'compare2D3DFEM'; print('-djpeg',filename);

# GROUND TILT *************************************************************
# load FE model
# X1, dWdx=textread('FEM/dwdx-0m.txt','%f %f','headerlines',1,nargout=2)
# X2 ,dWdy=textread('FEM/dwdy-0m.txt','%f %f','headerlines',1,nargout=2)
DX = pd.read_csv('FEM/dwdx-0m.txt', delim_whitespace=True, header=0,
                 names=['X', 'Y'])
DY = pd.read_csv('FEM/dwdy-0m.txt', delim_whitespace=True, header=0,
                 names=['X', 'Y'])

# compare ground tilt against FE model
# figure('Position',concat([10,scrsz(4) / 5,scrsz(3) / 3,scrsz(4) / 4]))
# subplot(1,2,1)
fig, axs = plt.subplots(1, 2, sharex=True, sharey=True)
ax = axs[0]
ax.plot(x, dwdx, 'r.', label='2D')
ax.plot(DX['X'], DX['Y'], 'g--', label='FEM')
ax.set_xlabel('x (m)')
ax.set_ylabel('ground tilt')
ax.set_title('E')
ax.legend(frameon=False)
# legend('boxoff')

ax = axs[1]
ax.plot(x, dwdy, 'r.')
ax.plot(DY['X'], DY['Y'], 'g--')
ax.set_xlabel('x (m)')
ax.set_title('N')
plt.tight_layout()

plt.show()
# filename = 'comparetiltFEM'; print('-djpeg',filename);

# PRINT TABLE *************************************************************
table = pd.DataFrame({'x_km': x/1.e3,
                      'y_km': y/1.e3,
                      'u_mm': 1.e3*u,
                      'v_mm': 1.e3*v,
                      'w_mm': 1.e3*w,
                      'dwdx_M': 1.e6*dwdx,
                      'dwdy_M': 1.e6*dwdy,
                      'u3_mm': 1.e3*u3.flatten(),
                      'v3_mm': 1.e3*v3.flatten(),
                      'w3_mm': 1.e3*w3.flatten()})
table.to_csv('Table1.txt', index=False, float_format='%6.4f')