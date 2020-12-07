"""
Test of the 3D sill-like source

Test tilt and strain
"""
from numpy import sqrt, arange, dot, finfo
import matplotlib.pyplot as plt
import pandas as pd

from dmodels_py import sill


eps = finfo('float').eps

# penny-crack parameters
x0 = 0
y0 = 0
z0 = 1000
P_G = 0.01
a = 1000
nu = 0.25

# profile
# profile
x = 1000. * (arange(eps, 5, 0.05))
y = dot(2, x)
r = sqrt(x**2 + y**2)
XLim = 10
csv_kwargs = {'delim_whitespace': True, 'names': ['X', 'Y'], 'comment': '%'}
# TILT ********************************************************************
_, _, _, _, dwdx, dwdy, *_ = sill(x0, y0, z0, P_G, a, nu, x, y, 0)
# plot ground tilt
# figure('Position',concat([10,scrsz(4) / 5,dot(2,scrsz(3)) / 3,scrsz(4) / 3]))

# read FEM models for radial and vertical deformations
DWDX = pd.read_csv('FEM/sill_dwdx_x.txt', **csv_kwargs)
DWDY = pd.read_csv('FEM/sill_dwdy_y.txt', **csv_kwargs)
# X,DWDX=textread('FEM/sill_dwdx_x.txt','%f %f','headerlines',1,nargout=2)
# Y,DWDY=textread('FEM/sill_dwdy_y.txt','%f %f','headerlines',1,nargout=2)

fig, axs = plt.subplots(1, 2, sharex=True, sharey=True)

ax = axs[0]
ax.plot(0.001*x, dwdx, 'r.')
ax.plot(0.001*DWDX['X'], DWDX['Y'], 'b')
ax.set_title('North')
ax.set_xlabel('X, in kilometers')
ax.set_ylabel('Tilt')
ax.set_xlim([0, XLim])

ax = axs[1]
ax.plot(0.001*y, dwdy, 'r.', label='Analytical')
ax.plot(0.001*DWDY['X'], DWDY['Y'], 'b', label='FEM')
ax.set_title('East')
ax.set_xlabel('Y, in kilometers')
ax.legend(frameon=False)

fig.suptitle('SILL-LIKE SOURCE: $z$=0, $a$=1000 m, $z_0$=1000 m')
fig.tight_layout()

# STRAIN ******************************************************************
fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)
XLim = 10
EEA = pd.read_csv('FEM/sill_eea_r_500m.txt', **csv_kwargs)
GAMMA1 = pd.read_csv('FEM/sill_gamma1_r_500m.txt', **csv_kwargs)
GAMMA2 = pd.read_csv('FEM/sill_gamma2_r_500m.txt', **csv_kwargs)
_, _, _, _, _, _, eea, gamma1, gamma2 = sill(x0, y0, z0, P_G, a, nu,
                                                 x, y, 500)
# plot strain
# figure('Position',concat([10,scrsz(4)/5,dot(2,scrsz(3))/3,scrsz(4)/3]))
# R,EEA=textread('FEM/sill_eea_r_500m.txt','%f %f','headerlines',1)

ax = axs[0]
ax.plot(dot(0.001, r), eea, 'r.')
ax.plot(dot(0.001, EEA['X']), EEA['Y'], 'b')
ax.set_title('$e_a$')
ax.set_ylabel('Strain')
ax.set_xlim([0, XLim])

# R,GAMMA1=textread('FEM/sill_gamma1_r_500m.txt','%f %f','headerlines',1)
ax = axs[1]
ax.plot(dot(0.001, r), gamma1, 'r.')
ax.plot(dot(0.001, GAMMA1['X']), GAMMA1['Y'], 'b')
ax.set_title('$\\gamma_1$')
ax.set_xlabel('Radial distance, in kilometers')

# R,GAMMA2=textread('FEM/sill_gamma2_r_500m.txt','%f %f','headerlines',1)
ax = axs[2]
ax.plot(dot(0.001, r), gamma2, 'r.', label='analytical')
ax.plot(dot(0.001, GAMMA2['X']), GAMMA2['Y'], 'b', label='FEM')
ax.set_title('$\\gamma_2$')
ax.legend(frameon=False)

fig.suptitle('SILL-LIKE SOURCE: $z$=500 m, $a$=1000 m, $z_0$=1000 m')
fig.tight_layout()

plt.show()
print("That's all folks!")
