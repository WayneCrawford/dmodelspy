# Generated with SMOP  0.41-beta
from smop.libsmop import *
# mctigue2D_test.m
"""
Test of the 2D finite spherical source
"""
import sys

import matplotlib.pyplot as plt
import pandas as pd
# from numpy import arange

# sys.path.append('..')
from mctigue3Ddispl import mctigue3D as mctigue3Ddispl
from mctigue2D import mctigue2D

# clear('all')
# close_('all')
# clc
# scrsz=get(0,'ScreenSize')

# read input data for flat half-space model (McTigue, 1987)
with open('inputdata.txt') as f:
    z0 = float(f.readline().split()[-1])
    P_G = float(f.readline().split()[-1])
    a = float(f.readline().split()[-1])
    nu = float(f.readline().split()[-1])
# __, z0 = textread('inputdata.txt','%s %f',1,nargout=2)
# __,P_G=textread('inputdata.txt','%s %f',1,'headerlines',1,nargout=2)
# __,a=textread('inputdata.txt','%s %f',1,'headerlines',2,nargout=2)
# __,nu=textread('inputdata.txt','%s %f',1,'headerlines',3,nargout=2)

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
u3, v3, w3 = mctigue3Ddispl(0, 0, z0, P_G, a, nu, x, y, z, nargout=5)
u, v, w, dwdx, dwdy = mctigue2D(0, 0, z0, P_G, a, nu, x, y, nargout=5)

# compare 3D, 2D and FEM models
# figure('Position',concat([10,scrsz(4) / 5,dot(2,scrsz(3)) / 3,scrsz(4) / 3]))
# subplot(1,3,1)
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
ax.plot(x, w3)
ax.plot(x, w, 'r.', label='2D')
ax.plot(W['X'], W['Y'], 'g--', label='FEM')
ax.set_xlabel('x (m)')
ax.set_title('Z')
ax.legend()
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
ax.legend()
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
# fid=fopen('Table1.txt','w+')
table = pd.DataFrame([x/1.e3, y/1.e3,
                      1.e3*u, 1.e3*v, 1.e3*w,
                      1.e6*dwdx, 1.e6*dwdy,
                      1.e3*u3, 1.e3*v3, 1.e3*w3])
table.to_csv('Table1.txt')
# with open('Table1.txt', 'w+') as f:
# M=concat([x.T / 1000, y.T / 1000, dot(1000,u.T), dot(1000,v.T),
#           dot(1000,w.T), dot(1000000.0,dwdx.T), dot(1000000.0,dwdy.T),
#           dot(1000,u3.T), dot(1000,v3.T), dot(1000,w3.T)])
# fprintf(fid,'%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f,%6.4f
#  ,%6.4f\n',M.T)
