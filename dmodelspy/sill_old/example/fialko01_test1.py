# Generated with SMOP  0.41-beta
# from smop.libsmop import *
# fialko01_test1.m
from numpy import sqrt, arange, dot, finfo
from scipy.interpolate import CubicSpline as spline
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.append('..')
from fialko01 import fialko01
"""
Test of the 3D sill-like source
"""    
# clear('all')
# close_('all')
# clc
# scrsz=get(0,'ScreenSize')
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
    u, v, w, dV, *_ = fialko01(x0, y0, z0, P_G, a, nu, x, y, z)
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
# 
# plot displacement at z=0;
# figure('Position',concat([10,scrsz(4) / 5,dot(2,scrsz(3)) / 3,scrsz(4) / 3]))
# u,v,w,dV=fialko01(x0,y0,z0,P_G,a,nu,x,y,0,nargout=4)
# rr,ur=textread('FEM/sill_ur_r.txt','%f %f','headerlines',1,nargout=2)
# rw,W=textread('FEM/sill_w_r.txt','%f %f','headerlines',1,nargout=2)
#
# # plot displacement at z=500 m;
# u,v,w,dV=fialko01(x0,y0,z0,P_G,a,nu,x,y,500,nargout=4)
# figure('Position',concat([10,scrsz(4) / 5,dot(2,scrsz(3)) / 3,scrsz(4) / 3]))
# XLim=15
# rr,ur=textread('FEM/sill_ur_r_500m.txt','%f %f','headerlines',1,nargout=2)
# rw,W=textread('FEM/sill_w_r_500m.txt','%f %f','headerlines',1,nargout=2)
# 
# UR=spline(rr,ur,r)
# subplot(1,3,1)
# plot(dot(0.001,r),u,'r.',dot(0.001,r),multiply(UR,x) / r,'b')
# title('SILL-LIKE SOURCE')
# xlabel('radial distance, in kilometers')
# ylabel('U, in meters')
# xlim(concat([0,XLim]))
# ylim(concat([0,Inf]))
# subplot(1,3,2)
# plot(dot(0.001,r),v,'r.',dot(0.001,r),multiply(UR,y) / r,'b')
# title('z=500 m, a=1000 m, z_0=1000 m')
# xlabel('radial distance, in kilometers')
# ylabel('V, in meters')
# xlim(concat([0,XLim]))
# ylim(concat([0,Inf]))
# subplot(1,3,3)
# plot(dot(0.001,r),w,'r.',dot(0.001,rw),W,'b')
# xlabel('radial distance, in kilometers')
# ylabel('W, in meters')
# xlim(concat([0,XLim]))
# ylim(concat([0,Inf]))
# legend('analytical','FEM')
# legend('Location','Best')
# legend('boxoff')
plt.show()
print("That's all folks!")