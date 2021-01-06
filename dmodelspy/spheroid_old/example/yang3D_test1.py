# Generated with SMOP  0.41-beta
# from smop.libsmop import *
# yang3D_test1.m
from numpy import arange, zeros, dot
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.append('..')
from yang import yang
"""
Test of the 3D finite spheroidal source
TEST1 - vertical spheroid, free surface
"""
# clear('all')
# close_('all')
# clc
# scrsz=get(0,'ScreenSize')
# ellipsoid parameters
x0=0
y0=0
z0=3000

a=1000
P_G=0.1
mu=9600000000.0
nu=0.25
theta=90
phi=0

x=(arange(- 30000,30000,500))
y=x.copy()
z=zeros(x.shape)
csv_kwargs = {'delim_whitespace': True, 'names': ['X', 'Y'], 'comment': '%'}

for A, Atxt in zip([0.99, 0.75, 0.5, 0.25],
                   ['100', '075', '050', '025']):
    u,v,w,*_=yang(x0,y0,z0,a,A,P_G,mu,nu,theta,phi,x,y,z)
    # plot displacement
    # figure('Position',concat([10,scrsz(4) / 5,dot(2,scrsz(3)) / 3,scrsz(4) / 3]))
    XLim=10
    fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)

    U = pd.read_csv(f'FEM/Vertical spheroid/u_A{Atxt}.txt', **csv_kwargs)
    ax = axs[0]
    ax.plot(0.001*x,u,'r.', label='analytical')
    ax.plot(0.001*U['X'],U['Y'],'b', label='FEM')
    ax.set_title('U')
    ax.set_xlabel('X, in kilometers')
    ax.set_ylabel('Displacement, in meters')
    ax.set_xlim([-XLim, XLim])
    ax.legend(frameon=False)

    V = pd.read_csv(f'FEM/Vertical spheroid/v_A{Atxt}.txt', **csv_kwargs)
    ax = axs[1]
    ax.plot(dot(0.001,x),v,'r.', label='analytical')
    ax.plot(0.001*V['X'],V['Y'],'b', label='FEM')
    ax.set_title('V')
    ax.set_xlabel('X, in kilometers')
    ax.set_xlim([-XLim, XLim])
    ax.legend(frameon=False)

    W = pd.read_csv(f'FEM/Vertical spheroid/w_A{Atxt}.txt', **csv_kwargs)
    ax = axs[2]
    ax.plot(dot(0.001,x),w,'r.', label='analytical')
    ax.plot(dot(0.001,W['X']),W['Y'],'b', label='FEM')
    ax.set_title('W')
    ax.set_xlabel('X, in kilometers')
    ax.set_xlim([-XLim,XLim])
    #ax.legend(frameon=False)
    
    fig.suptitle(f'z=0, $\\theta = 90^\circ $, $\phi = 0^\circ $, A={A:g}')

# A=0.99
# 
# u,v,w=yang(x0,y0,z0,a,A,P_G,mu,nu,theta,phi,x,y,z)
# # plot displacement
# # figure('Position',concat([10,scrsz(4) / 5,dot(2,scrsz(3)) / 3,scrsz(4) / 3]))
# XLim=10
# csv_kwargs = {'delim_whitespace': True, 'names': ['X', 'Y'], 'comment': '%'}
# 
# U = pd.read_csv('FEM/Vertical spheroid/u_A100.txt', **csv_kwargs)
# # X,U=textread('FEM/Vertical spheroid/u_A100.txt','%f %f','headerlines',1,nargout=2)
# fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)
# ax = axs[0]
# ax.plot(dot(0.001,x),u,'r.', label='analytical')
# ax.plot(dot(0.001,X),U,'b', label='FEM')
# ax.set_title('z=0, \theta=90\deg, \phi=0\deg, A=1')
# ax.set_xlabel('X, in kilometers')
# ax.set_ylabel('U, in meters')
# ax.set_xlim([-XLim, XLim])
# ax.legend(frameon=False)
# 
# V = pd.read_csv('FEM/Vertical spheroid/v_A100.txt', **csv_kwargs)
# ax = axs[1]
# ax.plot(dot(0.001,x),v,'r.', label='analytical')
# ax.plot(0.001*V['X'],V['Y'],'b', label='FEM')
# ax.set_title('z=0, \theta=90\deg, \phi=0\deg, A=1')
# ax.set_xlabel('X, in kilometers')
# ax.set_ylabel('V, in meters')
# ax.set_xlim([-XLim, XLim])
# ax.legend(frameon=False)
# 
# # X,W=textread('FEM/Vertical spheroid/w_A100.txt','%f %f','headerlines',1,nargout=2)
# W = pd.read_csv('FEM/Vertical spheroid/w_A100.txt', **csv_kwargs)
# ax = axs[2]
# ax.plot(dot(0.001,x),w,'r.', label='analytical')
# ax.plot(dot(0.001,X),W,'b', label='FEM')
# ax.set_title('z=0, \theta=90\deg, \phi=0\deg, A=1')
# ax.set_xlabel('X, in kilometers')
# ax.set_ylabel('W, in meters')
# ax.set_xlim([- XLim,XLim])
# #ax.legend(frameon=False)
# A=0.75
# 
# u,v,w=yang(x0,y0,z0,a,A,P_G,mu,nu,theta,phi,x,y,z)
# # plot displacement
# # figure('Position',concat([10,scrsz(4) / 5,dot(2,scrsz(3)) / 3,scrsz(4) / 3]))
# fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)
# XLim=10
# X,U=textread('FEM/Vertical spheroid/u_A075.txt','%f %f','headerlines',1,nargout=2)
# 
# ax = axs[0]
# ax.plot(dot(0.001,x),u,'r.', label='analytical')
# ax.plot(dot(0.001,X),U,'b', label='FEM')
# ax.set_title('z=0, \theta=90\deg, \phi=0\deg, A=0.75')
# ax.set_xlabel('X, in kilometers')
# ax.set_ylabel('U, in meters')
# ax.set_xlim(concat([- XLim,XLim]))
# # ax.legend(frameon=False)
# X,V=textread('FEM/Vertical spheroid/v_A075.txt','%f %f','headerlines',1,nargout=2)
# 
# ax = axs[1]
# ax.plot(dot(0.001,x),v,'r.', label='analytical')
# ax.plot(dot(0.001,X),V,'b', label='FEM')
# ax.set_title('z=0, \theta=90\deg, \phi=0\deg, A=0.75')
# ax.set_xlabel('X, in kilometers')
# ax.set_ylabel('V, in meters')
# ax.set_xlim(concat([- XLim,XLim]))
# ax.legend(frameon=False)
# X,W=textread('FEM/Vertical spheroid/w_A075.txt','%f %f','headerlines',1,nargout=2)
# 
# ax = axs[2]
# ax.plot(dot(0.001,x),w,'r.', label='analytical')
# ax.plot(dot(0.001,X),W,'b', label='FEM')
# ax.set_title('z=0, \theta=90\deg, \phi=0\deg, A=0.75')
# ax.set_xlabel('X, in kilometers')
# ax.set_ylabel('W, in meters')
# ax.set_xlim([- XLim,XLim])
# # ax.legend(frameon=False)
# 
# A=0.5
# 
# u,v,w=yang(x0,y0,z0,a,A,P_G,mu,nu,theta,phi,x,y,z,nargout=3)
# # plot displacement
# figure('Position',concat([10,scrsz(4) / 5,dot(2,scrsz(3)) / 3,scrsz(4) / 3]))
# XLim=10
# X,U=textread('FEM/Vertical spheroid/u_A050.txt','%f %f','headerlines',1,nargout=2)
# 
# ax = axs[0]
# ax.plot(dot(0.001,x),u,'r.',dot(0.001,X),U,'b')
# ax.title('z=0, \theta=90\deg, \phi=0\deg, A=0.50')
# ax.xlabel('X, in kilometers')
# ax.ylabel('U, in meters')
# ax.xlim(concat([- XLim,XLim]))
# # ax.legend(frameon=False)
# X,V=textread('FEM/Vertical spheroid/v_A050.txt','%f %f','headerlines',1,nargout=2)
# 
# ax = axs[1]
# ax.plot(dot(0.001,x),v,'r.',dot(0.001,X),V,'b')
# ax.title('z=0, \theta=90\deg, \phi=0\deg, A=0.50')
# ax.xlabel('X, in kilometers')
# ax.ylabel('V, in meters')
# ax.xlim(concat([- XLim,XLim]))
# ax.legend(frameon=False)
# X,W=textread('FEM/Vertical spheroid/w_A050.txt','%f %f','headerlines',1,nargout=2)
# 
# ax = axs[2]
# ax.plot(dot(0.001,x),w,'r.',dot(0.001,X),W,'b')
# ax.title('z=0, \theta=90\deg, \phi=0\deg, A=0.50')
# ax.xlabel('X, in kilometers')
# ax.ylabel('W, in meters')
# ax.xlim(concat([- XLim,XLim]))
# # ax.legend(frameon=False)
# # fclose('all')
# disp('That's all folks!')
# A=0.25
# 
# u,v,w=yang(x0,y0,z0,a,A,P_G,mu,nu,theta,phi,x,y,z,nargout=3)
# # plot displacement
# figure('Position',concat([10,scrsz(4) / 5,dot(2,scrsz(3)) / 3,scrsz(4) / 3]))
# XLim=10
# X,U=textread('FEM/Vertical spheroid/u_A025.txt','%f %f','headerlines',1,nargout=2)
# 
# ax = axs[0]
# ax.plot(dot(0.001,x),u,'r.',dot(0.001,X),U,'b')
# ax.title('z=0, \theta=90\deg, \phi=0\deg, A=0.25')
# ax.xlabel('X, in kilometers')
# ax.ylabel('U, in meters')
# ax.xlim(concat([- XLim,XLim]))
# # ax.legend(frameon=False)
# X,V=textread('FEM/Vertical spheroid/v_A025.txt','%f %f','headerlines',1,nargout=2)
# 
# ax = axs[1]
# ax.plot(dot(0.001,x),v,'r.',dot(0.001,X),V,'b')
# ax.title('z=0, \theta=90\deg, \phi=0\deg, A=0.25')
# ax.xlabel('X, in kilometers')
# ax.ylabel('V, in meters')
# ax.xlim(concat([- XLim,XLim]))
# ax.legend(frameon=False)
# X,W=textread('FEM/Vertical spheroid/w_A025.txt','%f %f','headerlines',1,nargout=2)
# 
# ax = axs[2]
# ax.plot(dot(0.001,x),w,'r.',dot(0.001,X),W,'b')
# ax.title('z=0, \theta=90\deg, \phi=0\deg, A=0.25')
# ax.xlabel('X, in kilometers')
# ax.ylabel('W, in meters')
# ax.xlim(concat([- XLim,XLim]))
# # ax.legend(frameon=False)

plt.show()
# fclose('all')
print("That's all folks!")