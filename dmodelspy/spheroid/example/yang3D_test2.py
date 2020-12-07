# Generated with SMOP  0.41-beta
# from smop.libsmop import *
# yang3D_test2.m
from numpy import arange, zeros, dot
import matplotlib.pyplot as plt
import pandas as pd
import sys
sys.path.append('..')
from yang import yang
"""
Test of the 3D finite spheroidal source
TEST2 - tilted spheroid, free surface
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
A=0.5
P_G=0.1
mu=9600000000.0
nu=0.25
phi=0

x=(arange(- 30000,30000,500))
y=x.copy()
z=zeros(x.shape)
csv_kwargs = {'delim_whitespace': True, 'names': ['X', 'Y'], 'comment': '%'}

for theta, Ttxt in zip([90, 60, 30, 0.01],
                       ['90', '60', '30', '00']):
    u,v,w,*_=yang(x0,y0,z0,a,A,P_G,mu,nu,theta,phi,x,y,z)
    # plot displacement
    # figure('Position',concat([10,scrsz(4) / 5,dot(2,scrsz(3)) / 3,scrsz(4) / 3]))
    XLim=10
    fig, axs = plt.subplots(1, 3, sharex=True, sharey=True)

    U = pd.read_csv(f'FEM/Tilted spheroid/u_A05_t{Ttxt}_p00.txt', **csv_kwargs)
    ax = axs[0]
    ax.plot(0.001*x,u,'r.', label='analytical')
    ax.plot(0.001*U['X'],U['Y'],'b', label='FEM')
    ax.set_title('U')
    ax.set_xlabel('X, in kilometers')
    ax.set_ylabel('Displacement, in meters')
    ax.set_xlim([-XLim, XLim])
    ax.legend(frameon=False)

    V = pd.read_csv(f'FEM/Tilted spheroid/v_A05_t{Ttxt}_p00.txt', **csv_kwargs)
    ax = axs[1]
    ax.plot(dot(0.001,x),v,'r.', label='analytical')
    ax.plot(0.001*V['X'],V['Y'],'b', label='FEM')
    ax.set_title('V')
    ax.set_xlabel('X, in kilometers')
    ax.set_xlim([-XLim, XLim])
    ax.legend(frameon=False)

    W = pd.read_csv(f'FEM/Tilted spheroid/w_A05_t{Ttxt}_p00.txt', **csv_kwargs)
    ax = axs[2]
    ax.plot(dot(0.001,x),w,'r.', label='analytical')
    ax.plot(dot(0.001,W['X']),W['Y'],'b', label='FEM')
    ax.set_title('W')
    ax.set_xlabel('X, in kilometers')
    ax.set_xlim([-XLim,XLim])
    #ax.legend(frameon=False)
    
    fig.suptitle(f'z=0, $\\theta = {theta:g}^\circ $, $\phi = 0^\circ $, A=0.5')
# theta=90
# 
# u,v,w=yang(x0,y0,z0,a,A,P_G,mu,nu,theta,phi,x,y,z,nargout=3)
# # plot displacement
# figure('Position',concat([10,scrsz(4) / 5,dot(2,scrsz(3)) / 3,scrsz(4) / 3]))
# XLim=10
# X,U=textread('FEM/Tilted spheroid/u_A05_t90_p00.txt','%f %f','headerlines',1,nargout=2)
# 
# subplot(1,3,1)
# plot(dot(0.001,x),u,'r.',dot(0.001,X),U,'b')
# title('z=0, \theta=90\deg, \phi=0\deg, A=0.5')
# xlabel('X (km)')
# ylabel('U (m)')
# xlim(concat([- XLim,XLim]))
# #legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');
# X,V=textread('FEM/Tilted spheroid/v_A05_t90_p00.txt','%f %f','headerlines',1,nargout=2)
# 
# subplot(1,3,2)
# plot(dot(0.001,x),v,'r.',dot(0.001,X),V,'b')
# title('z=0, \theta=90\deg, \phi=0\deg, A=0.5')
# xlabel('X (km)')
# ylabel('V (m)')
# xlim(concat([- XLim,XLim]))
# legend('analytical','FEM')
# legend('Location','Best')
# legend('boxoff')
# X,W=textread('FEM/Tilted spheroid/w_A05_t90_p00.txt','%f %f','headerlines',1,nargout=2)
# 
# subplot(1,3,3)
# plot(dot(0.001,x),w,'r.',dot(0.001,X),W,'b')
# title('z=0, \theta=90\deg, \phi=0\deg, A=0.5')
# xlabel('X (km)')
# ylabel('W (m)')
# xlim(concat([- XLim,XLim]))
# #legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');
# 
# theta=60
# 
# u,v,w=yang(x0,y0,z0,a,A,P_G,mu,nu,theta,phi,x,y,z,nargout=3)
# # plot displacement
# figure('Position',concat([10,scrsz(4) / 5,dot(2,scrsz(3)) / 3,scrsz(4) / 3]))
# XLim=10
# X,U=textread('FEM/Tilted spheroid/u_A05_t60_p00.txt','%f %f','headerlines',1,nargout=2)
# 
# subplot(1,3,1)
# plot(dot(0.001,x),u,'r.',dot(0.001,X),U,'b')
# title('z=0, \theta=60\deg, \phi=0\deg, A=0.5')
# xlabel('X (km)')
# ylabel('U (m)')
# xlim(concat([- XLim,XLim]))
# #legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');
# X,V=textread('FEM/Tilted spheroid/v_A05_t60_p00.txt','%f %f','headerlines',1,nargout=2)
# 
# subplot(1,3,2)
# plot(dot(0.001,x),v,'r.',dot(0.001,X),V,'b')
# title('z=0, \theta=60\deg, \phi=0\deg, A=0.5')
# xlabel('X (km)')
# ylabel('V (m)')
# xlim(concat([- XLim,XLim]))
# legend('analytical','FEM')
# legend('Location','Best')
# legend('boxoff')
# X,W=textread('FEM/Tilted spheroid/w_A05_t60_p00.txt','%f %f','headerlines',1,nargout=2)
# 
# subplot(1,3,3)
# plot(dot(0.001,x),w,'r.',dot(0.001,X),W,'b')
# title('z=0, \theta=60\deg, \phi=0\deg, A=0.5')
# xlabel('X (km)')
# ylabel('W (m)')
# xlim(concat([- XLim,XLim]))
# #legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');
# 
# theta=30
# 
# u,v,w=yang(x0,y0,z0,a,A,P_G,mu,nu,theta,phi,x,y,z,nargout=3)
# # plot displacement
# figure('Position',concat([10,scrsz(4) / 5,dot(2,scrsz(3)) / 3,scrsz(4) / 3]))
# XLim=10
# X,U=textread('FEM/Tilted spheroid/u_A05_t30_p00.txt','%f %f','headerlines',1,nargout=2)
# 
# subplot(1,3,1)
# plot(dot(0.001,x),u,'r.',dot(0.001,X),U,'b')
# title('z=0, \theta=30\deg, \phi=0\deg, A=0.50')
# xlabel('X (km)')
# ylabel('U (m)')
# xlim(concat([- XLim,XLim]))
# #legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');
# X,V=textread('FEM/Tilted spheroid/v_A05_t30_p00.txt','%f %f','headerlines',1,nargout=2)
# 
# subplot(1,3,2)
# plot(dot(0.001,x),v,'r.',dot(0.001,X),V,'b')
# title('z=0, \theta=30\deg, \phi=0\deg, A=0.50')
# xlabel('X (km)')
# ylabel('V (m)')
# xlim(concat([- XLim,XLim]))
# legend('analytical','FEM')
# legend('Location','Best')
# legend('boxoff')
# X,W=textread('FEM/Tilted spheroid/w_A05_t30_p00.txt','%f %f','headerlines',1,nargout=2)
# 
# subplot(1,3,3)
# plot(dot(0.001,x),w,'r.',dot(0.001,X),W,'b')
# title('z=0, \theta=30\deg, \phi=0\deg, A=0.50')
# xlabel('X (km)')
# ylabel('W (m)')
# xlim(concat([- XLim,XLim]))
# #legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');
# fclose('all')
# disp('That's all folks!')
# theta=0.01
# 
# u,v,w=yang(x0,y0,z0,a,A,P_G,mu,nu,theta,phi,x,y,z,nargout=3)
# # plot displacement
# figure('Position',concat([10,scrsz(4) / 5,dot(2,scrsz(3)) / 3,scrsz(4) / 3]))
# XLim=10
# X,U=textread('FEM/Tilted spheroid/u_A05_t00_p00.txt','%f %f','headerlines',1,nargout=2)
# 
# subplot(1,3,1)
# plot(dot(0.001,x),u,'r.',dot(0.001,X),U,'b')
# title('z=0, \theta=0\deg, \phi=0\deg, A=0.5')
# xlabel('X (km)')
# ylabel('U (m)')
# xlim(concat([- XLim,XLim]))
# #legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');
# X,V=textread('FEM/Tilted spheroid/v_A05_t00_p00.txt','%f %f','headerlines',1,nargout=2)
# 
# subplot(1,3,2)
# plot(dot(0.001,x),v,'r.',dot(0.001,X),V,'b')
# title('z=0, \theta=0\deg, \phi=0\deg, A=0.5')
# xlabel('X (km)')
# ylabel('V (m)')
# xlim(concat([- XLim,XLim]))
# legend('analytical','FEM')
# legend('Location','Best')
# legend('boxoff')
# X,W=textread('FEM/Tilted spheroid/w_A05_t00_p00.txt','%f %f','headerlines',1,nargout=2)
# 
# subplot(1,3,3)
# plot(dot(0.001,x),w,'r.',dot(0.001,X),W,'b')
# title('z=0, \theta=0\deg, \phi=0\deg, A=0.5')
# xlabel('X (km)')
# ylabel('W (m)')
# xlim(concat([- XLim,XLim]))
# #legend('analytical','FEM'); legend('Location','Best'); legend('boxoff');

plt.show()
print("That's all folks!")