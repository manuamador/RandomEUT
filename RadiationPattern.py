#!/usr/bin/env python
# -*- coding: utf8 -*-
"""
This program computes the radiation of a spherical equipment under test (EUT) with n arbitrary Hertzian dipoles on its surface
- The directivity of the EUT is computed
- The radiation pattern is printed in a picture file
"""

from __future__ import division
from numpy import *
from numpy.random import *
from pylab import *
from pylab import rcParams
from mpl_toolkits.mplot3d import Axes3D
import os
from Efield import Efarfield

c = 299792458.0

a=1. #EUT radius in m
f= [10e6,100e6,1e9]#Frequency in Hz
#ka=2*pi*f/c*a #electric size, ka=2*pi/lambda*a

#Radiation pattern
np=360  #number of points along phi
nt=180  #number of points along theta

phi=linspace(0,2*pi,np)
theta=arccos(2*linspace(0,1,nt)-1)#linspace(0,pi,nt)
TH,PH=meshgrid(theta,phi)
t=reshape(TH,(np*nt))
p=reshape(PH,(np*nt))
R = 100 #Measurement sphere radius in m



Ethac=zeros((len(phi),len(theta),len(f)),complex) #$E_\theta$
Ephac=zeros((len(phi),len(theta),len(f)),complex) #$E_\phi$



n=10 #number of radiating dipoles on the EUT
#generate the dipoles
theta_eut=arccos(2*rand(n,1)-1) #uniformly random along theta
phi_eut=2*pi*rand(n,1) #uniformly random along phi
x=a*cos(phi_eut)*sin(theta_eut)  
y=a*sin(phi_eut)*sin(theta_eut) 
z=a*cos(theta_eut) 
tilt=arccos(2*rand(n,1)-1)
azimut=2*pi*rand(n,1)
ld=.1 #taille des dipôles en m     
amplitude=ones((n,1))*ld #random amplitude of the currents
phas=2*pi*rand(n,1) #random phase

#Currents matrix
I=concatenate((x,y,z,tilt,azimut,amplitude,phas), axis=1)

#Efield computation
Etheta,Ephi=Efarfield(R,t,p,I,f)
Ethac=reshape(Etheta,(np,nt,len(f)))
Ephac=reshape(Ephi,(np,nt,len(f)))

P=real(Ethac)**2+real(Ephac)**2 #power

#Radiation diagram   
for u in range(0,len(f)):
    D=(P[:,:,u]).max()/(P[:,:,u]).mean() #Directivity $D_{\max}$
    fig = figure(figsize=(12, 12), dpi=50)
    ax = fig.add_subplot(111, projection='3d', frame_on=False)
    ax._axis3don = False
    R = P[:,:,u]/(P[:,:,u]).max()
    x = R  * outer(cos(phi), sin(theta))
    y = R  * outer(sin(phi), sin(theta))
    z = R  * outer(ones_like(phi), cos(theta))
    ##B&W 3D rendering    
    ax.plot_surface(x, y, z,  rstride=2, cstride=2,color='w',\
        linewidth=0.6,shade=False)
    ##Color 3D rendering, normalized power
    #ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=cm.jet(R),\
    #    linewidth=1,shade=True,antialiased=False)
    max_radius = 0.7
    print('f= %2.3f GHz, D = %2.1f' %(f[u]/1e9,D))
    title(r'$f= %2.3f$ GHz, $D = %2.1f$' %(f[u]/1e9,D),fontsize=20)
    for axis in 'xyz':
        getattr(ax, 'set_{}lim'.format(axis))((-max_radius, max_radius))
    fname = 'f%2.3f_%2.1f' %(f[u]/1e9,D)
    print 'Saving frame', fname
    fig.savefig(fname+'.png',bbox='tight')
    #fig.savefig(fname+'.svg',bbox='tight')
    #fig.savefig(fname+'.pdf',bbox='tight')
    close()
