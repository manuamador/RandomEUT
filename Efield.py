#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Functions to compute the radiation of a set of Hertzian dipoles
"""
from __future__ import division
from numpy import *

c = 299792458.0
  

def Efield(x,y,z,I,f):
    """
    Computes the E-field radiated by a set of Hertzian dipoles
    Parameters:
    - x,y,z, rectangular coordinates of the measuring point
    - I Dipole matrix (x|y|z|tilt|azimut|amplitude|phase)
    - f frequencies
    """
    N=len(x)
    Ex=zeros((N,len(f)))
    Ey=zeros((N,len(f)))
    Ez=zeros((N,len(f)))
    #Erac=zeros((N,len(f)))
    for i in range(0,N): 
        X=R*cos(phi[i])*sin(theta[i])
        Y=R*sin(phi[i])*sin(theta[i])
        Z=R*cos(theta[i])
        DX = X-I[:,0]
        DY = Y-I[:,1]
        DZ = Z-I[:,2]
        dist = sqrt(DX**2+DY**2+DZ**2)
        dp=tile(dist, (len(f),1))
        fp=tile(f,(len(dist),1))
        phaseI=tile(I[:,6],(len(f),1))
        phase=2*pi*dp*fp.T/c+phaseI
        ca    = cos(I[:,3])
        sa    = sin(I[:,3])
        cb    = cos(I[:,4])
        sb    = sin(I[:,4])
        distx = ((-sb)**2+(1-(-sb)**2)*ca)*DX+(-sb*cb*(1-ca))*DY\
            +(cb*sa)*DZ
        disty = (-sb*cb*(1-ca))*DX+((cb)**2+(1-cb**2)*ca)*DY+(sb*sa)*DZ
        distz = (-cb*sa)*DX+(-sb*sa)*DY+ca*DZ
        distxy = sqrt(distx**2+disty**2)
        costheta = distz/dist
        sintheta = distxy/dist
        cosphi   = distx/distxy
        sinphi   = disty/distxy
        L =tile(I[:,5],(len(f),1))*1/dp*(fp.T/c)**2*377
        Ex = sum(exp(1j*phase)*L*tile(((((-sb)**2+(1-(-sb)**2)*ca)\
            *(-sintheta*costheta*cosphi)+(-sb*cb*(1-ca))\
            *(-sintheta*costheta*sinphi)+(-cb*sa)\
            *(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
        Ey = sum(exp(1j*phase)*L*tile((((-sb*cb*(1-ca))\
            *(-sintheta*costheta*cosphi)+((cb)**2+(1-(cb)**2)*ca)\
            *(-sintheta*costheta*sinphi)+(-sb*sa)\
            *(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
        Ez = sum(exp(1j*phase)*L*tile((((cb*sa)\
            *(-sintheta*costheta*cosphi)+(sb*sa)\
            *(-sintheta*costheta*sinphi)\
            +ca*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)           
    return Ex,Ey,Ez


def Efarfield(R,theta,phi,I,f):
    """
    Computes the E-field radiated by a set of Hertzian dipoles
    Parameters:
    - R measurement distance
    - x,y,z, rectangular coordinates of the measuring point
    - I Dipole matrix (x|y|z|tilt|azimut|amplitude|phase)
    - f frequencies
    """
    N=len(theta)
    Eth=zeros((N,len(f)))
    Eph=zeros((N,len(f)))
    #Er=zeros((N,len(f)))
    for i in range(0,N): 
        X=R*cos(phi[i])*sin(theta[i])
        Y=R*sin(phi[i])*sin(theta[i])
        Z=R*cos(theta[i])
        DX = X-I[:,0]
        DY = Y-I[:,1]
        DZ = Z-I[:,2]
        dist = sqrt(DX**2+DY**2+DZ**2)
        dp=tile(dist, (len(f),1))
        fp=tile(f,(len(dist),1))
        phaseI=tile(I[:,6],(len(f),1))
        phase=2*pi*dp*fp.T/c+phaseI
        ca    = cos(I[:,3])
        sa    = sin(I[:,3])
        cb    = cos(I[:,4])
        sb    = sin(I[:,4])
        distx = ((-sb)**2+(1-(-sb)**2)*ca)*DX+(-sb*cb*(1-ca))*DY\
            +(cb*sa)*DZ
        disty = (-sb*cb*(1-ca))*DX+((cb)**2+(1-cb**2)*ca)*DY+(sb*sa)*DZ
        distz = (-cb*sa)*DX+(-sb*sa)*DY+ca*DZ
        distxy = sqrt(distx**2+disty**2)
        costheta = distz/dist
        sintheta = distxy/dist
        cosphi   = distx/distxy
        sinphi   = disty/distxy
        L =tile(I[:,5],(len(f),1))*1/dp*(fp.T/c)**2*377
        Exx = sum(exp(1j*phase)*L*tile(((((-sb)**2+(1-(-sb)**2)*ca)\
            *(-sintheta*costheta*cosphi)+(-sb*cb*(1-ca))\
            *(-sintheta*costheta*sinphi)+(-cb*sa)\
            *(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
        Eyy = sum(exp(1j*phase)*L*tile((((-sb*cb*(1-ca))\
            *(-sintheta*costheta*cosphi)+((cb)**2+(1-(cb)**2)*ca)\
            *(-sintheta*costheta*sinphi)+(-sb*sa)\
            *(-sintheta*(-sintheta)))),(len(f),1)),axis=1)
        Ezz = sum(exp(1j*phase)*L*tile((((cb*sa)\
            *(-sintheta*costheta*cosphi)+(sb*sa)\
            *(-sintheta*costheta*sinphi)\
            +ca*(-sintheta*(-sintheta)))),(len(f),1)),axis=1)           
        Eth[i,:]= abs(Exx*cos(theta[i])*cos(phi[i])\
            +Eyy*cos(theta[i])*sin(phi[i])-Ezz*sin(theta[i]))
        Eph[i,:]= abs(-Exx*sin(phi[i])+Eyy*cos(phi[i]))
        #Er[i,:] = abs(Exx*sin(theta[i])*cos(phi[i])\
        #   +Eyy*sin(theta[i])*sin(phi[i])+Ezz*cos(theta[i]))
    return Eth,Eph