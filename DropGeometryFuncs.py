# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 09:11:40 2023

@author: laplaud

"""

import numpy as np

import VallapFunc_Sim as vf

####
# 1. Compute cone geometrical transformation parameters

""" For a specific cone, of diameter Rc and angle Alpha with the vertical, we compute :
the radius Ro of the original circle, 
the angle Beta of the sector removed from the circle, 
the fraction Fr of circle surface in the cone"""

def coneGeometry(Rcone,Alpha):
    # perimeter of cone 
    Pcone = 2*np.pi*Rcone # (in mm)

    # Cone Side length
    Lcone = Rcone/np.sin(Alpha)

    """ Unfolded cone as a cropped circle of radius Lcone """

    Rcircle = Lcone.copy()

    # Perimter of full circle
    Pcircle = 2*np.pi*Rcircle

    # Fraction of the circle remaining to form the cone
    Fr = Pcone/Pcircle # analytically Fr = sin(Alpha)

    # Angle of the removed circle part
    Beta = 2*np.pi*(1-Fr) # analytically Beta = 2*pi*(1-sin(Alpha))
        
    return(Rcircle,Beta,Fr)



####
# 2. Transformation of 2D coordinate from cone configuration to circle configuration and vice versa 

def Cone2Circle(X,Y,Alpha):
    # (X,Y) points to transform, Alpha angle of the cone,
    central = (Y == 0) & (X<0)
    A,R = vf.ToCirc(X,Y,angle = 'rad')
    if np.size(A)== 1:
        if central:
            A = np.abs(A)
    else:
        A[central] = np.abs(A[central])
    Xnew,Ynew = vf.ToCart((A)*np.sin(Alpha),R/np.sin(Alpha),angle = 'rad')

    return(Xnew,Ynew)


def Circle2Cone(X,Y,Alpha): 
    # (X,Y) points to transform, Alpha angle of the cone
    A,R = vf.ToCirc(X,Y,angle = 'rad')
    Xnew,Ynew = vf.ToCart((A)/np.sin(Alpha),R*np.sin(Alpha),angle = 'rad')
    
    return(Xnew,Ynew)


def Cone2CircleZ(X,Y,Z,Alpha):
    # (X,Y) points to transform, Alpha angle of the cone, Ad = 0 angle of removed sector bissecant (cone config)
    
    A,R = vf.ToCirc(X,Y,angle = 'rad') 
    
    R1 = Z*np.sin(Alpha)*np.cos(Alpha) + R*np.square(np.sin(Alpha))
    
    X1,Y1 = vf.ToCart(A, R1, angle = 'rad')
    
    Xnew,Ynew = Cone2Circle(X1,Y1,Alpha)
    
    Znew = (Z-R)*np.sin(Alpha)

    return(Xnew,Ynew,Znew)

    


def VelCone2Circle(VX,VY,X,Y,Alpha):
    # (VX,VY) velocities to transform, X,Y origin of the velovcity vector, Alpha angle of the cone
    Theta,R = vf.ToCirc(X,Y,angle = 'rad')
    
    central = (Y == 0) & (X<0)
    if np.size(Theta)== 1:
        if central:
            Theta = np.abs(Theta)

    a = (np.cos(Theta)*np.cos(Theta*np.sin(Alpha))/np.sin(Alpha))+ np.sin(Theta)*np.sin(Theta*np.sin(Alpha))
    b = (np.sin(Theta)*np.cos(Theta*np.sin(Alpha))/np.sin(Alpha))- np.cos(Theta)*np.sin(Theta*np.sin(Alpha))
    c = (np.cos(Theta)*np.sin(Theta*np.sin(Alpha))/np.sin(Alpha))- np.sin(Theta)*np.cos(Theta*np.sin(Alpha))
    d = (np.sin(Theta)*np.sin(Theta*np.sin(Alpha))/np.sin(Alpha))+ np.cos(Theta)*np.cos(Theta*np.sin(Alpha))
    
    VXnew = VX*a + VY*b
    VYnew = VX*c + VY*d

    return(VXnew,VYnew)


def VelCircle2Cone(VX,VY,X,Y,Alpha):
    # (VX,VY) velocities to transform, X,Y origin of the velovcity vector, Alpha angle of the cone
    Theta,R = vf.ToCirc(X,Y,angle = 'rad')
    a = (np.cos(Theta)*np.cos(Theta*np.sin(Alpha))/np.sin(Alpha))+ np.sin(Theta)*np.sin(Theta*np.sin(Alpha))
    b = (np.sin(Theta)*np.cos(Theta*np.sin(Alpha))/np.sin(Alpha))- np.cos(Theta)*np.sin(Theta*np.sin(Alpha))
    c = (np.cos(Theta)*np.sin(Theta*np.sin(Alpha))/np.sin(Alpha))- np.sin(Theta)*np.cos(Theta*np.sin(Alpha))
    d = (np.sin(Theta)*np.sin(Theta*np.sin(Alpha))/np.sin(Alpha))+ np.cos(Theta)*np.cos(Theta*np.sin(Alpha))
    
    VXnew = (d*VX - b*VY)/(a*d-b*c)
    VYnew = (c*VX - a*VY)/(b*c-a*d)

    return(VXnew,VYnew)

####
# 3. Compute the height of a column at (X,Y) in a sphere of radius R centered on (Xs,0,0)
def SphereH(R,X,Y,Xs):

    
    H = 2*np.sqrt(R**2-(X-Xs)**2-Y**2)
    
    
    H[~np.isnan(H)] = 0.001*R + H[~np.isnan(H)]
    
    H[np.isnan(H)] = 0
    
    return(H)

###
# 4. Fraction of the drop volume impacting the cone for a impact at distance r from the center
## The computation is done using the formula for the intersection of a sphere and a cylinder 

def volFrac(r,Rd,Rc):
    F = np.empty(len(r))
    for rr,ir in zip(r,range(len(r))):
        if rr+Rd<=Rc:
            F[ir] = 1
        else:                
            F[ir] = vf.interVolSC(Rd,Rc,rr)/(4/3*np.pi*Rd**3)
    return(F*100)



















