# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 09:11:40 2023

@author: laplaud

"""

import numpy as np

import VallapFunc_DP as vf

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

def Cone2Circle(X,Y,Alpha,Ad):
    # (X,Y) points to transform, Alpha angle of the cone, Ad angle of removed sector bissecant (cone config)
    A,R = vf.ToCirc(X,Y,angle = 'rad')
    Xnew,Ynew = vf.ToCart((A-Ad)*np.sin(Alpha)+Ad,R/np.sin(Alpha),angle = 'rad')

    return(Xnew,Ynew)


def Circle2Cone(X,Y,Alpha,Ad): 
    # (X,Y) points to transform, Alpha angle of the cone, Ad angle of removed sector bissecant (circle config)
    A,R = vf.ToCirc(X,Y,angle = 'rad')
    Xnew,Ynew = vf.ToCart((A-Ad)/np.sin(Alpha)+Ad,R*np.sin(Alpha),angle = 'rad')
    
    return(Xnew,Ynew)

####
# 3. Compute the height of a column at (X,Y) in a sphere of radius R centered on (Xs,0,0)
def SphereH(R,X,Y,Xs):

    
    H = 2*np.sqrt(R**2-(X-Xs)**2-Y**2)
    
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



















