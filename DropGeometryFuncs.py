# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 09:11:40 2023

@author: laplaud

"""

import numpy as np
import numpy.pi, numpy.cos, numpy.sin as pi,cos,sin

# 1. Compute cone geometrical transformation parameters

""" For a specific cone, of diameter Rc and angle Alpha with the vertical, we compute :
the radius Ro of the original circle, 
the angle Beta of the sector removed from the circle, 
the fraction Fr of circle surface in the cone"""

def coneGeometry(Rcone,Alpha):
    # perimeter of cone 
    Pcone = 2*pi*Rcone # (in mm)

    # Cone Side length
    Lcone = Rcone/sin(Alpha)

    """ Unfolded cone as a cropped circle of radius Lcone """

    Rcircle = Lcone.copy()

    # Perimter of full circle
    Pcircle = 2*pi*Rcircle

    # Fraction of the circle remaining to form the cone
    Fr = Pcone/Pcircle # analytically Fr = sin(Alpha)

    # Angle of the removed circle part
    Beta = 2*pi*(1-Fr) # analytically Beta = 2*pi*(1-sin(Alpha))
        
    return(Rcircle,Beta,Fr)