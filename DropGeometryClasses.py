# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 09:11:18 2023

@author: laplaud
"""

import DropGeometryFuncs as dgf


# 1. A class for the cone where the drop is going to impact.
class Cone:

    def __init__(self, radius, angle):

        self.Rcone = radius

        self.Alpha = angle
        
        self.Rcircle,self.Beta,self.Frac = dgf.coneGeometry(radius,angle)
        