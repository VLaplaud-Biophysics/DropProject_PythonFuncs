# -*- coding: utf-8 -*-
"""
Created on Mon May 13 09:14:37 2024

@author: laplaud
"""
import numpy as np 

import ImpactSimulationFuncs as isf


coneSurface = 7.3 # [mm²] polymorpha :7.3 globosa : 12.1

npts = 11 # resolution for the diagram

coneAngles = np.linspace(0.1,89.9,npts)/360*2*np.pi

ndrops = 300
dropScaling = np.linspace(0.1,10, npts)

label = 'SuperlargeScale'

isf.OptiDiagrams(coneSurface,coneAngles,npts,ndrops,dropScaling,label)  