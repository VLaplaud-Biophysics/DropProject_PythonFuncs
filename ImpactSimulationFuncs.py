# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 08:04:43 2023

@author: laplaud
"""

import numpy as np

import matplotlib.pyplot as plt

import DropGeometryClasses as dgc

import VallapFunc_DP as vf


###
# 1. Plotting function for different volume fractions
def plotFracs(Angle,npts,DropDiam,ConeDiams):
    
    DropSpeed = 5 # [mm/ms] 

    fig,ax = plt.subplots(dpi =200)
    fig.suptitle('Cone of ' + str(round(Angle*180/np.pi)) + '° angle')
    ax.set_xlabel('Off-centering/Cone radius')
    ax.set_ylabel('Drop fraction in the jet (%)')

    fig1,ax1 = plt.subplots(dpi =200)
    fig1.suptitle('Cone of ' + str(round(Angle*180/np.pi)) + '° angle')
    ax1.set_xlabel('Off-centering/Cone radius')
    ax1.set_ylabel('Impact fraction in the jet (%)')

    fig2,ax2 = plt.subplots(dpi =200)
    fig2.suptitle('Cone of ' + str(round(Angle*180/np.pi)) + '° angle')
    ax2.set_xlabel('Off-centering/Cone radius')
    ax2.set_ylabel('Sheet fraction/Jet fraction')

    cpt = 0
    
    nsim = len(ConeDiams)

    for cd in ConeDiams:
        
        cpt += 1

        print(f'Computing for cone n°{cpt:03d} of {nsim:03d}.', end='\r')
        
        cr = cd/2 
    
        Cone = dgc.Cone(cr,Angle)
        
        OffCmax = 0.495*(cd+DropDiam)

        OffCents = np.linspace(0.1,OffCmax,npts)
        Drops = [dgc.Drop(DropDiam/2,x,50,DropSpeed) for x in OffCents]
        Impacts = [Cone.impact(D) for D in Drops]
        JetFracs = [I.JetFrac for I in Impacts]
        DropFracs = [I.JetFrac*I.VolFrac/100 for I in Impacts]

        lab = 'Drop size / cone size = ' + str(np.round(100*DropDiam/cd)/100)
        
        ax.plot(OffCents/cr,DropFracs,'-*',label=lab)
        ax1.plot(OffCents/cr,JetFracs,'-*',label=lab)
        ax2.plot(OffCents/cr,np.divide(np.subtract(100,JetFracs),JetFracs),'-*',label=lab)

        ax.legend()
        ax1.legend()
        ax2.legend()
        
###
# 2. Phase diagrams

def PhaseDiagrams(Angle,ConeDiam,RelDropDiams,RelOffCents):
    
    meshDD,meshOC = np.meshgrid(RelDropDiams,RelOffCents)
    
    meshDD = meshDD.flatten()
    meshOC = meshOC.flatten()
    
    DropSpeed = 5 # [mm/ms]
    
    JetFracs = np.empty(np.shape(meshDD))
    
    Cone = dgc.Cone(ConeDiam/2,Angle)
    
    nsim = len(meshDD)
    
    cpt = 0
    
    for rdd,roc,i in zip(meshDD,meshOC,range(nsim)):
        
        cpt += 1

        print(f'Computing impacts n°{cpt:1d} of {nsim:1d}.', end='\r')
            
        dd = rdd*ConeDiam
        
        OffCmax = 0.495*(ConeDiam+dd)
        
        oc = roc*ConeDiam/2
            
        if oc < OffCmax:

            Drop = dgc.Drop(dd/2,oc,50,DropSpeed)
            Impact = Cone.impact(Drop)
            JetFracs[i] = Impact.JetFrac
            
        else:
            
            JetFracs[i] = np.nan
            
    
    # Impact fraction in jet
    fig0,ax0 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax0.set_title('Cone angle = ' + str(Angle/(2*np.pi)*360) + '°')
    ax0.set_xlabel('Offcent/ConeRadius')
    ax0.set_ylabel('DropSize/ConeSize')
    
    
    sc0 = ax0.scatter(meshOC,meshDD,c=JetFracs,cmap='jet',s=60)

    cbar0 = plt.colorbar(sc0)
    cbar0.set_label('Impact fraction in the jet (%)')
    fig0.tight_layout()

    # Sheet/Jet volume ratio
    fig1,ax1 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax1.set_title('Cone angle = ' + str(Angle/(2*np.pi)*360) + '°')
    ax1.set_xlabel('Offcent/ConeRadius')
    ax1.set_ylabel('DropSize/ConeSize')
    
    
    sc1 = ax1.scatter(meshOC,meshDD,c=np.divide(np.subtract(100,JetFracs),JetFracs),cmap='RdYlBu', vmin=0, vmax=1.8,s=60)
    
    cbar1 = plt.colorbar(sc1)
    cbar1.set_label('Sheet/Jet volume ratio')
    fig1.tight_layout()

    return(fig0,fig1)