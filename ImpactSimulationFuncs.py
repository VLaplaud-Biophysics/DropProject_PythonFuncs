# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 08:04:43 2023

@author: laplaud
"""

import numpy as np

import matplotlib.pyplot as plt

import DropGeometryClasses as dgc

import VallapFunc_DP as vf

from IPython import get_ipython


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
# 2. Fixed angle diagrams

def FixedAngleDiagrams(Angle,ConeDiam,RelDropDiams,RelOffCents,oriType):
    
    npts = len(RelDropDiams)
    
    meshDD,meshOC = np.meshgrid(RelDropDiams,RelOffCents)
    
    meshDD = meshDD.flatten()
    meshOC = meshOC.flatten()
    
    DropSpeed = 5 # [mm/ms]
    rho = 997e-9 # [kg/mm3]
    
    JetFracs = np.empty(np.shape(meshDD))
    SheetWide = np.empty(np.shape(meshDD))
    JetNRJ = np.empty(np.shape(meshDD))
    
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
            Impact = Cone.impact(Drop,oriType)
            JetFracs[i] = Impact.compute_JetFrac()
            JetNRJ[i] = Impact.VolFrac*Impact.JetFrac/100*4/3*np.pi*(dd/2)**3*rho*DropSpeed**2 # [J]
            # JetFrac[%] * VolFrac[U] * VolDrop[mm3] * WaterDensity[kg/mm3] * DropSpeed²[mm²/ms²] = JetKineticNRJ[kg.mm²/ms² = J]
            SheetWide[i] = Impact.SheetOpening()[0]
            
        else:
            
            JetFracs[i] = np.nan
            SheetWide[i] = np.nan
            JetNRJ[i] = np.nan
       
    pointSize = 250000/npts**2
    
    # Impact fraction in jet
    fig0,ax0 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax0.set_title('Cone angle = ' + str(Angle/(2*np.pi)*360))
    ax0.set_xlabel('Offcent/ConeRadius')
    ax0.set_ylabel('DropSize/ConeSize')
    
    
    sc0 = ax0.scatter(meshOC,meshDD,c=JetFracs,cmap='PuOr',s=pointSize)

    cbar0 = plt.colorbar(sc0)
    cbar0.set_label('Impact fraction in the jet (%)')
    fig0.tight_layout()

    # Sheet/Jet volume ratio
    fig1,ax1 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax1.set_title('Cone angle = ' + str(Angle/(2*np.pi)*360))
    ax1.set_xlabel('Offcent/ConeRadius')
    ax1.set_ylabel('DropSize/ConeSize')
    
    
    sc1 = ax1.scatter(meshOC,meshDD,c=np.divide(np.subtract(100,JetFracs),JetFracs),cmap='plasma',s=pointSize)
    
    cbar1 = plt.colorbar(sc1)
    cbar1.set_label('Sheet/Jet volume ratio')
    fig1.tight_layout()
    
    
    # Sheet/Jet volume ratio
    fig2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax2.set_title('Cone angle = ' + str(Angle/(2*np.pi)*360))
    ax2.set_xlabel('Offcent/ConeRadius')
    ax2.set_ylabel('DropSize/ConeSize')
    
    sc2 = ax2.scatter(meshOC,meshDD,c=SheetWide*360/(2*np.pi),cmap='cividis',s=pointSize)
    
    cbar2 = plt.colorbar(sc2)
    cbar2.set_label('Sheet opening [°]')
    fig2.tight_layout()
    
    # kinetic energy in the jet
    fig3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax3.set_title('Cone angle = ' + str(Angle/(2*np.pi)*360))
    ax3.set_xlabel('Offcent/ConeRadius')
    ax3.set_ylabel('DropSize/ConeSize')
    
    sc3 = ax3.scatter(meshOC,meshDD,c=JetNRJ*1000,cmap='jet',s=pointSize)
    
    cbar3 = plt.colorbar(sc3)
    cbar3.set_label('Maximum kinetic energy in jet [mJ]')
    fig3.tight_layout()

    return(fig0,fig1,fig2,fig3)

###
# 2. Fixed drop diagrams

def FixedDropDiagrams(RelDropDiam,ConeDiam,Angles,RelOffCents,oriType):
    
    npts = len(Angles)
    
    meshA,meshOC = np.meshgrid(Angles,RelOffCents)
    
    meshA = meshA.flatten()
    meshOC = meshOC.flatten()
    
    DropSpeed = 5 # [mm/ms]
    rho = 997e-9 # [kg/mm3]
    
    JetFracs = np.empty(np.shape(meshA))
    SheetWide = np.empty(np.shape(meshA))
    JetNRJ = np.empty(np.shape(meshA))
    JetNRJproba = np.empty(np.shape(meshA))
    

    dd = RelDropDiam*ConeDiam
    
    nsim = len(meshOC)
    
    cpt = 0
    
    for ra,roc,i in zip(meshA,meshOC,range(nsim)):
        
        cpt += 1

        print(f'Computing impacts n°{cpt:1d} of {nsim:1d}.', end='\r')
            
        
        OffCmax = 0.495*(ConeDiam+dd)
        
        oc = roc*ConeDiam/2
            
        if oc < OffCmax:
            
            Cone = dgc.Cone(ConeDiam/2,ra)

            Drop = dgc.Drop(dd/2,oc,50,DropSpeed)
            Impact = Cone.impact(Drop,oriType)
            JetFracs[i] = Impact.compute_JetFrac()
            JetNRJ[i] = Impact.VolFrac*Impact.JetFrac/100*4/3*np.pi*(dd/2)**3*rho*DropSpeed**2 # [J]
            # JetFrac[%] * VolFrac[U] * VolDrop[mm3] * WaterDensity[kg/mm3] * DropSpeed²[mm²/ms²] = JetKineticNRJ[kg.mm²/ms² = J]
            JetNRJproba[i] = Impact.VolFrac*Impact.JetFrac/100*4/3*np.pi*(dd/2)**3*rho*DropSpeed**2*oc/OffCmax # 
            SheetWide[i] = Impact.SheetOpening()[0]
            
        else:
            
            JetFracs[i] = np.nan
            SheetWide[i] = np.nan
            JetNRJ[i] = np.nan
            JetNRJproba[i] = np.nan
       
    pointSize = 250000/npts**2
    
    # Impact fraction in jet
    fig0,ax0 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax0.set_title('DropSize/ConeSize : ' + str(RelDropDiam))
    ax0.set_xlabel('Offcent/ConeRadius')
    ax0.set_ylabel('Cone angle [°]')
    ax0.set_xlim([0,np.max(RelOffCents)])
    
    
    sc0 = ax0.scatter(meshOC,meshA/(2*np.pi)*360,c=JetFracs,cmap='PuOr',s=pointSize)

    cbar0 = plt.colorbar(sc0)
    cbar0.set_label('Impact fraction in the jet (%)')
    fig0.tight_layout()

    # Sheet/Jet volume ratio
    fig1,ax1 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax1.set_title('DropSize/ConeSize : ' + str(RelDropDiam))
    ax1.set_xlabel('Offcent/ConeRadius')
    ax1.set_ylabel('Cone angle [°]')
    ax1.set_xlim([0,np.max(RelOffCents)])
    
    
    sc1 = ax1.scatter(meshOC,meshA/(2*np.pi)*360,c=np.divide(np.subtract(100,JetFracs),JetFracs),cmap='plasma',s=pointSize)
    
    cbar1 = plt.colorbar(sc1)
    cbar1.set_label('Sheet/Jet volume ratio')
    fig1.tight_layout()
    
    
    # Sheet/Jet volume ratio
    fig2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax2.set_title('DropSize/ConeSize : ' + str(RelDropDiam))
    ax2.set_xlabel('Offcent/ConeRadius')
    ax2.set_ylabel('Cone angle [°]')
    ax2.set_xlim([0,np.max(RelOffCents)])
    
    sc2 = ax2.scatter(meshOC,meshA/(2*np.pi)*360,c=SheetWide*360/(2*np.pi),cmap='cividis',s=pointSize)
    
    cbar2 = plt.colorbar(sc2)
    cbar2.set_label('Sheet opening [°]')
    fig2.tight_layout()
    
    # kinetic energy in the jet
    fig3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax3.set_title('DropSize/ConeSize : ' + str(RelDropDiam))
    ax3.set_xlabel('Offcent/ConeRadius')
    ax3.set_ylabel('Cone angle [°]')
    ax3.set_xlim([0,np.max(RelOffCents)])
    
    sc3 = ax3.scatter(meshOC,meshA/(2*np.pi)*360,c=JetNRJ*1000,cmap='jet',s=pointSize)
    
    cbar3 = plt.colorbar(sc3)
    cbar3.set_label('Maximum kinetic energy in jet [mJ]')
    fig3.tight_layout()
    
    # kinetic energy in the jet
    fig4,ax4 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax4.set_title('DropSize/ConeSize : ' + str(RelDropDiam))
    ax4.set_xlabel('Offcent/ConeRadius')
    ax4.set_ylabel('Cone angle [°]')
    ax4.set_xlim([0,np.max(RelOffCents)])
    
    sc4 = ax4.scatter(meshOC,meshA/(2*np.pi)*360,c=JetNRJproba*1000,cmap='jet',s=pointSize)
    
    cbar4 = plt.colorbar(sc4)
    cbar4.set_label('Probable kinetic energy in a jet [mJ]')
    fig4.tight_layout()

    return(fig0,fig1,fig2,fig3,fig4)


###
# 3. Fixed offcent diagrams

def FixedDistDiagrams(RelOffCent,ConeDiam,Angles,RelDropDiams,oriType):
    
    npts = len(Angles)
    
    DD = RelDropDiams*ConeDiam
    
    meshA,meshDD = np.meshgrid(Angles,DD)
    
    meshA = meshA.flatten()
    meshDD = meshDD.flatten()
    
    DropSpeed = 5 # [mm/ms]
    rho = 997e-9 # [kg/mm3]
    
    JetFracs = np.empty(np.shape(meshA))
    SheetWide = np.empty(np.shape(meshA))
    JetNRJ = np.empty(np.shape(meshA))
    JetNRJproba = np.empty(np.shape(meshA))
    
    
    nsim = len(meshDD)
    
    cpt = 0
    
    for a,dd,i in zip(meshA,meshDD,range(nsim)):
        
        cpt += 1

        print(f'Computing impacts n°{cpt:1d} of {nsim:1d}.', end='\r')
            
        
        OffCmax = 0.495*(ConeDiam+dd)
        
        oc = RelOffCent*ConeDiam/2
        
            
        if oc < OffCmax:
            
            Cone = dgc.Cone(ConeDiam/2,a)

            Drop = dgc.Drop(dd/2,oc,50,DropSpeed)
            Impact = Cone.impact(Drop,oriType)
            JetFracs[i] = Impact.compute_JetFrac()
            JetNRJ[i] = Impact.VolFrac*Impact.JetFrac/100*4/3*np.pi*(dd/2)**3*rho*DropSpeed**2 # [J]
            # JetFrac[%] * VolFrac[U] * VolDrop[mm3] * WaterDensity[kg/mm3] * DropSpeed²[mm²/ms²] = JetKineticNRJ[kg.mm²/ms² = J]
            JetNRJproba[i] = Impact.VolFrac*Impact.JetFrac/100*4/3*np.pi*(dd/2)**3*rho*DropSpeed**2*oc/OffCmax # 
            SheetWide[i] = Impact.SheetOpening()[0]
            
        else:
            
            JetFracs[i] = np.nan
            SheetWide[i] = np.nan
            JetNRJ[i] = np.nan
            JetNRJproba[i] = np.nan
       
    pointSize = 250000/npts**2
    
    # Impact fraction in jet
    fig0,ax0 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax0.set_title('Offcent/ConeRadius: ' + str(RelOffCent))
    ax0.set_xlabel('DropSize/ConeSize')
    ax0.set_ylabel('Cone angle [°]')
    ax0.set_xlim([0,np.max(RelDropDiams)])
    
    
    sc0 = ax0.scatter(meshDD,meshA/(2*np.pi)*360,c=JetFracs,cmap='PuOr',s=pointSize)

    cbar0 = plt.colorbar(sc0)
    cbar0.set_label('Impact fraction in the jet (%)')
    fig0.tight_layout()

    # Sheet/Jet volume ratio
    fig1,ax1 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax1.set_title('Offcent/ConeRadius: ' + str(RelOffCent))
    ax1.set_xlabel('DropSize/ConeSize')
    ax1.set_ylabel('Cone angle [°]')
    ax1.set_xlim([0,np.max(RelDropDiams)])
    
    
    sc1 = ax1.scatter(meshDD,meshA/(2*np.pi)*360,c=np.divide(np.subtract(100,JetFracs),JetFracs),cmap='plasma',s=pointSize)
    
    cbar1 = plt.colorbar(sc1)
    cbar1.set_label('Sheet/Jet volume ratio')
    fig1.tight_layout()
    
    
    # Sheet/Jet volume ratio
    fig2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax2.set_title('Offcent/ConeRadius: ' + str(RelOffCent))
    ax2.set_xlabel('DropSize/ConeSize')
    ax2.set_ylabel('Cone angle [°]')
    ax2.set_xlim([0,np.max(RelDropDiams)])
    
    sc2 = ax2.scatter(meshDD,meshA/(2*np.pi)*360,c=SheetWide*360/(2*np.pi),cmap='cividis',s=pointSize)
    
    cbar2 = plt.colorbar(sc2)
    cbar2.set_label('Sheet opening [°]')
    fig2.tight_layout()
    
    # kinetic energy in the jet
    fig3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax3.set_title('Offcent/ConeRadius: ' + str(RelOffCent))
    ax3.set_xlabel('DropSize/ConeSize')
    ax3.set_ylabel('Cone angle [°]')
    ax3.set_xlim([0,np.max(RelDropDiams)])
    
    sc3 = ax3.scatter(meshDD,meshA/(2*np.pi)*360,c=JetNRJ*1000,cmap='jet',s=pointSize)
    
    cbar3 = plt.colorbar(sc3)
    cbar3.set_label('Maximum kinetic energy in jet [mJ]')
    fig3.tight_layout()
    
    # kinetic energy in the jet
    fig4,ax4 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax4.set_title('Offcent/ConeRadius: ' + str(RelOffCent))
    ax4.set_xlabel('DropSize/ConeSize')
    ax4.set_ylabel('Cone angle [°]')
    ax4.set_xlim([0,np.max(RelDropDiams)])
    
    sc4 = ax4.scatter(meshDD,meshA/(2*np.pi)*360,c=JetNRJproba*1000,cmap='jet',s=pointSize)
    
    cbar4 = plt.colorbar(sc4)
    cbar4.set_label('Probable kinetic energy in a jet [mJ]')
    fig4.tight_layout()

    return(fig0,fig1,fig2,fig3,fig4)




###
# 3. Fixed offcent diagrams

def PhaseDiagrams(RelOffCents,ConeDiam,Angles,RelDropDiams,oriType):
    
    npts = len(Angles)
    
    DD = RelDropDiams*ConeDiam
    
    OC = RelOffCents*ConeDiam/2
    
    meshA,meshDD,meshOC = np.meshgrid(Angles,DD,OC)
    
    DropSpeed = 5 # [mm/ms]
    rho = 997e-9 # [kg/mm3]
    
    JetFracs = np.empty(np.shape(meshA))
    SheetWide = np.empty(np.shape(meshA))
    JetNRJ = np.empty(np.shape(meshA))
    JetNRJproba = np.empty(np.shape(meshA))
    
    
    nsim = meshDD.size
    
    cpt = 0
    
    for a,dd,oc in zip(meshA.flatten(),meshDD.flatten(),meshOC.flatten()):
        
        idx = np.argwhere((meshA==a)&(meshDD==dd)&(meshOC==oc))
        
        cpt += 1

        print(f'Computing impacts n°{cpt:1d} of {nsim:1d}.', end='\r')
            
        
        OffCmax = 0.495*(ConeDiam+dd)
        
            
        if oc < OffCmax:
            
            Cone = dgc.Cone(ConeDiam/2,a)

            Drop = dgc.Drop(dd/2,oc,50,DropSpeed)
            Impact = Cone.impact(Drop,oriType)
            JetFracs[idx[0][0],idx[0][1],idx[0][2]] = Impact.compute_JetFrac()
            JetNRJ[idx[0][0],idx[0][1],idx[0][2]] = Impact.VolFrac*Impact.JetFrac/100*4/3*np.pi*(dd/2)**3*rho*DropSpeed**2 # [J]
            # JetFrac[%] * VolFrac[U] * VolDrop[mm3] * WaterDensity[kg/mm3] * DropSpeed²[mm²/ms²] = JetKineticNRJ[kg.mm²/ms² = J]
            JetNRJproba[idx[0][0],idx[0][1],idx[0][2]] = Impact.VolFrac*Impact.JetFrac/100*4/3*np.pi*(dd/2)**3*rho*DropSpeed**2*oc/OffCmax # 
            SheetWide[idx[0][0],idx[0][1],idx[0][2]] = Impact.SheetOpening()[0]
            
        else:
            
            JetFracs[idx[0][0],idx[0][1],idx[0][2]] = np.nan
            SheetWide[idx[0][0],idx[0][1],idx[0][2]] = np.nan
            JetNRJ[idx[0][0],idx[0][1],idx[0][2]] = np.nan
            JetNRJproba[idx[0][0],idx[0][1],idx[0][2]] = np.nan
       
    pointSize = 250000/npts**2
    pointSize3d = 2500/npts**2
    
    get_ipython().run_line_magic('matplotlib', 'qt')
    fig = plt.figure(figsize=(5,5),dpi = 200)
    ax = fig.add_subplot(projection='3d')

    sc = ax.scatter(meshA,meshDD,meshOC,c=JetFracs,cmap='PuOr',s=pointSize3d, alpha = 0.5)
    
    ax.set_title('Impact fraction in the jet (%)')
    ax.set_xlabel('Cone angle [°]')
    ax.set_ylabel('DropSize/ConeSize')
    ax.set_zlabel('Offcent/ConeRadius')
    
    cbar = plt.colorbar(sc)
    cbar.set_label('Impact fraction in the jet (%)')
    
    fig.tight_layout()
    
    # # Impact fraction in jet
    # fig0,ax0 = plt.subplots(dpi=150,figsize = (7,6)) 
    # ax0.set_title('Offcent/ConeRadius: ' + str(RelOffCent))
    # ax0.set_xlabel('DropSize/ConeSize')
    # ax0.set_ylabel('Cone angle [°]')
    # ax0.set_xlim([0,np.max(RelDropDiams)])
    
    
    # sc0 = ax0.scatter(meshDD,meshA/(2*np.pi)*360,c=JetFracs,cmap='PuOr',s=pointSize)

    # cbar0 = plt.colorbar(sc0)
    # cbar0.set_label('Impact fraction in the jet (%)')
    # fig0.tight_layout()

    # # Sheet/Jet volume ratio
    # fig1,ax1 = plt.subplots(dpi=150,figsize = (7,6)) 
    # ax1.set_title('Offcent/ConeRadius: ' + str(RelOffCent))
    # ax1.set_xlabel('DropSize/ConeSize')
    # ax1.set_ylabel('Cone angle [°]')
    # ax1.set_xlim([0,np.max(RelDropDiams)])
    
    
    # sc1 = ax1.scatter(meshDD,meshA/(2*np.pi)*360,c=np.divide(np.subtract(100,JetFracs),JetFracs),cmap='plasma',s=pointSize)
    
    # cbar1 = plt.colorbar(sc1)
    # cbar1.set_label('Sheet/Jet volume ratio')
    # fig1.tight_layout()
    
    
    # # Sheet/Jet volume ratio
    # fig2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
    # ax2.set_title('Offcent/ConeRadius: ' + str(RelOffCent))
    # ax2.set_xlabel('DropSize/ConeSize')
    # ax2.set_ylabel('Cone angle [°]')
    # ax2.set_xlim([0,np.max(RelDropDiams)])
    
    # sc2 = ax2.scatter(meshDD,meshA/(2*np.pi)*360,c=SheetWide*360/(2*np.pi),cmap='cividis',s=pointSize)
    
    # cbar2 = plt.colorbar(sc2)
    # cbar2.set_label('Sheet opening [°]')
    # fig2.tight_layout()
    
    # # kinetic energy in the jet
    # fig3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
    # ax3.set_title('Offcent/ConeRadius: ' + str(RelOffCent))
    # ax3.set_xlabel('DropSize/ConeSize')
    # ax3.set_ylabel('Cone angle [°]')
    # ax3.set_xlim([0,np.max(RelDropDiams)])
    
    # sc3 = ax3.scatter(meshDD,meshA/(2*np.pi)*360,c=JetNRJ*1000,cmap='jet',s=pointSize)
    
    # cbar3 = plt.colorbar(sc3)
    # cbar3.set_label('Maximum kinetic energy in jet [mJ]')
    # fig3.tight_layout()
    
    # # kinetic energy in the jet
    # fig4,ax4 = plt.subplots(dpi=150,figsize = (7,6)) 
    # ax4.set_title('Offcent/ConeRadius: ' + str(RelOffCent))
    # ax4.set_xlabel('DropSize/ConeSize')
    # ax4.set_ylabel('Cone angle [°]')
    # ax4.set_xlim([0,np.max(RelDropDiams)])
    
    # sc4 = ax4.scatter(meshDD,meshA/(2*np.pi)*360,c=JetNRJproba*1000,cmap='jet',s=pointSize)
    
    # cbar4 = plt.colorbar(sc4)
    # cbar4.set_label('Probable kinetic energy in a jet [mJ]')
    # fig4.tight_layout()

    # return(fig0,fig1,fig2,fig3,fig4)