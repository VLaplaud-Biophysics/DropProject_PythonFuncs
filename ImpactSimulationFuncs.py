# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 08:04:43 2023

@author: laplaud
"""

import numpy as np

import matplotlib.pyplot as plt

import DropGeometryClasses as dgc

import VallapFunc_DP as vf

import os

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
# 2. Parameter space diagrams

def PhaseDiagrams(RelOffCents,ConeDiam,Angles,RelDropDiams,oriType,DiagDim,label):
    
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
      
        
    if DiagDim == 2:
        
        print('\n\nMaking and saving figures',end='\r')
        
        pointSize = 250000/npts**2 
        
        
        savepath = r'd:\Users\laplaud\Desktop\PostDoc\Code\DropProject_WithAna\Figures/' + label
        
                
        if not os.path.exists(savepath):
            os.mkdir(savepath) # create folder
        
        if not os.path.exists(savepath + '\JetFrac'):
            os.mkdir(savepath + '\JetFrac') # create folder
            
            
        if not os.path.exists(savepath + '\VolRatio'):
            os.mkdir(savepath + '\VolRatio') # create folder
            
            
        if not os.path.exists(savepath + '\SheetOpening'):
            os.mkdir(savepath + '\SheetOpening') # create folder
            
            
        if not os.path.exists(savepath + '\JetNRJ'):
            os.mkdir(savepath + '\JetNRJ') # create folder
        
        ########### Fixed angle diagrams ###########       

        
        for a,ia in zip(Angles,range(len(Angles))):
            
            # Impact fraction in jet
            f0,ax0 = plt.subplots(dpi=150,figsize = (7,6)) 
            ax0.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
            ax0.set_xlabel('Offcent/ConeRadius')
            ax0.set_ylabel('DropSize/ConeSize')
            
            
            sc0 = ax0.scatter(meshOC[:,ia,:]*2/ConeDiam,meshDD[:,ia,:]/ConeDiam,c=JetFracs[:,ia,:],cmap='PuOr',s=pointSize)

            cbar0 = plt.colorbar(sc0)
            cbar0.set_label('Impact fraction in the jet (%)')
            f0.tight_layout()
            
            f0.savefig(savepath + '\JetFrac\FxdADgm_'
                       + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_JetFrac.png')

            plt.close(f0)
            
            
            print('Making and saving figures     ',end='\r')

            # Sheet/Jet volume ratio
            f1,ax1 = plt.subplots(dpi=150,figsize = (7,6)) 
            ax1.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
            ax1.set_xlabel('Offcent/ConeRadius')
            ax1.set_ylabel('DropSize/ConeSize')
            
            
            sc1 = ax1.scatter(meshOC[:,ia,:]*2/ConeDiam,meshDD[:,ia,:]/ConeDiam,c=np.divide(np.subtract(100,JetFracs[:,ia,:]),JetFracs[:,ia,:]),cmap='plasma',s=pointSize)
            
            cbar1 = plt.colorbar(sc1)
            cbar1.set_label('Sheet/Jet volume ratio')
            f1.tight_layout()
            
            f1.savefig(savepath + '\VolRatio\FxdADgm_'
                       + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_VolRatio.png')

            plt.close(f1)
            
            
            print('Making and saving figures.     ',end='\r')
            
            
            # Sheet/Jet volume ratio
            f2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
            ax2.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
            ax2.set_xlabel('Offcent/ConeRadius')
            ax2.set_ylabel('DropSize/ConeSize')
            
            sc2 = ax2.scatter(meshOC[:,ia,:]*2/ConeDiam,meshDD[:,ia,:]/ConeDiam,c=SheetWide[:,ia,:]*360/(2*np.pi),cmap='cividis',s=pointSize)
            
            cbar2 = plt.colorbar(sc2)
            cbar2.set_label('Sheet opening [°]')
            f2.tight_layout()
            
            f2.savefig(savepath + '\SheetOpening\FxdADgm_'
                       + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_SheetOpen.png')

            plt.close(f2)
            
            
            print('Making and saving figures..     ',end='\r')
            
            # kinetic energy in the jet
            f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
            ax3.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
            ax3.set_xlabel('Offcent/ConeRadius')
            ax3.set_ylabel('DropSize/ConeSize')
            
            sc3 = ax3.scatter(meshOC[:,ia,:]*2/ConeDiam,meshDD[:,ia,:]/ConeDiam,c=JetNRJ[:,ia,:]*1000,cmap='jet',s=pointSize)
            
            cbar3 = plt.colorbar(sc3)
            cbar3.set_label('Maximum kinetic energy in jet [mJ]')
            f3.tight_layout()
            
            f3.savefig(savepath + '\JetNRJ\FxdADgm_'
                       + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_JetNRJ.png')

            plt.close(f3)
            
            
            print('Making and saving figures...     ',end='\r')

        
        ########### Fixed drop diagrams ###########
        
        for rdd,idd in zip(RelDropDiams,range(len(RelDropDiams))):
            
            # Impact fraction in jet
            f0,ax0 = plt.subplots(dpi=150,figsize = (7,6)) 
            ax0.set_title('DropSize/ConeSize : ' + str(round(rdd*10)/10))
            ax0.set_xlabel('Offcent/ConeRadius')
            ax0.set_ylabel('Cone angle [°]')
            ax0.set_xlim([0,np.max(RelOffCents)])
            
            
            sc0 = ax0.scatter(meshOC[idd,:,:]*2/ConeDiam,meshA[idd,:,:]/(2*np.pi)*360,c=JetFracs[idd,:,:],cmap='PuOr',s=pointSize)

            cbar0 = plt.colorbar(sc0)
            cbar0.set_label('Impact fraction in the jet (%)')
            f0.tight_layout()
            
            f0.savefig(savepath + '\JetFrac\FxdDrDgm_'
               + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SizeRatio_JetFrac.png')
    
            plt.close(f0)
            
            
            print('Making and saving figures     ',end='\r')
    

            # Sheet/Jet volume ratio
            f1,ax1 = plt.subplots(dpi=150,figsize = (7,6)) 
            ax1.set_title('DropSize/ConeSize : ' + str(round(rdd*10)/10))
            ax1.set_xlabel('Offcent/ConeRadius')
            ax1.set_ylabel('Cone angle [°]')
            ax1.set_xlim([0,np.max(RelOffCents)])
            
            
            sc1 = ax1.scatter(meshOC[idd,:,:]*2/ConeDiam,meshA[idd,:,:]/(2*np.pi)*360,c=np.divide(np.subtract(100,JetFracs[idd,:,:]),JetFracs[idd,:,:]),cmap='plasma',s=pointSize)
            
            cbar1 = plt.colorbar(sc1)
            cbar1.set_label('Sheet/Jet volume ratio')
            f1.tight_layout()
            
            f1.savefig(savepath + '\VolRatio\FxdDrDgm_'
                       + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SizeRatio_VolRatio.png')
            
            plt.close(f1)
            
            
            print('Making and saving figures.     ',end='\r')
            

            # Sheet/Jet volume ratio
            f2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
            ax2.set_title('DropSize/ConeSize : ' + str(round(rdd*10)/10))
            ax2.set_xlabel('Offcent/ConeRadius')
            ax2.set_ylabel('Cone angle [°]')
            ax2.set_xlim([0,np.max(RelOffCents)])
            
            sc2 = ax2.scatter(meshOC[idd,:,:]*2/ConeDiam,meshA[idd,:,:]/(2*np.pi)*360,c=SheetWide[idd,:,:]*360/(2*np.pi),cmap='cividis',s=pointSize)
            
            cbar2 = plt.colorbar(sc2)
            cbar2.set_label('Sheet opening [°]')
            f2.tight_layout()
            
            f2.savefig(savepath + '\SheetOpening\FxdDrDgm_'
                       + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SizeRatio_SheetOpen.png')
            
            plt.close(f2)
            
            
            print('Making and saving figures..     ',end='\r')
            
            
            # kinetic energy in the jet
            f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
            ax3.set_title('DropSize/ConeSize : ' + str(round(rdd*10)/10))
            ax3.set_xlabel('Offcent/ConeRadius')
            ax3.set_ylabel('Cone angle [°]')
            ax3.set_xlim([0,np.max(RelOffCents)])
            
            sc3 = ax3.scatter(meshOC[idd,:,:]*2/ConeDiam,meshA[idd,:,:]/(2*np.pi)*360,c=JetNRJ[idd,:,:]*1000,cmap='jet',s=pointSize)
            
            cbar3 = plt.colorbar(sc3)
            cbar3.set_label('Maximum kinetic energy in jet [mJ]')
            f3.tight_layout()
            
            f3.savefig(savepath + '\JetNRJ\FxdDrDgm_'
                       + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SizeRatio_JetNRJ.png')
            
            plt.close(f3)
            
            
            print('Making and saving figures...     ',end='\r')
            
            # kinetic energy in the jet
            f4,ax4 = plt.subplots(dpi=150,figsize = (7,6)) 
            ax4.set_title('DropSize/ConeSize : ' + str(round(rdd*10)/10))
            ax4.set_xlabel('Offcent/ConeRadius')
            ax4.set_ylabel('Cone angle [°]')
            ax4.set_xlim([0,np.max(RelOffCents)])
            
            sc4 = ax4.scatter(meshOC[idd,:,:]*2/ConeDiam,meshA[idd,:,:]/(2*np.pi)*360,c=JetNRJproba[idd,:,:]*1000,cmap='jet',s=pointSize)
            
            cbar4 = plt.colorbar(sc4)
            cbar4.set_label('Probable kinetic energy in a jet [mJ]')
            f4.tight_layout()
            
            f4.savefig(savepath + '\JetNRJ\FxdDrDgm_'
                       + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SizeRatio_JetNRJproba.png')
            
            plt.close(f4)
            
            
            print('Making and saving figures....     ',end='\r')
            
            
            
        ########### Fixed Offcent diagrams ###########
        
        for roc,ioc in zip(RelOffCents,range(len(RelOffCents))):
        
            # Impact fraction in jet
            f0,ax0 = plt.subplots(dpi=150,figsize = (7,6)) 
            ax0.set_title('Offcent/ConeRadius: ' + str(round(roc*10)/10))
            ax0.set_xlabel('DropSize/ConeSize')
            ax0.set_ylabel('Cone angle [°]')
            ax0.set_xlim([0,np.max(RelDropDiams)])
            
            
            sc0 = ax0.scatter(meshDD[:,:,ioc]/ConeDiam,meshA[:,:,ioc]/(2*np.pi)*360,c=JetFracs[:,:,ioc],cmap='PuOr',s=pointSize)
        
            cbar0 = plt.colorbar(sc0)
            cbar0.set_label('Impact fraction in the jet (%)')
            f0.tight_layout()
            
            f0.savefig(savepath + '\JetFrac\FxdDiDgm_'
             + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OffCent_JetFrac.png')
  
            plt.close(f0)
            
            
            print('Making and saving figures     ',end='\r')
    
            # Sheet/Jet volume ratio
            f1,ax1 = plt.subplots(dpi=150,figsize = (7,6)) 
            ax1.set_title('Offcent/ConeRadius: ' + str(round(roc*10)/10))
            ax1.set_xlabel('DropSize/ConeSize')
            ax1.set_ylabel('Cone angle [°]')
            ax1.set_xlim([0,np.max(RelDropDiams)])
            
            
            sc1 = ax1.scatter(meshDD[:,:,ioc]/ConeDiam,meshA[:,:,ioc]/(2*np.pi)*360,c=np.divide(np.subtract(100,JetFracs[:,:,ioc]),JetFracs[:,:,ioc]),cmap='plasma',s=pointSize)
            
            cbar1 = plt.colorbar(sc1)
            cbar1.set_label('Sheet/Jet volume ratio')
            f1.tight_layout()
            
            f1.savefig(savepath + '\VolRatio\FxdDiDgm_'
                       + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OffCent_VolRatio.png')
            
            plt.close(f1)
            
            
            print('Making and saving figures.     ',end='\r')
            
            # Sheet/Jet volume ratio
            f2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
            ax2.set_title('Offcent/ConeRadius: ' + str(round(roc*10)/10))
            ax2.set_xlabel('DropSize/ConeSize')
            ax2.set_ylabel('Cone angle [°]')
            ax2.set_xlim([0,np.max(RelDropDiams)])
            
            sc2 = ax2.scatter(meshDD[:,:,ioc]/ConeDiam,meshA[:,:,ioc]/(2*np.pi)*360,c=SheetWide[:,:,ioc]*360/(2*np.pi),cmap='cividis',s=pointSize)
            
            cbar2 = plt.colorbar(sc2)
            cbar2.set_label('Sheet opening [°]')
            f2.tight_layout()
            
            f2.savefig(savepath + '\SheetOpening\FxdDiDgm_'
                       + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OffCent_SheetOpen.png')
            
            plt.close(f2)
            
            
            print('Making and saving figures..     ',end='\r')
            
            # kinetic energy in the jet
            f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
            ax3.set_title('Offcent/ConeRadius: ' + str(round(roc*10)/10))
            ax3.set_xlabel('DropSize/ConeSize')
            ax3.set_ylabel('Cone angle [°]')
            ax3.set_xlim([0,np.max(RelDropDiams)])
            
            sc3 = ax3.scatter(meshDD[:,:,ioc]/ConeDiam,meshA[:,:,ioc]/(2*np.pi)*360,c=JetNRJ[:,:,ioc]*1000,cmap='jet',s=pointSize)
            
            cbar3 = plt.colorbar(sc3)
            cbar3.set_label('Maximum kinetic energy in jet [mJ]')
            f3.tight_layout()
            
            f3.savefig(savepath + '\JetNRJ\FxdDiDgm_'
                       + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OffCent_JetNRJ.png')
            
            plt.close(f3)
            
            
            print('Making and saving figures...     ',end='\r')
            
            # kinetic energy in the jet
            f4,ax4 = plt.subplots(dpi=150,figsize = (7,6)) 
            ax4.set_title('Offcent/ConeRadius: ' + str(round(roc*10)/10))
            ax4.set_xlabel('DropSize/ConeSize')
            ax4.set_ylabel('Cone angle [°]')
            ax4.set_xlim([0,np.max(RelDropDiams)])
            
            sc4 = ax4.scatter(meshDD[:,:,ioc]/ConeDiam,meshA[:,:,ioc]/(2*np.pi)*360,c=JetNRJproba[:,:,ioc]*1000,cmap='jet',s=pointSize)
            
            cbar4 = plt.colorbar(sc4)
            cbar4.set_label('Probable kinetic energy in a jet [mJ]')
            f4.tight_layout()
            
            f4.savefig(savepath + '\JetNRJ\FxdDiDgm_'
                       + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OffCent_JetNRJproba.png')
            
            plt.close(f4)
            
            
            print('Making and saving figures....     ',end='\r')
            
        print('Making and saving figures :       \nAll figures done and saved !!')
              
    
    if DiagDim == 3:
        pointSize3d = 2500/npts**2
        
        get_ipython().run_line_magic('matplotlib', 'qt')
        
        
        fig0 = plt.figure(figsize=(5,5),dpi = 200)
        ax0 = fig0.add_subplot(projection='3d')
    
        sc0 = ax0.scatter(meshA,meshDD,meshOC,c=JetFracs,cmap='PuOr',s=pointSize3d, alpha = 0.5)
        
        ax0.set_title('Impact fraction in the jet (%)')
        ax0.set_xlabel('Cone angle [°]')
        ax0.set_ylabel('DropSize/ConeSize')
        ax0.set_zlabel('Offcent/ConeRadius')
        
        cbar0 = plt.colorbar(sc0)
        cbar0.set_label('Impact fraction in the jet (%)')
        
        fig0.tight_layout()
        
        
        fig1 = plt.figure(figsize=(5,5),dpi = 200)
        ax1 = fig1.add_subplot(projection='3d')
    
        sc1 = ax1.scatter(meshA,meshDD,meshOC,c=JetNRJ,cmap='Jet',s=pointSize3d, alpha = 0.5)
        
        ax1.set_title('Maximum kinetic energy in jet [mJ]')
        ax1.set_xlabel('Cone angle [°]')
        ax1.set_ylabel('DropSize/ConeSize')
        ax1.set_zlabel('Offcent/ConeRadius')
        
        cbar1 = plt.colorbar(sc1)
        cbar1.set_label('Impact fraction in the jet (%)')
        
        fig1.tight_layout()
        
        
        # fig = plt.figure(figsize=(5,5),dpi = 200)
        # ax = fig.add_subplot(projection='3d')
    
        # sc = ax.scatter(meshA,meshDD,meshOC,c=NONE,cmap='jet',s=pointSize3d, alpha = 0.5)
        
        # ax.set_title('TBD')
        # ax.set_xlabel('Cone angle [°]')
        # ax.set_ylabel('DropSize/ConeSize')
        # ax.set_zlabel('Offcent/ConeRadius')
        
        # cbar = plt.colorbar(sc)
        # cbar.set_label('Impact fraction in the jet (%)')
        
        # fig.tight_layout()
        
        return(fig0,fig1)

