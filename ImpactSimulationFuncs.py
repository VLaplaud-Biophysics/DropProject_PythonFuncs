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

import shutil

from IPython import get_ipython


###
# 1. Plotting function for different volume fractions
def plotFracs(Angle,npts,DropDiam,ConeDiams,oriType,velType,velIni,meshType):
    
    DropSpeed = 5 # [mm/ms] 


    fig1,ax1 = plt.subplots(dpi =200)
    fig1.suptitle('Cone of ' + str(round(Angle*180/np.pi)) + '° angle')
    ax1.set_xlabel('Off-centering/Cone radius')
    ax1.set_ylabel('Impact fraction in the jet (%)')
    
    
    fig2,ax2 = plt.subplots(dpi =200)
    fig2.suptitle('Cone of ' + str(round(Angle*180/np.pi)) + '° angle')
    ax2.set_xlabel('Off-centering/Cone radius')
    ax2.set_ylabel('Sheet fraction/Jet fraction')

    
    fig3,ax3 = plt.subplots(dpi =200)
    fig3.suptitle('Cone of ' + str(round(Angle*180/np.pi)) + '° angle')
    ax3.set_xlabel('Off-centering/Cone radius')
    ax3.set_ylabel('Drop fraction in the jet (%)')


    cpt = 0
    
    nsim = len(ConeDiams)

    for cd in ConeDiams:
        
        cpt += 1

        print(f'Computing for cone n°{cpt:03d} of {nsim:03d}.', end='\r')
        
        cr = cd/2 
    
        Cone = dgc.Cone(cr,Angle)
        
        OffCmax = 0.45*(cd+DropDiam)

        OffCents = np.linspace(0.1,OffCmax,npts)
        Drops = [dgc.Drop(DropDiam/2,x,71,DropSpeed) for x in OffCents]
        Impacts = [Cone.impact(D,oriType,velIni,meshType) for D in Drops]
        JetFracs = [I.compute_JetFrac(velType) for I in Impacts]
        DropFracs = [I.compute_JetFrac(velType)*I.VolFrac/100 for I in Impacts]

        lab = 'Drop size / cone size = ' + str(np.round(100*DropDiam/cd)/100)
        
        ax3.plot(OffCents/cr,DropFracs,'-*',label=lab)
        ax1.plot(OffCents/cr,JetFracs,'-*',label=lab)
        ax2.plot(OffCents/cr,np.divide(np.subtract(100,JetFracs),JetFracs),'-*',label=lab)

        ax3.legend()
        ax1.legend()
        ax2.legend()
     
###
# 2. Parameter space diagrams

def PhaseDiagrams(RelOffCents,ConeSize,ConeSizeType,Angles,RelDropDiams,oriType,velType,velIni,meshType,label):
    
    npts = len(Angles)
    
    if ConeSizeType == 'surface':        
        
        Normalisation = np.sqrt(ConeSize/np.pi)
        
        NormStr = '_Norm'
        
    elif ConeSizeType == 'radius':   
        
        Normalisation = ConeSize
        NormStr = '/ConeSize'
        
    
    
    
    DD = RelDropDiams*Normalisation
    OC = RelOffCents*Normalisation/2
        
    meshA,meshDD,meshOC = np.meshgrid(Angles,DD,OC)    
    
    DropSpeed = 5 # [mm/ms]
    
    JetFracs = np.empty(np.shape(meshA))
    JetMass = np.empty(np.shape(meshA))
    JetNRJ = np.empty(np.shape(meshA))
    JetNRJ_Bal = np.empty(np.shape(meshA))
    SheetWide = np.empty(np.shape(meshA))
    DispertionDist = np.empty(np.shape(meshA))
    DispertionDist_Var = np.empty(np.shape(meshA))
    
    
    nsim = meshDD.size
    
    savepath = r'd:\Users\laplaud\Desktop\PostDoc\Code\DropProject_WithAna\Figures/' + label + '_' + str(npts) + 'npts'
    
    
    if os.path.exists(savepath+ '\JetFrac'):
        
        print('Clearing previous figure folder...', end = '')
        
        shutil.rmtree(savepath + '\JetFrac\\')   
        shutil.rmtree(savepath + '\VolRatio\\')   
        shutil.rmtree(savepath + '\EnRJ_Bal\\')   
        shutil.rmtree(savepath + '\EnRJ\\')   
        shutil.rmtree(savepath + '\DispertionDist\\')   
        shutil.rmtree(savepath + '\DispertionDist_Var\\')   
        shutil.rmtree(savepath + '\SheetOpening\\') 
        
        print('Done.\n')
        

    print('Creating figures folders...', end = '')
    
    os.makedirs(savepath + '\JetFrac\FixedAngle',exist_ok=True) # create folder
    os.makedirs(savepath + '\JetFrac\FixedDrop',exist_ok=True) # create folder
    os.makedirs(savepath + '\JetFrac\FixedDist',exist_ok=True) # create folder

    os.makedirs(savepath + '\VolRatio\FixedAngle',exist_ok=True) # create folder
    os.makedirs(savepath + '\VolRatio\FixedDrop',exist_ok=True) # create folder
    os.makedirs(savepath + '\VolRatio\FixedDist',exist_ok=True) # create folder

    os.makedirs(savepath + '\EnRJ_Bal\FixedAngle',exist_ok=True) # create folder
    os.makedirs(savepath + '\EnRJ_Bal\FixedDrop',exist_ok=True) # create folder
    os.makedirs(savepath + '\EnRJ_Bal\FixedDist',exist_ok=True) # create folder

    os.makedirs(savepath + '\EnRJ\FixedAngle',exist_ok=True) # create folder
    os.makedirs(savepath + '\EnRJ\FixedDrop',exist_ok=True) # create folder
    os.makedirs(savepath + '\EnRJ\FixedDist',exist_ok=True) # create folder
     
    os.makedirs(savepath + '\DispertionDist\FixedAngle',exist_ok=True) # create folder
    os.makedirs(savepath + '\DispertionDist\FixedDrop',exist_ok=True) # create folder
    os.makedirs(savepath + '\DispertionDist\FixedDist',exist_ok=True) # create folder
     
    os.makedirs(savepath + '\DispertionDist_Var\FixedAngle',exist_ok=True) # create folder
    os.makedirs(savepath + '\DispertionDist_Var\FixedDrop',exist_ok=True) # create folder
    os.makedirs(savepath + '\DispertionDist_Var\FixedDist',exist_ok=True) # create folder
     
    os.makedirs(savepath + '\SheetOpening\FixedAngle',exist_ok=True) # create folder
    os.makedirs(savepath + '\SheetOpening\FixedDrop',exist_ok=True) # create folder
    os.makedirs(savepath + '\SheetOpening\FixedDist',exist_ok=True) # create folder
        
            
    print('Done.\n')
    
    if os.path.exists(savepath) & os.path.exists(savepath + '\Data_DispertionDist_Var_' + str(npts) + 'npts.npy'):
        
        print('Loading previous simulations results...', end = '')
        
        JetFracs = np.load(savepath + '\Data_JetFracs_' + str(npts) + 'npts.npy')
        JetMass = np.load(savepath + '\Data_JetMass_' + str(npts) + 'npts.npy')
        JetNRJ = np.load(savepath + '\Data_JetNRJ_' + str(npts) + 'npts.npy')
        JetNRJ_Bal = np.load(savepath + '\Data_JetNRJ_Bal_' + str(npts) + 'npts.npy')
        SheetWide = np.load(savepath + '\Data_SheetWide_' + str(npts) + 'npts.npy')
        DispertionDist = np.load(savepath + '\Data_DispertionDist_' + str(npts) + 'npts.npy')   
        DispertionDist_Var = np.load(savepath + '\Data_DispertionDist_Var_' + str(npts) + 'npts.npy') 
        
    

        print('Done !')
         
    else:
      
        
        cpt = 0
        
        for a,dd,oc in zip(meshA.flatten(),meshDD.flatten(),meshOC.flatten()):
            
            idx = np.argwhere((meshA==a)&(meshDD==dd)&(meshOC==oc))
            
            cpt += 1
            if ConeSizeType == 'surface':        
            
                ConeDiam = 2*np.sqrt(ConeSize*np.sin(a)/np.pi)
                
            elif ConeSizeType == 'radius': 
                
                ConeDiam = ConeSize
            
    
            print(f'Computing impacts n°{cpt:1d} of {nsim:1d}.', end='\r')
                
            
            OffCmax = 0.495*(ConeDiam+dd)
            
                
            if oc < OffCmax:
                
                Cone = dgc.Cone(ConeDiam/2,a)                
                
                Drop = dgc.Drop(dd/2,oc,71,DropSpeed) # ~5000 points in the drop
                Impact = Cone.impact(Drop,oriType,velIni,meshType)
                JetFracs[idx[0][0],idx[0][1],idx[0][2]] = Impact.compute_JetFrac(velType)
                NRJ_tmp = Impact.compute_JetNRJ(velType)
                JetNRJ[idx[0][0],idx[0][1],idx[0][2]] = NRJ_tmp[0]
                JetNRJ_Bal[idx[0][0],idx[0][1],idx[0][2]] = NRJ_tmp[1]
                Dist_tmp = Impact.compute_DispertionDist(velType)
                DispertionDist[idx[0][0],idx[0][1],idx[0][2]] = Dist_tmp[0] 
                DispertionDist_Var[idx[0][0],idx[0][1],idx[0][2]] = Dist_tmp[1] 
                SheetWide[idx[0][0],idx[0][1],idx[0][2]] = Impact.SheetOpening(velType)[0]
                
                JetMass[idx[0][0],idx[0][1],idx[0][2]] = Impact.VolFrac/100*Impact.compute_JetFrac(velType)/100*Drop.Mass
            else:
                
                JetFracs[idx[0][0],idx[0][1],idx[0][2]] = np.nan
                SheetWide[idx[0][0],idx[0][1],idx[0][2]] = np.nan
                DispertionDist[idx[0][0],idx[0][1],idx[0][2]] = np.nan
                DispertionDist_Var[idx[0][0],idx[0][1],idx[0][2]] = np.nan
                JetMass[idx[0][0],idx[0][1],idx[0][2]] = np.nan
                JetNRJ[idx[0][0],idx[0][1],idx[0][2]] = np.nan
                JetNRJ_Bal[idx[0][0],idx[0][1],idx[0][2]] = np.nan
          
         
       
         
        np.save(savepath + '\Data_JetFracs_' + str(npts) + 'npts.npy',JetFracs)
        np.save(savepath + '\Data_JetMass_' + str(npts) + 'npts.npy',JetMass)
        np.save(savepath + '\Data_JetNRJ_' + str(npts) + 'npts.npy',JetNRJ)
        np.save(savepath + '\Data_JetNRJ_Bal_' + str(npts) + 'npts.npy',JetNRJ_Bal)
        np.save(savepath + '\Data_SheetWide_' + str(npts) + 'npts.npy',SheetWide)
        np.save(savepath + '\Data_DispertionDist_' + str(npts) + 'npts.npy',DispertionDist)
        np.save(savepath + '\Data_DispertionDist_Var_' + str(npts) + 'npts.npy',DispertionDist_Var)
              
        
    print('\n\nMaking and saving figures',end='\r')
    
    pointSize = 250000/npts**2  
                
    
    ########### Fixed angle diagrams ###########       

    
    for a,ia in zip(Angles,range(len(Angles))):
        
        # Impact fraction in jet
        f0,ax0 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax0.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
        ax0.set_xlabel('Offcent' + NormStr)
        ax0.set_ylabel('DropDiam' + NormStr)
    

        sc0 = ax0.scatter(2*meshOC[:,ia,:]/Normalisation,meshDD[:,ia,:]/Normalisation,c=JetFracs[:,ia,:],vmin=0, vmax = 100,
                          cmap='PuOr',marker='s',s=pointSize,zorder=2)

        cbar0 = plt.colorbar(sc0)
        cbar0.set_label('Impact fraction in the jet (%)')
        
        # ax0.set_aspect('equal')
        

        
        
        f0.tight_layout()
        
        f0.savefig(savepath + '\JetFrac\FixedAngle\FxdADgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_JetFrac.png')

        plt.close(f0)
        
        
        print('Making and saving figures     ',end='\r')

        # Sheet/Jet volume ratio
        f1,ax1 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax1.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
        ax1.set_xlabel('Offcent' + NormStr)
        ax1.set_ylabel('DropDiam' + NormStr)
        
        
        # sc1 = ax1.scatter(2*meshOC[:,ia,:]/Normalisation,meshDD[:,ia,:]/Normalisation,c=np.divide(np.subtract(100,JetFracs[:,ia,:]),JetFracs[:,ia,:])
        #                   ,cmap='plasma',s=pointSize)
        sc1 = ax1.scatter(2*meshOC[:,ia,:]/Normalisation,meshDD[:,ia,:]/Normalisation,c=np.divide(np.subtract(100,JetFracs[:,ia,:]),JetFracs[:,ia,:])
                          ,vmin=0, vmax = 3,cmap='plasma',s=pointSize,marker='s')
        
        cbar1 = plt.colorbar(sc1)
        cbar1.set_label('Sheet/Jet volume ratio')
        
        # ax1.set_aspect('equal')
        
        f1.tight_layout()
        
        f1.savefig(savepath + '\VolRatio\FixedAngle\FxdADgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_VolRatio.png')

        plt.close(f1)
        
        
        print('Making and saving figures.     ',end='\r')
        
        
        # Sheet opening
        f2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax2.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
        ax2.set_xlabel('Offcent' + NormStr)
        ax2.set_ylabel('DropDiam' + NormStr)
        
        sc2 = ax2.scatter(2*meshOC[:,ia,:]/Normalisation,meshDD[:,ia,:]/Normalisation,c=SheetWide[:,ia,:]*360/(2*np.pi),vmin=0, vmax = 360,
                          cmap='cividis',s=pointSize,marker='s')
        
        cbar2 = plt.colorbar(sc2)
        cbar2.set_label('Sheet opening [°]')
        # ax2.set_aspect('equal')
        f2.tight_layout()
        
        f2.savefig(savepath + '\SheetOpening\FixedAngle\FxdADgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_SheetOpen.png')

        plt.close(f2)
        
        
        print('Making and saving figures..     ',end='\r')
        
        # Dispertion distance
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
        ax3.set_xlabel('Offcent' + NormStr)
        ax3.set_ylabel('DropDiam' + NormStr)
        
        sc3 = ax3.scatter(2*meshOC[:,ia,:]/Normalisation,meshDD[:,ia,:]/Normalisation,c=DispertionDist[:,ia,:]/1000,vmin=0, cmap='inferno',s=pointSize,marker='s')
        
        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Average dispertion distance [m]')
        # ax3.set_aspect('equal')
        f3.tight_layout()
        
        f3.savefig(savepath + '\DispertionDist\FixedAngle\FxdADgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_DispertionDist.png')

        plt.close(f3)
        
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
        ax3.set_xlabel('Offcent' + NormStr)
        ax3.set_ylabel('DropDiam' + NormStr)

        sc3 = ax3.scatter(2*meshOC[:,ia,:]/Normalisation,meshDD[:,ia,:]/Normalisation,c=DispertionDist_Var[:,ia,:]/1000,cmap='plasma',
                          s=pointSize,marker='s')
        
        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Variability of dispertion dist [m]')
        # ax3.set_aspect('equal')
        f3.tight_layout()
        
        f3.savefig(savepath + '\DispertionDist_Var\FixedAngle\FxdADgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_DispertionDist_Var.png')

        plt.close(f3)
        
        
        print('Making and saving figures...     ',end='\r')
        
        # Kinetic energy
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
        ax3.set_xlabel('Offcent' + NormStr)
        ax3.set_ylabel('DropDiam' + NormStr)

        sc3 = ax3.scatter(2*meshOC[:,ia,:]/Normalisation,meshDD[:,ia,:]/Normalisation,c=JetNRJ[:,ia,:],cmap='rainbow',
                          s=pointSize,marker='s')
        
        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Kinetic energy [mJ]')
        # ax3.set_aspect('equal')
        f3.tight_layout()
        
        f3.savefig(savepath + '\EnRJ\FixedAngle\FxdADgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_KineticEnergy.png')

        plt.close(f3)
        
        
        # Balistic energy        
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
        ax3.set_xlabel('Offcent' + NormStr)
        ax3.set_ylabel('DropDiam' + NormStr)

        sc3 = ax3.scatter(2*meshOC[:,ia,:]/Normalisation,meshDD[:,ia,:]/Normalisation,c=JetNRJ_Bal[:,ia,:],cmap='jet',
                          s=pointSize,marker='s')
        
        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Balistic energy [mJ]')
        # ax3.set_aspect('equal')
        f3.tight_layout()
        
        f3.savefig(savepath + '\EnRJ_Bal\FixedAngle\FxdADgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_BalisticEnergy.png')

        plt.close(f3)
        

    
    ########### Fixed drop diagrams ###########
    
    for rdd,idd in zip(RelDropDiams,range(len(RelDropDiams))):
        
        # Impact fraction in jet
        f0,ax0 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax0.set_title('DropDiam_Norm : ' + str(round(rdd*10)/10))
        ax0.set_xlabel('Offcent' + NormStr)
        ax0.set_ylabel('Cone angle [°]')
        ax0.set_xlim([0,np.max(RelOffCents)])
        
        
        # sc0 = ax0.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=JetFracs[idd,:,:],cmap='PuOr',s=pointSize)
# 
        sc0 = ax0.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=JetFracs[idd,:,:],vmin=0, vmax = 100,cmap='PuOr'
                          ,marker='s',s=pointSize)

        cbar0 = plt.colorbar(sc0)
        cbar0.set_label('Impact fraction in the jet (%)')
        f0.tight_layout()
        
        f0.savefig(savepath + '\JetFrac\FixedDrop\FxdDrDgm_'
            + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SizeRatio_JetFrac.png')

        plt.close(f0)
        
        
        print('Making and saving figures     ',end='\r')


        # Sheet/Jet volume ratio
        f1,ax1 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax1.set_title('DropDiam_Norm : ' + str(round(rdd*10)/10))
        ax1.set_xlabel('Offcent' + NormStr)
        ax1.set_ylabel('Cone angle [°]')
        ax1.set_xlim([0,np.max(RelOffCents)])
        
        
        # sc1 = ax1.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=np.divide(np.subtract(100,JetFracs[idd,:,:]),JetFracs[idd,:,:]),
                          # cmap='plasma',s=pointSize)
        
        sc1 = ax1.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=np.divide(np.subtract(100,JetFracs[idd,:,:]),JetFracs[idd,:,:]),
                          vmin=0, vmax = 3,cmap='plasma',s=pointSize,marker='s')
        
        cbar1 = plt.colorbar(sc1)
        cbar1.set_label('Sheet/Jet volume ratio')
        f1.tight_layout()
        
        f1.savefig(savepath + '\VolRatio\FixedDrop\FxdDrDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SizeRatio_VolRatio.png')
        
        plt.close(f1)
        
        
        print('Making and saving figures.     ',end='\r')
        

        # Sheet opening
        f2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax2.set_title('DropDiam_Norm : ' + str(round(rdd*10)/10))
        ax2.set_xlabel('Offcent' + NormStr)
        ax2.set_ylabel('Cone angle [°]')
        ax2.set_xlim([0,np.max(RelOffCents)])
        
        # sc2 = ax2.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=SheetWide[idd,:,:]*360/(2*np.pi),cmap='cividis',s=pointSize)
        
        sc2 = ax2.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=SheetWide[idd,:,:]*360/(2*np.pi),vmin=0, 
                          vmax = 360,cmap='cividis',s=pointSize,marker='s')
        
        cbar2 = plt.colorbar(sc2)
        cbar2.set_label('Sheet opening [°]')
        f2.tight_layout()
        
        f2.savefig(savepath + '\SheetOpening\FixedDrop\FxdDrDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SizeRatio_SheetOpen.png')
        
        plt.close(f2)
        
        
        print('Making and saving figures..     ',end='\r')
        
        
        # dispersal distance
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('DropDiam_Norm : ' + str(round(rdd*10)/10))
        ax3.set_xlabel('Offcent' + NormStr)
        ax3.set_ylabel('Cone angle [°]')
        ax3.set_xlim([0,np.max(RelOffCents)])
        # 
        # sc3 = ax3.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=DispertionDist[idd,:,:]*1000,cmap='jet',s=pointSize)
        
        sc3 = ax3.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=DispertionDist[idd,:,:]/1000,vmin=0,cmap='inferno',s=pointSize,marker='s')
        
        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Average dispertion distance [m]')
        f3.tight_layout()
        
        f3.savefig(savepath + '\DispertionDist\FixedDrop\FxdDrDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SizeRatio_DispertionDist.png')
        
        plt.close(f3)
        
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('DropDiam_Norm : ' + str(round(rdd*10)/10))
        ax3.set_xlabel('Offcent' + NormStr)
        ax3.set_ylabel('Cone angle [°]')
        ax3.set_xlim([0,np.max(RelOffCents)])
        # 
        # sc3 = ax3.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=DispertionDist[idd,:,:]*1000,cmap='jet',s=pointSize)
        
        sc3 = ax3.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=DispertionDist_Var[idd,:,:]/1000,vmin=0,
                          cmap='plasma',s=pointSize,marker='s')
        
        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Variability of dispertion dist [m]')
        f3.tight_layout()
        
        f3.savefig(savepath + '\DispertionDist_Var\FixedDrop\FxdDrDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SizeRatio_DispertionDist_Var.png')
        
        plt.close(f3)
        
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('DropDiam_Norm : ' + str(round(rdd*10)/10))
        ax3.set_xlabel('Offcent' + NormStr)
        ax3.set_ylabel('Cone angle [°]')
        ax3.set_xlim([0,np.max(RelOffCents)])
        # 
        # sc3 = ax3.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=DispertionDist[idd,:,:]*1000,cmap='jet',s=pointSize)
        
        sc3 = ax3.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=JetNRJ[idd,:,:],vmin=0,
                          cmap='rainbow',s=pointSize,marker='s')
        
        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Kinetic energy [mJ]')
        f3.tight_layout()
        
        f3.savefig(savepath + '\EnRJ\FixedDrop\FxdDrDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SizeRatio_NRJ.png')
        
        plt.close(f3)
        
        
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('DropDiam_Norm : ' + str(round(rdd*10)/10))
        ax3.set_xlabel('Offcent' + NormStr)
        ax3.set_ylabel('Cone angle [°]')
        ax3.set_xlim([0,np.max(RelOffCents)])
        # 
        # sc3 = ax3.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=DispertionDist[idd,:,:]*1000,cmap='jet',s=pointSize)
        
        sc3 = ax3.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=JetNRJ_Bal[idd,:,:],vmin=0,
                          cmap='jet',s=pointSize,marker='s')
        
        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Balistic energy [mJ]')
        f3.tight_layout()
        
        f3.savefig(savepath + '\EnRJ_Bal\FixedDrop\_1_FxdDrDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SizeRatio_EnRJ_Bal.png')
        
        plt.close(f3)
        
        
        print('Making and saving figures...     ',end='\r')
        
        
        
    ########### Fixed Offcent diagrams ###########
    
    for roc,ioc in zip(RelOffCents,range(len(RelOffCents))):
    
        # Impact fraction in jet
        f0,ax0 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax0.set_title('Offcent_Norm: ' + str(round(roc*10)/10))
        ax0.set_xlabel('DropDiam' + NormStr)
        ax0.set_ylabel('Cone angle [°]')
        ax0.set_xlim([0,np.max(RelDropDiams)])
        
        
        # sc0 = ax0.scatter(meshDD[:,:,ioc]/Normalisation,meshA[:,:,ioc]/(2*np.pi)*360,c=JetFracs[:,:,ioc],cmap='PuOr',s=pointSize)
    
        sc0 = ax0.scatter(meshDD[:,:,ioc]/Normalisation,meshA[:,:,ioc]/(2*np.pi)*360,c=JetFracs[:,:,ioc],vmin=0, vmax = 100,cmap='PuOr',marker='s',s=pointSize)
    
        cbar0 = plt.colorbar(sc0)
        cbar0.set_label('Impact fraction in the jet (%)')
        f0.tight_layout()
        
        f0.savefig(savepath + '\JetFrac\FixedDist\FxdDiDgm_'
          + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OffCent_JetFrac.png')
  
        plt.close(f0)
        
        
        print('Making and saving figures     ',end='\r')

        # Sheet/Jet volume ratio
        f1,ax1 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax1.set_title('Offcent_Norm: ' + str(round(roc*10)/10))
        ax1.set_xlabel('DropDiam' + NormStr)
        ax1.set_ylabel('Cone angle [°]')
        ax1.set_xlim([0,np.max(RelDropDiams)])
        
        
        # sc1 = ax1.scatter(meshDD[:,:,ioc]/Normalisation,meshA[:,:,ioc]/(2*np.pi)*360,c=np.divide(np.subtract(100,JetFracs[:,:,ioc]),JetFracs[:,:,ioc]),
        #                   cmap='plasma',s=pointSize)
        
        sc1 = ax1.scatter(meshDD[:,:,ioc]/Normalisation,meshA[:,:,ioc]/(2*np.pi)*360,c=np.divide(np.subtract(100,JetFracs[:,:,ioc]),JetFracs[:,:,ioc]),
                          vmin=0, vmax = 3,cmap='plasma',marker='s',s=pointSize)
        
        cbar1 = plt.colorbar(sc1)
        cbar1.set_label('Sheet/Jet volume ratio')
        f1.tight_layout()
        
        f1.savefig(savepath + '\VolRatio\FixedDist\FxdDiDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OffCent_VolRatio.png')
        
        plt.close(f1)
        
        
        print('Making and saving figures.     ',end='\r')
        
        # Sheet/Jet volume ratio
        f2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax2.set_title('Offcent_Norm: ' + str(round(roc*10)/10))
        ax2.set_xlabel('DropDiam' + NormStr)
        ax2.set_ylabel('Cone angle [°]')
        ax2.set_xlim([0,np.max(RelDropDiams)])
        
        # sc2 = ax2.scatter(meshDD[:,:,ioc]/Normalisation,meshA[:,:,ioc]/(2*np.pi)*360,c=SheetWide[:,:,ioc]*360/(2*np.pi),cmap='cividis',s=pointSize)
        
        sc2 = ax2.scatter(meshDD[:,:,ioc]/Normalisation,meshA[:,:,ioc]/(2*np.pi)*360,c=SheetWide[:,:,ioc]*360/(2*np.pi),vmin=0, vmax = 360,
                          cmap='cividis',s=pointSize,marker='s')
        
        cbar2 = plt.colorbar(sc2)
        cbar2.set_label('Sheet opening [°]')
        f2.tight_layout()
        
        f2.savefig(savepath + '\SheetOpening\FixedDist\FxdDiDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OffCent_SheetOpen.png')
        
        plt.close(f2)
        
        
        print('Making and saving figures..     ',end='\r')
        
        # kinetic energy in the jet
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('Offcent_Norm: ' + str(round(roc*10)/10))
        ax3.set_xlabel('DropDiam' + NormStr)
        ax3.set_ylabel('Cone angle [°]')
        ax3.set_xlim([0,np.max(RelDropDiams)])
        
        sc3 = ax3.scatter(meshDD[:,:,ioc]/Normalisation,meshA[:,:,ioc]/(2*np.pi)*360,c=DispertionDist[:,:,ioc]/1000,vmin = 0,cmap='inferno',s=pointSize,marker='s')
        
        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Average dispertion distance [m]')
        f3.tight_layout()
        
        f3.savefig(savepath + '\DispertionDist\FixedDist\FxdDiDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OffCent_DispertionDist.png')
        
        plt.close(f3)
        
        # kinetic energy ratio in the jet
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('Offcent_Norm: ' + str(round(roc*10)/10))
        ax3.set_xlabel('DropDiam' + NormStr)
        ax3.set_ylabel('Cone angle [°]')
        ax3.set_xlim([0,np.max(RelDropDiams)])
        
        sc3 = ax3.scatter(meshDD[:,:,ioc]/Normalisation,meshA[:,:,ioc]/(2*np.pi)*360,c=DispertionDist_Var[:,:,ioc]/1000,
                          cmap='plasma',s=pointSize,marker='s')
        

        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Variability of dispertion dist [m]')
        f3.tight_layout()
        
        f3.savefig(savepath + '\DispertionDist_Var\FixedDist\FxdDiDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OffCent_DispertionDist_Var.png')
        
        plt.close(f3)
        
        # kinetic energy ratio in the jet
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('Offcent_Norm: ' + str(round(roc*10)/10))
        ax3.set_xlabel('DropDiam' + NormStr)
        ax3.set_ylabel('Cone angle [°]')
        ax3.set_xlim([0,np.max(RelDropDiams)])
        
        sc3 = ax3.scatter(meshDD[:,:,ioc]/Normalisation,meshA[:,:,ioc]/(2*np.pi)*360,c=JetNRJ[:,:,ioc],
                          cmap='rainbow',s=pointSize,marker='s')
        

        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Kinetic energy [mJ]')
        f3.tight_layout()
        
        f3.savefig(savepath + '\EnRJ\FixedDist\FxdDiDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OffCent_NRJ.png')
        
        plt.close(f3)
        
        # balistic energy ratio in the jet
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('Offcent_Norm: ' + str(round(roc*10)/10))
        ax3.set_xlabel('DropDiam' + NormStr)
        ax3.set_ylabel('Cone angle [°]')
        ax3.set_xlim([0,np.max(RelDropDiams)])
        
        sc3 = ax3.scatter(meshDD[:,:,ioc]/Normalisation,meshA[:,:,ioc]/(2*np.pi)*360,c=JetNRJ_Bal[:,:,ioc]/10000/2*9.81e-3,
                          cmap='jet',s=pointSize,marker='s')
        

        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Balistic energy [mJ]')
        f3.tight_layout()
        
        f3.savefig(savepath + '\EnRJ_Bal\FixedDist\_1_FxdDiDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OffCent_EnRJ_Bal.png')
        
        plt.close(f3)
        
    print('\nDone.')
    
    return
        
        



###
# 3. Cone optimization diagrams

def OptiDiagrams(ConeSizes,sizeType,ConeAngles,oriType,velIni,meshType,npts,ndrops,dropDist,label):
    
    savepath = r'd:\Users\laplaud\Desktop\PostDoc\Code\DropProject_WithAna\Figures\Optimization/' + label + '_' + str(npts) + '_' + str(ndrops)

    coneAngles = np.linspace(ConeAngles[0],ConeAngles[1],npts)/360*2*np.pi

    coneSizes = np.linspace(ConeSizes[0],ConeSizes[1],npts)
    
    if sizeType == 'area':
        
        coneLabel = 'Cone surface area (mm²)'
        
    elif sizeType == 'side':
        
        coneLabel = 'Cone side length (mm)'
            
    elif sizeType == 'diameter':
        
        coneLabel = 'Cone diameter (mm)'
        

    meshCA,meshCS = np.meshgrid(coneAngles,coneSizes)
    
    if os.path.exists(savepath):
        
        print('Loading previous simulations results...', end = '')
        
        
        impactVolumes = np.load(savepath + '\Data_impactVolumes.npy')
        jetVolumes = np.load(savepath + '\Data_jetVolumes.npy')
        efficiency = np.load(savepath + '\Data_efficiency.npy')
        NRJefficiency = np.load(savepath + '\Data_NRJefficiency.npy')
        jetNRJs = np.load(savepath + '\Data_jetNRJs.npy')
        jetNRJsBalis = np.load(savepath + '\Data_jetNRJsBalis.npy')
        TotalVolume = np.load(savepath + '\Data_TotalVolume.npy')

        print('Done !')
         
    else:
        os.makedirs(savepath) # create folder

        jetVolumes = np.empty(np.shape(meshCA))
        impactVolumes = np.empty(np.shape(meshCA))
        efficiency = np.empty(np.shape(meshCA))
        NRJefficiency = np.empty(np.shape(meshCA))
        jetNRJs = np.empty(np.shape(meshCS))
        jetNRJsBalis = np.empty(np.shape(meshCS))

        # Drop distribution

        rho = 1000 # in [kg/m^3]

        if dropDist == 'exp':
            
            dropSizes = np.random.exponential(size=ndrops,scale= 1) # 
            
        
        if dropDist == 'gamma':
            
            dropSizes = np.random.gamma(3,0.5, ndrops)

        dropSizes[dropSizes>5] = dropSizes[dropSizes>5]-5 # Max radius [mm]

        dropVels = np.sqrt(8/3*1000/1.3*10*dropSizes/500) # in m/s or mm/ms

        if sizeType == 'area':
            
            dropRs = np.sqrt(np.random.rand(ndrops)*(np.sqrt(np.max(coneSizes)/np.pi*np.sin(np.max(coneAngles)))+np.max(dropSizes))**2)
            
        elif sizeType == 'side':
            
            dropRs = np.sqrt(np.random.rand(ndrops)*(np.max(coneSizes)*np.sin(np.max(coneAngles))+np.max(dropSizes))**2)
                
        elif sizeType == 'diameter':

            dropRs = np.sqrt(np.random.rand(ndrops)*(np.max(coneSizes)+np.max(dropSizes))**2)
        

        dropAs = np.random.rand(ndrops)*2*np.pi
        
        TotalVolume = np.sum(4/3*np.pi*dropSizes**3)/1000 # in [cm3] 

        # Drop and cones plot
        fig,ax = plt.subplots(dpi=200)
        dropXs,dropYs = vf.ToCart(dropAs,dropRs,angle='rad')
        tx = np.linspace(0,2*np.pi,30)
        for cr in coneSizes:
            ax.plot(np.cos(tx)*cr*np.tan(np.pi/2),np.sin(tx)*cr*np.tan(np.pi/2),'g')
        ax.plot(0,0,'go')
        for x,y,r in zip(dropXs[0:50:],dropYs[0:50:],dropSizes[0:50:]):
            ax.plot(x,y,'b.',ms=1)
            ax.plot(x+r*np.cos(tx),y+r*np.sin(tx),'c-',lw=1)

        ax.set_aspect('equal')
        
        fig.savefig(savepath + '\DropAndCone.png')

        plt.close(fig)
        
        f,ax = plt.subplots(dpi=150,figsize = (7,6)) 
        ax.hist(dropSizes,density = True,bins=20,label='exponential corrected (max = 5)')
        ax.set_title('PDF of drop radius')
        ax.set_xlabel('Drop radius (mm)')
        ax.legend()
        
        f.savefig(savepath + '\DropSizes.png')

        plt.close(f)

        f,ax = plt.subplots(dpi=150,figsize = (7,6)) 
        ax.plot(dropSizes,dropVels,'o')
        ax.set_xlabel('Drop radius (mm)')
        ax.set_ylabel('Drop speed (m/s)')
        
        f.savefig(savepath + '\DropSizeVsSpeed.png')

        plt.close(f)

        f,ax = plt.subplots(dpi=150,figsize = (7,6)) 
        ax.hist(dropVels)
        ax.set_title('Drop speeds (m/s)')
        ax.set_xlabel('Drop speed (m/s)')
        
        f.savefig(savepath + '\DropSpeeds.png')

        plt.close(f)


        ###### Simulations


        for a,ia in zip(coneAngles,range(len(coneAngles))):
            for cs,ics in zip(coneSizes,range(len(coneSizes))):

                coneNum = ia*len(coneSizes) + ics + 1

                jetVolume = np.empty(np.shape(dropSizes))
                impactVolume = np.empty(np.shape(dropSizes))
                jetNRJ = np.empty(np.shape(dropSizes))
                TotalTouchingNRJ = 0
                
                if sizeType == 'area':
                    
                    cr = np.sqrt(cs/np.pi*np.sin(a)) # cone Radius
                    
                elif sizeType == 'side':
                    
                    cr = cs*np.sin(a)
                        
                elif sizeType == 'diameter':
                    
                    cr = cs/2
                

                for ds,dr,dv,di in zip(dropSizes,dropRs,dropVels,range(len(dropVels))):

                    print("Computing impacts for cone " + str(coneNum) + "/" + str(len(coneSizes)*len(coneAngles)) + " : Impact n°"+str(di)+"/"+str(len(dropSizes))+".".ljust(10),end='\r')
                    

                    
                    if dr<(cr+ds)*0.95:
                        
                        I = dgc.Cone(cr,a).impact(dgc.Drop(ds,dr,71,dv),oriType,velIni,meshType)
                        
                        TotalTouchingNRJ = TotalTouchingNRJ + 4/3*np.pi*(ds/1000)**3*rho/2*dv**2
    
                        impactVolume[di] = I.VolFrac/100*4/3*np.pi*(ds/1000)**3
    
                        jetVolume[di] = impactVolume[di]*I.compute_JetFrac('full_div0')/100
    
                        jetNRJ[di] = jetVolume[di]*rho/2*dv**2
                        
                        
                        del I
                        
                    else:
                        
                        impactVolume[di] = 0
    
                        jetVolume[di] = 0
    
                        jetNRJ[di] = 0

                    

                impactVolumes[ics,ia] = np.sum(impactVolume) # in [m^3]

                jetVolumes[ics,ia] = np.sum(jetVolume) # in [m^3]

                efficiency[ics,ia] = np.sum(jetVolume)/np.sum(impactVolume)*100 # in %

                jetNRJs[ics,ia] = np.sum(jetNRJ) # [J]
                
                NRJefficiency[ics,ia] = np.sum(jetNRJ)/TotalTouchingNRJ # [J]
                
                jetNRJsBalis[ics,ia] = np.sum(jetNRJ)*np.sin(2*a) # [balistic J]
                
                
                
        np.save(savepath + '\Data_impactVolumes.npy',impactVolumes)
        np.save(savepath + '\Data_jetVolumes.npy',jetVolumes)
        np.save(savepath + '\Data_efficiency.npy',efficiency)
        np.save(savepath + '\Data_NRJefficiency.npy',NRJefficiency)
        np.save(savepath + '\Data_jetNRJs.npy',jetNRJs)
        np.save(savepath + '\Data_TotalVolume.npy',TotalVolume)

    # ############# Optimization diagrams

    pointSize = 250000/(npts**2) 

     # Impact volume
    f0,ax0 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax0.set_title('Impact volume for random rain of ' + str(ndrops) + ' drops')
    ax0.set_xlabel('Cone angle (°)')
    ax0.set_ylabel(coneLabel)

    sc0 = ax0.scatter(meshCA/(2*np.pi)*360,meshCS,c=impactVolumes*1e6/TotalVolume*100,vmin=0,cmap='viridis',s=pointSize,marker='s')

    cbar0 = plt.colorbar(sc0)
    cbar0.set_label('Total volume impacting the cone (% of total rain volume)')
    f0.tight_layout()

    f0.savefig(savepath + '\OptiDgm_'
                + label + '_'+str(int(npts))+'npts_ImpactVolume.png')

    plt.close(f0)

     # Jet volume
    f1,ax1 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax1.set_title('Jet volume for random rain of ' + str(ndrops) + ' drops')
    ax1.set_xlabel('Cone angle (°)')
    ax1.set_ylabel(coneLabel)

    sc1 = ax1.scatter(meshCA/(2*np.pi)*360,meshCS,c=jetVolumes*1e6/TotalVolume*100,vmin=0,cmap='PuOr',s=pointSize,marker='s')

    cbar1 = plt.colorbar(sc1)
    cbar1.set_label('Total volume ejected as jets (% of total rain volume)')
    f1.tight_layout()

    f1.savefig(savepath + '\OptiDgm_'
                + label + '_'+str(int(npts))+'npts_JetVolume.png')

    plt.close(f1)

    # Jet/Impact volume
    f2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax2.set_title('Efficiency [jet/impact volumes] for random rain of ' + str(ndrops) + ' drops')
    ax2.set_xlabel('Cone angle (°)')
    ax2.set_ylabel(coneLabel)

    sc2 = ax2.scatter(meshCA/(2*np.pi)*360,meshCS,c=efficiency,vmin=0,vmax = 100,cmap='jet',s=pointSize,marker='s')
        

    cbar2 = plt.colorbar(sc2)
    cbar2.set_label('Efficiency [jet/impact volumes] (%)')
    f2.tight_layout()

    f2.savefig(savepath + '\OptiDgm_'
                + label + '_'+str(int(npts))+'npts_Efficiency.png')

    plt.close(f2)
    
    # Jet/touching nrj
    f2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax2.set_title('NRJ efficiency [jet/touching] for random rain of ' + str(ndrops) + ' drops')
    ax2.set_xlabel('Cone angle (°)')
    ax2.set_ylabel(coneLabel)

    sc2 = ax2.scatter(meshCA/(2*np.pi)*360,meshCS,c=NRJefficiency*100,vmin=0,vmax = 100,cmap='jet',s=pointSize,marker='s')
        

    cbar2 = plt.colorbar(sc2)
    cbar2.set_label('NRJ efficiency [jet/touching] (%)')
    f2.tight_layout()

    f2.savefig(savepath + '\OptiDgm_'
                + label + '_'+str(int(npts))+'npts_NRJefficiency.png')

    plt.close(f2)

     # Kinetic energy
    f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax3.set_title('Maximum kinetic energy in jets\n for random rain of ' + str(ndrops) + ' drops')
    ax3.set_xlabel('Cone angle (°)')
    ax3.set_ylabel(coneLabel)

    sc3 = ax3.scatter(meshCA/(2*np.pi)*360,meshCS,c=jetNRJs,cmap='jet',s=pointSize,marker='s')

    cbar3 = plt.colorbar(sc3)
    cbar3.set_label('Kinetic energy of the jets (J)')
    f3.tight_layout()

    f3.savefig(savepath + '\OptiDgm_'
                + label + '_'+str(int(npts))+'npts_KineticEnergy.png')

    plt.close(f3)
    
    
     # Balistic energy
    f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax3.set_title('Maximum balistic energy in jets\n for random rain of ' + str(ndrops) + ' drops')
    ax3.set_xlabel('Cone angle (°)')
    ax3.set_ylabel(coneLabel)

    sc3 = ax3.scatter(meshCA/(2*np.pi)*360,meshCS,c=jetNRJsBalis,cmap='jet',s=pointSize,marker='s')

    cbar3 = plt.colorbar(sc3)
    cbar3.set_label('Balistic energy of the jets (J)')
    f3.tight_layout()

    f3.savefig(savepath + '\OptiDgm_'
                + label + '_'+str(int(npts))+'npts_BalisticEnergy.png')

    plt.close(f3)
    
    
    print('\n\nFigures ploted and saved !')

