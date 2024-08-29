# -*- coding: utf-8 -*-
"""
Created on Mon Oct 30 08:04:43 2023

@author: laplaud
"""

import numpy as np

import matplotlib.pyplot as plt

import DropGeometryClasses as dgc

import os

import shutil


###
# 1. Plotting function for different volume fractions
def plotFracs(Angle,npts,DropDiam,ConeDiams,oriType,velIni,meshType):
    
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
        JetFracs = [I.compute_JetFrac() for I in Impacts]
        DropFracs = [I.compute_JetFrac()*I.VolFrac/100 for I in Impacts]

        lab = 'Drop size / cone size = ' + str(np.round(100*DropDiam/cd)/100)
        
        ax3.plot(OffCents/cr,DropFracs,'-*',label=lab)
        ax1.plot(OffCents/cr,JetFracs,'-*',label=lab)
        ax2.plot(OffCents/cr,np.divide(np.subtract(100,JetFracs),JetFracs),'-*',label=lab)

        ax3.legend()
        ax1.legend()
        ax2.legend()
     
###
# 2. Parameter space diagrams

def PhaseDiagrams(RelOffCents,ConeSize,ConeSizeType,Angles,RelDropDiams,oriType,velIni,velType,meshType,path,label):
    
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
    SecJetFracs = np.empty(np.shape(meshA))
    JetMass = np.empty(np.shape(meshA))
    JetNRJ = np.empty(np.shape(meshA))
    JetVel = np.empty(np.shape(meshA))
    SheetWide = np.empty(np.shape(meshA))
    ShapeFactor = np.empty(np.shape(meshA))
    DispertionDist = np.empty(np.shape(meshA))
    DispertionDistVar = np.empty(np.shape(meshA))
    
    
    nsim = meshDD.size
    
    savepath = path + label + '_' + str(npts) + 'npts'
    
    
    if os.path.exists(savepath+ '\JetFrac'):
        
        print('Clearing previous figure folder...', end = '')
        
        shutil.rmtree(savepath + '\JetFrac\\')   
        shutil.rmtree(savepath + '\SecJetFrac\\')   
        shutil.rmtree(savepath + '\VolRatio\\')   
        shutil.rmtree(savepath + '\JetVel\\')   
        shutil.rmtree(savepath + '\EnRJ\\')   
        shutil.rmtree(savepath + '\DispertionDist\\')   
        shutil.rmtree(savepath + '\DispertionDistVar\\')   
        shutil.rmtree(savepath + '\SheetOpening\\')  
        shutil.rmtree(savepath + '\ShapeFactor\\') 
        
        print('Done.\n')
        

    print('Creating figures folders...', end = '')
    
    os.makedirs(savepath + '\JetFrac\FixedAngle',exist_ok=True) # create folder
    os.makedirs(savepath + '\JetFrac\FixedDrop',exist_ok=True) # create folder
    os.makedirs(savepath + '\JetFrac\FixedDist',exist_ok=True) # create folder
    
    os.makedirs(savepath + '\SecJetFrac\FixedAngle',exist_ok=True) # create folder
    os.makedirs(savepath + '\SecJetFrac\FixedDrop',exist_ok=True) # create folder
    os.makedirs(savepath + '\SecJetFrac\FixedDist',exist_ok=True) # create folder
    
    os.makedirs(savepath + '\BothJetFrac\FixedAngle',exist_ok=True) # create folder
    os.makedirs(savepath + '\BothJetFrac\FixedDrop',exist_ok=True) # create folder
    os.makedirs(savepath + '\BothJetFrac\FixedDist',exist_ok=True) # create folder

    os.makedirs(savepath + '\VolRatio\FixedAngle',exist_ok=True) # create folder
    os.makedirs(savepath + '\VolRatio\FixedDrop',exist_ok=True) # create folder
    os.makedirs(savepath + '\VolRatio\FixedDist',exist_ok=True) # create folder

    os.makedirs(savepath + '\JetVel\FixedAngle',exist_ok=True) # create folder
    os.makedirs(savepath + '\JetVel\FixedDrop',exist_ok=True) # create folder
    os.makedirs(savepath + '\JetVel\FixedDist',exist_ok=True) # create folder

    os.makedirs(savepath + '\EnRJ\FixedAngle',exist_ok=True) # create folder
    os.makedirs(savepath + '\EnRJ\FixedDrop',exist_ok=True) # create folder
    os.makedirs(savepath + '\EnRJ\FixedDist',exist_ok=True) # create folder
     
    os.makedirs(savepath + '\DispertionDist\FixedAngle',exist_ok=True) # create folder
    os.makedirs(savepath + '\DispertionDist\FixedDrop',exist_ok=True) # create folder
    os.makedirs(savepath + '\DispertionDist\FixedDist',exist_ok=True) # create folder
     
    os.makedirs(savepath + '\DispertionDistVar\FixedAngle',exist_ok=True) # create folder
    os.makedirs(savepath + '\DispertionDistVar\FixedDrop',exist_ok=True) # create folder
    os.makedirs(savepath + '\DispertionDistVar\FixedDist',exist_ok=True) # create folder
     
    os.makedirs(savepath + '\SheetOpening\FixedAngle',exist_ok=True) # create folder
    os.makedirs(savepath + '\SheetOpening\FixedDrop',exist_ok=True) # create folder
    os.makedirs(savepath + '\SheetOpening\FixedDist',exist_ok=True) # create folder
     
    os.makedirs(savepath + '\ShapeFactor\FixedAngle',exist_ok=True) # create folder
    os.makedirs(savepath + '\ShapeFactor\FixedDrop',exist_ok=True) # create folder
    os.makedirs(savepath + '\ShapeFactor\FixedDist',exist_ok=True) # create folder
        
            
    print('Done.\n')
    
    if os.path.exists(savepath) & os.path.exists(savepath + '\Data_DispertionDistVar_' + str(npts) + 'npts.npy'):
        
        print('Loading previous simulations results...', end = '')
        
        JetFracs = np.load(savepath + '\Data_JetFracs_' + str(npts) + 'npts.npy')
        SecJetFracs = np.load(savepath + '\Data_SecJetFracs_' + str(npts) + 'npts.npy')
        JetMass = np.load(savepath + '\Data_JetMass_' + str(npts) + 'npts.npy')
        JetNRJ = np.load(savepath + '\Data_JetNRJ_' + str(npts) + 'npts.npy')
        JetVel = np.load(savepath + '\Data_JetVel_' + str(npts) + 'npts.npy')
        SheetWide = np.load(savepath + '\Data_SheetWide_' + str(npts) + 'npts.npy')
        ShapeFactor = np.load(savepath + '\Data_ShapeFactor_' + str(npts) + 'npts.npy')
        DispertionDist = np.load(savepath + '\Data_DispertionDist_' + str(npts) + 'npts.npy')   
        DispertionDistVar = np.load(savepath + '\Data_DispertionDistVar_' + str(npts) + 'npts.npy') 

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
                
                Drop = dgc.Drop(dd/2,oc,72,DropSpeed) # ~5000 points in the drop
                Impact = Cone.impact(Drop,oriType,velIni,velType,meshType)
                
                
                JetFracs[idx[0][0],idx[0][1],idx[0][2]] = Impact.get_JetFrac()[0]
                SecJetFracs[idx[0][0],idx[0][1],idx[0][2]] = Impact.get_JetFrac()[1]
                NRJ_tmp = Impact.compute_JetNRJ()
                JetNRJ[idx[0][0],idx[0][1],idx[0][2]] = NRJ_tmp[0]
                JetVel[idx[0][0],idx[0][1],idx[0][2]] = Impact.compute_JetVel()
                Dist_tmp = Impact.compute_DispertionDist()
                DispertionDist[idx[0][0],idx[0][1],idx[0][2]] = Dist_tmp[0] 
                DispertionDistVar[idx[0][0],idx[0][1],idx[0][2]] = Dist_tmp[1] 
                SheetWide[idx[0][0],idx[0][1],idx[0][2]] = Impact.SheetOpening()[0]
                ShapeFactor[idx[0][0],idx[0][1],idx[0][2]] = Impact.compute_ShapeFactor()
                
                JetMass[idx[0][0],idx[0][1],idx[0][2]] = Impact.VolFrac/100*Impact.get_JetFrac()[0]/100*Drop.Mass
                

            else:
                
                JetFracs[idx[0][0],idx[0][1],idx[0][2]] = np.nan
                SecJetFracs[idx[0][0],idx[0][1],idx[0][2]] = np.nan
                SheetWide[idx[0][0],idx[0][1],idx[0][2]] = np.nan
                ShapeFactor[idx[0][0],idx[0][1],idx[0][2]] = np.nan
                DispertionDist[idx[0][0],idx[0][1],idx[0][2]] = np.nan
                DispertionDistVar[idx[0][0],idx[0][1],idx[0][2]] = np.nan
                JetMass[idx[0][0],idx[0][1],idx[0][2]] = np.nan
                JetNRJ[idx[0][0],idx[0][1],idx[0][2]] = np.nan
                JetVel[idx[0][0],idx[0][1],idx[0][2]] = np.nan
          
         
       
        ShapeFactor[np.isnan(ShapeFactor)]=0 
        np.save(savepath + '\Data_JetFracs_' + str(npts) + 'npts.npy',JetFracs)
        np.save(savepath + '\Data_SecJetFracs_' + str(npts) + 'npts.npy',SecJetFracs)
        np.save(savepath + '\Data_JetMass_' + str(npts) + 'npts.npy',JetMass)
        np.save(savepath + '\Data_JetNRJ_' + str(npts) + 'npts.npy',JetNRJ)
        np.save(savepath + '\Data_JetVel_' + str(npts) + 'npts.npy',JetVel)
        np.save(savepath + '\Data_SheetWide_' + str(npts) + 'npts.npy',SheetWide)
        np.save(savepath + '\Data_ShapeFactor_' + str(npts) + 'npts.npy',ShapeFactor)
        np.save(savepath + '\Data_DispertionDist_' + str(npts) + 'npts.npy',DispertionDist)
        np.save(savepath + '\Data_DispertionDistVar_' + str(npts) + 'npts.npy',DispertionDistVar)
              
        
    print('\n\nMaking and saving figures',end='\r')
    
    pointSize = 250000/npts**2  
                
    ########### Maximums by angle    ###########
    
    JFAnglesMaxes = np.nanmax(JetFracs,axis=(0,2))
    
    fa,axa = plt.subplots(dpi=200)
    axa.set_title('Maximum jet fraction per angle')
    axa.set_xlabel('Angle (°)')
    axa.set_ylabel('Max jetFrac (%)')
    axa.plot(Angles/(2*np.pi)*360,JFAnglesMaxes,'-ob',ms = 2,lw = 1.5)
    
    fa.tight_layout()
    
    fa.savefig(savepath + '\FxdA_'+ label + '_'+str(int(npts))+'npts_MaxCurve_JetFrac.png')

    plt.close(fa)
    
    SJFAnglesMaxes = np.nanmax(SecJetFracs,axis=(0,2))
    
    fa,axa = plt.subplots(dpi=200)
    axa.set_title('Maximum second jet fraction per angle')
    axa.set_xlabel('Angle (°)')
    axa.set_ylabel('Max second jetFrac (%)')
    axa.plot(Angles/(2*np.pi)*360,SJFAnglesMaxes,'-og',ms = 2,lw = 1.5)
    
    fa.tight_layout()
    
    fa.savefig(savepath + '\FxdA_'+ label + '_'+str(int(npts))+'npts_MaxCurve_SecJetFrac.png')

    plt.close(fa)
    
    
    
    DDAnglesMaxes = np.nanmax(DispertionDist/1000,axis=(0,2))
    
    fa,axa = plt.subplots(dpi=200)
    axa.set_title('Maximum dispertion distance per angle')
    axa.set_xlabel('Angle (°)')
    axa.set_ylabel('Max dispertion (m)')
    axa.plot(Angles/(2*np.pi)*360,DDAnglesMaxes,'-or',ms = 2,lw = 1.5)
    
    fa.tight_layout()
    
    fa.savefig(savepath + '\FxdA_'+ label + '_'+str(int(npts))+'npts_MaxCurve_DispDist.png')

    plt.close(fa)
    
    JNRJAnglesMaxes = np.nanmax(JetNRJ,axis=(0,2))
    
    fa,axa = plt.subplots(dpi=200)
    axa.set_title('Maximum kinetic energy per angle')
    axa.set_xlabel('Angle (°)')
    axa.set_ylabel('Max energy')
    axa.plot(Angles/(2*np.pi)*360,JNRJAnglesMaxes,'-ow',ms = 2,lw = 1.5)
    
    fa.tight_layout()
    
    fa.savefig(savepath + '\FxdA_'+ label + '_'+str(int(npts))+'npts_MaxCurve_JetNRJ.png')

    plt.close(fa)
    
    JVelAnglesMaxes = np.nanmax(JetVel,axis=(0,2))
    
    fa,axa = plt.subplots(dpi=200)
    axa.set_title('Maximum jet velocity per angle')
    axa.set_xlabel('Angle (°)')
    axa.set_ylabel('Max jet velocity [m/s]')
    axa.plot(Angles/(2*np.pi)*360,JVelAnglesMaxes,'-om',ms = 2,lw = 1.5)
    
    fa.tight_layout()
    
    fa.savefig(savepath + '\FxdA_'+ label + '_'+str(int(npts))+'npts_MaxCurve_JetVel.png')

    plt.close(fa)
    
    
    fa,axa = plt.subplots(dpi=200)
    axa.set_title('Maximum KNRJ, Dispertion, and jet fraction')
    axa.set_xlabel('Angle (°)')
    axa.set_ylabel('Common relative scale (%)')
    axa.plot(Angles/(2*np.pi)*360,100*JNRJAnglesMaxes/np.max(JNRJAnglesMaxes),'-ow',ms = 2,lw = 1.5,label = 'Max energy' )
    axa.plot(Angles/(2*np.pi)*360,100*DDAnglesMaxes/np.max(DDAnglesMaxes),'-or',ms = 2,lw = 1.5,label = 'Max dispertion')
    axa.plot(Angles/(2*np.pi)*360,100*JFAnglesMaxes/np.max(JFAnglesMaxes),'-ob',ms = 2,lw = 1.5,label = 'Max jetFrac')
    axa.plot(Angles/(2*np.pi)*360,100*SJFAnglesMaxes/np.max(SJFAnglesMaxes),'-og',ms = 2,lw = 1.5,label = 'Max SecjetFrac')
    axa.plot(Angles/(2*np.pi)*360,100*JVelAnglesMaxes/np.max(JVelAnglesMaxes),'-om',ms = 2,lw = 1.5, label = 'Max jetVel')
    plt.legend()
    
    fa.tight_layout()
    
    fa.savefig(savepath + '\FxdA_'+ label + '_'+str(int(npts))+'npts_MaxCurve_Together.png')

    plt.close(fa)
    
    ########### Fixed angle diagrams ###########       

    
    for a,ia in zip(Angles,range(len(Angles))):
        
        # Impact fraction in jet
        f0,ax0 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax0.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
        ax0.set_xlabel('Offcent' + NormStr)
        ax0.set_ylabel('DropDiam' + NormStr)
        
        ax0.set_xlim([0,np.max(RelOffCents)])
        ax0.set_ylim([0,np.max(RelDropDiams)])
    

        sc0 = ax0.scatter(2*meshOC[:,ia,:]/Normalisation,meshDD[:,ia,:]/Normalisation,c=JetFracs[:,ia,:],vmin=0, vmax = 100,
                          cmap='PuOr',marker='s',s=pointSize,zorder=2)

        cbar0 = plt.colorbar(sc0)
        cbar0.set_label('Impact fraction in the jet (%)')
        
        # ax0.set_aspect('equal')
        

        f0.tight_layout()
        
        f0.savefig(savepath + '\JetFrac\FixedAngle\FxdADgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_JetFrac.png')

        plt.close(f0)
            
        # Impact fraction in jet
        f0,ax0 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax0.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
        ax0.set_xlabel('Offcent' + NormStr)
        ax0.set_ylabel('DropDiam' + NormStr)
        
        ax0.set_xlim([0,np.max(RelOffCents)])
        ax0.set_ylim([0,np.max(RelDropDiams)])
    

        sc0 = ax0.scatter(2*meshOC[:,ia,:]/Normalisation,meshDD[:,ia,:]/Normalisation,c=SecJetFracs[:,ia,:],vmin=0, vmax = 100,
                          cmap='PuOr',marker='s',s=pointSize,zorder=2)

        cbar0 = plt.colorbar(sc0)
        cbar0.set_label('Impact fraction in the second jet (%)')
        
        # ax0.set_aspect('equal')
        

        f0.tight_layout()
        
        f0.savefig(savepath + '\SecJetFrac\FixedAngle\FxdADgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_SecJetFrac.png')

        plt.close(f0)
        
        
        # Impact fraction in jet
        f0,ax0 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax0.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
        ax0.set_xlabel('Offcent' + NormStr)
        ax0.set_ylabel('DropDiam' + NormStr)
        
        ax0.set_xlim([0,np.max(RelOffCents)])
        ax0.set_ylim([0,np.max(RelDropDiams)])
    
    
        sc0 = ax0.scatter(2*meshOC[:,ia,:]/Normalisation,meshDD[:,ia,:]/Normalisation,c=SecJetFracs[:,ia,:]+JetFracs[:,ia,:],vmin=0, vmax = 100,
                          cmap='PuOr',marker='s',s=pointSize,zorder=2)
    
        cbar0 = plt.colorbar(sc0)
        cbar0.set_label('Impact fraction in both jets (%)')
        
        # ax0.set_aspect('equal')
        
    
        f0.tight_layout()
        
        f0.savefig(savepath + '\BothJetFrac\FixedAngle\FxdADgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_BothJetFrac.png')
    
        plt.close(f0)
        

        
        
        print('Making and saving figures     ',end='\r')

        # Sheet/Jet volume ratio
        f1,ax1 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax1.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
        ax1.set_xlabel('Offcent' + NormStr)
        ax1.set_ylabel('DropDiam' + NormStr)
        
        ax1.set_xlim([0,np.max(RelOffCents)])
        ax1.set_ylim([0,np.max(RelDropDiams)])
        
        
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
        
        ax2.set_xlim([0,np.max(RelOffCents)])
        ax2.set_ylim([0,np.max(RelDropDiams)])
        
        sc2 = ax2.scatter(2*meshOC[:,ia,:]/Normalisation,meshDD[:,ia,:]/Normalisation,c=SheetWide[:,ia,:]*360/(2*np.pi),vmin=0, vmax = 360,
                          cmap='cividis',s=pointSize,marker='s')
        
        cbar2 = plt.colorbar(sc2)
        cbar2.set_label('Sheet opening [°]')
        # ax2.set_aspect('equal')
        f2.tight_layout()
        
        f2.savefig(savepath + '\SheetOpening\FixedAngle\FxdADgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_SheetOpen.png')

        plt.close(f2)
        
        
        # Shape factor
        f2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax2.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
        ax2.set_xlabel('Offcent' + NormStr)
        ax2.set_ylabel('DropDiam' + NormStr)
        
        ax2.set_xlim([0,np.max(RelOffCents)])
        ax2.set_ylim([0,np.max(RelDropDiams)])
        
        sc2 = ax2.scatter(2*meshOC[:,ia,:]/Normalisation,meshDD[:,ia,:]/Normalisation,c=ShapeFactor[:,ia,:],vmin=0,
                          cmap='viridis',s=pointSize,marker='s')
        
        cbar2 = plt.colorbar(sc2)
        cbar2.set_label('Shape factor')
        ax2.set_xlim([0,2])
        ax2.set_aspect('equal')
        f2.tight_layout()
        
        f2.savefig(savepath + '\ShapeFactor\FixedAngle\FxdADgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_ShapeFactor.png')

        plt.close(f2)
        
        
        print('Making and saving figures..     ',end='\r')
        
        # Dispertion distance
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
        ax3.set_xlabel('Offcent' + NormStr)
        ax3.set_ylabel('DropDiam' + NormStr)
        
        ax3.set_xlim([0,np.max(RelOffCents)])
        ax3.set_ylim([0,np.max(RelDropDiams)])
        
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
        
        ax3.set_xlim([0,np.max(RelOffCents)])
        ax3.set_ylim([0,np.max(RelDropDiams)])

        sc3 = ax3.scatter(2*meshOC[:,ia,:]/Normalisation,meshDD[:,ia,:]/Normalisation,c=DispertionDistVar[:,ia,:]/1000,cmap='plasma',
                          s=pointSize,marker='s')
        
        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Variability of dispertion dist [m]')
        # ax3.set_aspect('equal')
        f3.tight_layout()
        
        f3.savefig(savepath + '\DispertionDistVar\FixedAngle\FxdADgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_DispertionDistVar.png')

        

        plt.close(f3)
        
        
        print('Making and saving figures...     ',end='\r')
        
        # Kinetic energy
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
        ax3.set_xlabel('Offcent' + NormStr)
        ax3.set_ylabel('DropDiam' + NormStr)
        
        ax3.set_xlim([0,np.max(RelOffCents)])
        ax3.set_ylim([0,np.max(RelDropDiams)])

        sc3 = ax3.scatter(2*meshOC[:,ia,:]/Normalisation,meshDD[:,ia,:]/Normalisation,c=JetNRJ[:,ia,:],cmap='rainbow',
                          s=pointSize,marker='s')
        
        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Kinetic energy [mJ]')
        # ax3.set_aspect('equal')
        f3.tight_layout()
        
        f3.savefig(savepath + '\EnRJ\FixedAngle\FxdADgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_KineticEnergy.png')

        plt.close(f3)
        
        
        # Jet velocity [m/s]        
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('Cone angle = ' + str(round(a/(2*np.pi)*3600)/10))
        ax3.set_xlabel('Offcent' + NormStr)
        ax3.set_ylabel('DropDiam' + NormStr)
        
        ax3.set_xlim([0,np.max(RelOffCents)])
        ax3.set_ylim([0,np.max(RelDropDiams)])

        sc3 = ax3.scatter(2*meshOC[:,ia,:]/Normalisation,meshDD[:,ia,:]/Normalisation,c=JetVel[:,ia,:],cmap='viridis',vmin=0,vmax = 7.9,
                          s=pointSize,marker='s')
        
        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Jet velocity [m/s]')
        ax3.set_xlim([0,2])
        ax3.set_aspect('equal')
        f3.tight_layout()
        
        f3.savefig(savepath + '\JetVel\FixedAngle\FxdADgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(a/(2*np.pi)*3600)/10) + 'deg_JetVelocity.png')

        plt.close(f3)
        

    
    ########### Fixed drop diagrams ###########
    
    for rdd,idd in zip(RelDropDiams,range(len(RelDropDiams))):
        
        # Impact fraction in jet
        f0,ax0 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax0.set_title('DropDiam_Norm : ' + str(round(rdd*10)/10))
        ax0.set_xlabel('Offcent' + NormStr)
        ax0.set_ylabel('Cone angle [°]')
        ax0.set_xlim([0,np.max(RelOffCents)])
        

        sc0 = ax0.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=JetFracs[idd,:,:],vmin=0, vmax = 100,cmap='PuOr'
                          ,marker='s',s=pointSize)

        cbar0 = plt.colorbar(sc0)
        cbar0.set_label('Impact fraction in the jet (%)')
        f0.tight_layout()
        
        f0.savefig(savepath + '\JetFrac\FixedDrop\FxdDrDgm_'
            + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SR_JF.png')

        plt.close(f0)
        
        
        # Impact fraction in second jet
        f0,ax0 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax0.set_title('DropDiam_Norm : ' + str(round(rdd*10)/10))
        ax0.set_xlabel('Offcent' + NormStr)
        ax0.set_ylabel('Cone angle [°]')
        ax0.set_xlim([0,np.max(RelOffCents)])
        

        sc0 = ax0.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=SecJetFracs[idd,:,:],vmin=0, vmax = 100,cmap='PuOr'
                          ,marker='s',s=pointSize)

        cbar0 = plt.colorbar(sc0)
        cbar0.set_label('Impact fraction in the second jet (%)')
        f0.tight_layout()
        
        f0.savefig(savepath + '\SecJetFrac\FixedDrop\FxdDrDgm_'
            + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SR_SJF.png')

        plt.close(f0)
        
        
        print('Making and saving figures     ',end='\r')


        # Sheet/Jet volume ratio
        f1,ax1 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax1.set_title('DropDiam_Norm : ' + str(round(rdd*10)/10))
        ax1.set_xlabel('Offcent' + NormStr)
        ax1.set_ylabel('Cone angle [°]')
        ax1.set_xlim([0,np.max(RelOffCents)])
        
        sc1 = ax1.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=np.divide(np.subtract(100,JetFracs[idd,:,:]),JetFracs[idd,:,:]),
                          vmin=0, vmax = 3,cmap='plasma',s=pointSize,marker='s')
        
        cbar1 = plt.colorbar(sc1)
        cbar1.set_label('Sheet/Jet volume ratio')
        f1.tight_layout()
        
        f1.savefig(savepath + '\VolRatio\FixedDrop\FxdDrDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SR_VolRatio.png')
        
        plt.close(f1)
        
        
        print('Making and saving figures.     ',end='\r')
        

        # Sheet opening
        f2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax2.set_title('DropDiam_Norm : ' + str(round(rdd*10)/10))
        ax2.set_xlabel('Offcent' + NormStr)
        ax2.set_ylabel('Cone angle [°]')
        ax2.set_xlim([0,np.max(RelOffCents)])
        

        sc2 = ax2.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=SheetWide[idd,:,:]*360/(2*np.pi),vmin=0, 
                          vmax = 360,cmap='cividis',s=pointSize,marker='s')
        
        cbar2 = plt.colorbar(sc2)
        cbar2.set_label('Sheet opening [°]')
        f2.tight_layout()
        
        f2.savefig(savepath + '\SheetOpening\FixedDrop\FxdDrDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SR_SheetOpen.png')
        
        plt.close(f2)
        
        
        # Shape factor
        f2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax2.set_title('DropDiam_Norm : ' + str(round(rdd*10)/10))
        ax2.set_xlabel('Offcent' + NormStr)
        ax2.set_ylabel('Cone angle [°]')
        ax2.set_xlim([0,np.max(RelOffCents)])
        
        sc2 = ax2.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=ShapeFactor[idd,:,:],vmin=0,cmap='Blues_r',s=pointSize,marker='s')
        
        cbar2 = plt.colorbar(sc2)
        cbar2.set_label('Shape factor')
        f2.tight_layout()
        
        f2.savefig(savepath + '\ShapeFactor\FixedDrop\FxdDrDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SR_ShapeFactor.png')
        
        plt.close(f2)
        
        
        print('Making and saving figures..     ',end='\r')
        
        
        # dispersal distance
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('DropDiam_Norm : ' + str(round(rdd*10)/10))
        ax3.set_xlabel('Offcent' + NormStr)
        ax3.set_ylabel('Cone angle [°]')
        ax3.set_xlim([0,np.max(RelOffCents)])

        sc3 = ax3.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=DispertionDist[idd,:,:]/1000,vmin=0,cmap='inferno',s=pointSize,marker='s')
        
        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Average dispertion distance [m]')
        f3.tight_layout()
        
        f3.savefig(savepath + '\DispertionDist\FixedDrop\FxdDrDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SR_DispertionDist.png')
        
        plt.close(f3)
        
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('DropDiam_Norm : ' + str(round(rdd*10)/10))
        ax3.set_xlabel('Offcent' + NormStr)
        ax3.set_ylabel('Cone angle [°]')
        ax3.set_xlim([0,np.max(RelOffCents)])
       
        sc3 = ax3.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=DispertionDistVar[idd,:,:]/1000,vmin=0,
                          cmap='plasma',s=pointSize,marker='s')
        
        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Variability of dispertion dist [m]')
        f3.tight_layout()
        
        f3.savefig(savepath + '\DispertionDistVar\FixedDrop\FxdDrDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SR_DispertionDistVar.png')
        
        plt.close(f3)
        
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('DropDiam_Norm : ' + str(round(rdd*10)/10))
        ax3.set_xlabel('Offcent' + NormStr)
        ax3.set_ylabel('Cone angle [°]')
        ax3.set_xlim([0,np.max(RelOffCents)])
       
        sc3 = ax3.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=JetNRJ[idd,:,:],vmin=0,
                          cmap='rainbow',s=pointSize,marker='s')
        
        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Kinetic energy [mJ]')
        f3.tight_layout()
        
        f3.savefig(savepath + '\EnRJ\FixedDrop\FxdDrDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SR_NRJ.png')
        
        plt.close(f3)
        
        
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('DropDiam_Norm : ' + str(round(rdd*10)/10))
        ax3.set_xlabel('Offcent' + NormStr)
        ax3.set_ylabel('Cone angle [°]')
        ax3.set_xlim([0,np.max(RelOffCents)])
        
        sc3 = ax3.scatter(2*meshOC[idd,:,:]/Normalisation,meshA[idd,:,:]/(2*np.pi)*360,c=JetVel[idd,:,:],vmin=0,vmax = 7.9,
                          cmap='jet',s=pointSize,marker='s')
        
        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Jet velocity [m/s]')
        f3.tight_layout()
        
        f3.savefig(savepath + '\JetVel\FixedDrop\_1_FxdDrDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(rdd*10)/10) + 'SR_JetVel.png')
        
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
        
        sc0 = ax0.scatter(meshDD[:,:,ioc]/Normalisation,meshA[:,:,ioc]/(2*np.pi)*360,c=JetFracs[:,:,ioc],vmin=0, vmax = 100,cmap='PuOr',marker='s',s=pointSize)
    
        cbar0 = plt.colorbar(sc0)
        cbar0.set_label('Impact fraction in the jet (%)')
        f0.tight_layout()
        
        f0.savefig(savepath + '\JetFrac\FixedDist\FxdDiDgm_'
          + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OC_JetFrac.png')
  
        plt.close(f0)
        
        
    
        # Impact fraction in second jet
        f0,ax0 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax0.set_title('Offcent_Norm: ' + str(round(roc*10)/10))
        ax0.set_xlabel('DropDiam' + NormStr)
        ax0.set_ylabel('Cone angle [°]')
        ax0.set_xlim([0,np.max(RelDropDiams)])
        
        sc0 = ax0.scatter(meshDD[:,:,ioc]/Normalisation,meshA[:,:,ioc]/(2*np.pi)*360,c=SecJetFracs[:,:,ioc],vmin=0, vmax = 100,cmap='PuOr',marker='s',s=pointSize)
    
        cbar0 = plt.colorbar(sc0)
        cbar0.set_label('Impact fraction in the jet (%)')
        f0.tight_layout()
        
        f0.savefig(savepath + '\SecJetFrac\FixedDist\FxdDiDgm_'
          + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OC_SecJetFrac.png')
  
        plt.close(f0)
        
        
        print('Making and saving figures     ',end='\r')

        # Sheet/Jet volume ratio
        f1,ax1 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax1.set_title('Offcent_Norm: ' + str(round(roc*10)/10))
        ax1.set_xlabel('DropDiam' + NormStr)
        ax1.set_ylabel('Cone angle [°]')
        ax1.set_xlim([0,np.max(RelDropDiams)])
        
       
        sc1 = ax1.scatter(meshDD[:,:,ioc]/Normalisation,meshA[:,:,ioc]/(2*np.pi)*360,c=np.divide(np.subtract(100,JetFracs[:,:,ioc]),JetFracs[:,:,ioc]),
                          vmin=0, vmax = 3,cmap='plasma',marker='s',s=pointSize)
        
        cbar1 = plt.colorbar(sc1)
        cbar1.set_label('Sheet/Jet volume ratio')
        f1.tight_layout()
        
        f1.savefig(savepath + '\VolRatio\FixedDist\FxdDiDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OC_VolRatio.png')
        
        plt.close(f1)
        
        
        print('Making and saving figures.     ',end='\r')
        
        # Sheet/Jet volume ratio
        f2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax2.set_title('Offcent_Norm: ' + str(round(roc*10)/10))
        ax2.set_xlabel('DropDiam' + NormStr)
        ax2.set_ylabel('Cone angle [°]')
        ax2.set_xlim([0,np.max(RelDropDiams)])
        
        sc2 = ax2.scatter(meshDD[:,:,ioc]/Normalisation,meshA[:,:,ioc]/(2*np.pi)*360,c=SheetWide[:,:,ioc]*360/(2*np.pi),vmin=0, vmax = 360,
                          cmap='cividis',s=pointSize,marker='s')
        
        cbar2 = plt.colorbar(sc2)
        cbar2.set_label('Sheet opening [°]')
        f2.tight_layout()
        
        f2.savefig(savepath + '\SheetOpening\FixedDist\FxdDiDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OC_SheetOpen.png')
        
        plt.close(f2)
        
        # SheapeFactor
        f2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax2.set_title('Offcent_Norm: ' + str(round(roc*10)/10))
        ax2.set_xlabel('DropDiam' + NormStr)
        ax2.set_ylabel('Cone angle [°]')
        ax2.set_xlim([0,np.max(RelDropDiams)])
        
        sc2 = ax2.scatter(meshDD[:,:,ioc]/Normalisation,meshA[:,:,ioc]/(2*np.pi)*360,c=ShapeFactor[:,:,ioc],vmin=0,
                          cmap='Blues_r',s=pointSize,marker='s')
        
        cbar2 = plt.colorbar(sc2)
        cbar2.set_label('Shape Factor')
        f2.tight_layout()
        
        f2.savefig(savepath + '\ShapeFactor\FixedDist\FxdDiDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OC_ShapeFactor.png')
        
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
                    + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OC_DispertionDist.png')
        
        plt.close(f3)
        
        # kinetic energy ratio in the jet
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('Offcent_Norm: ' + str(round(roc*10)/10))
        ax3.set_xlabel('DropDiam' + NormStr)
        ax3.set_ylabel('Cone angle [°]')
        ax3.set_xlim([0,np.max(RelDropDiams)])
        
        sc3 = ax3.scatter(meshDD[:,:,ioc]/Normalisation,meshA[:,:,ioc]/(2*np.pi)*360,c=DispertionDistVar[:,:,ioc]/1000,
                          cmap='plasma',s=pointSize,marker='s')
        

        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Variability of dispertion dist [m]')
        f3.tight_layout()
        
        f3.savefig(savepath + '\DispertionDistVar\FixedDist\FxdDiDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OC_DispertionDistVar.png')
        
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
                    + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OC_NRJ.png')
        
        plt.close(f3)
        
        # Jet velocity [m/s] 
        f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
        ax3.set_title('Offcent_Norm: ' + str(round(roc*10)/10))
        ax3.set_xlabel('DropDiam' + NormStr)
        ax3.set_ylabel('Cone angle [°]')
        ax3.set_xlim([0,np.max(RelDropDiams)])
        
        sc3 = ax3.scatter(meshDD[:,:,ioc]/Normalisation,meshA[:,:,ioc]/(2*np.pi)*360,c=JetVel[:,:,ioc]/10000,vmin=0,vmax = 7.9,
                          cmap='jet',s=pointSize,marker='s')
        

        cbar3 = plt.colorbar(sc3)
        cbar3.set_label('Jet velocity [m/s]')
        f3.tight_layout()
        
        f3.savefig(savepath + '\JetVel\FixedDist\_1_FxdDiDgm_'
                    + label + '_'+str(int(npts))+'npts_' + str(round(roc*10)/10) + 'OC_JV.png')
        
        plt.close(f3)
        
    print('\nDone.')
    
    return
        
        



###
# 3. Cone optimization diagrams

def OptiDiagrams(ConeSurface,coneAngles,npts,ndrops,dropRadii,dropScaling,path,label):
    
    # ConeSurface = float, ConeAngles = np.linespace(x,y,z), dropScaling = np.linespace(a,b,c)
    
    if not len(coneAngles) == int(npts):
        raise ValueError('Invalid case ! npts should match the length of ''coneAngles'' and ''dropScaling'' ')
    
    oriType = 'Hgrad'
    velIni = 'Radial'
    meshType = 'zone'
    velType = 'full_div0'
    
    savepath = path + label + '_Res' + str(npts) + '_Drop' + str(ndrops)  
    
    MedianDropR = np.median(dropRadii)

    MedianDropS = 4*np.pi*MedianDropR**2
    
    DropNorm = ConeSurface/MedianDropS
    
    dropRadii = np.multiply(dropRadii,np.sqrt(DropNorm)) 
    
    MedianDropR = np.median(dropRadii)

    dropVels = np.sqrt(8/3*1000/1.3*10*dropRadii/500) # in m/s or mm/ms

    meshCA,meshDS = np.meshgrid(coneAngles,dropScaling)
    
    if os.path.exists(savepath+'\Data_impactVolumes.npy'):
        
        print('Loading previous simulations results...', end = '')
        
        
        impactVolumes = np.load(savepath + '\Data_impactVolumes.npy')
        jetVolumes = np.load(savepath + '\Data_jetVolumes.npy')
        SecjetVolumes = np.load(savepath + '\Data_SecjetVolumes.npy')
        efficiency = np.load(savepath + '\Data_efficiency.npy')
        jetNRJs = np.load(savepath + '\Data_jetNRJs.npy')
        jetVels = np.load(savepath + '\Data_jetVels.npy')
        DispersalDists = np.load(savepath + '\Data_DispersalDist.npy')
        DispersalHeights = np.load(savepath + '\Data_DispersalHeight.npy')
        DispersalVars = np.load(savepath + '\Data_DispersalVar.npy')
        TotalVolume = np.load(savepath + '\Data_TotalVolume.npy')

        print('Done !')
         
    else:
        os.makedirs(savepath,exist_ok=True) # create folder

        jetVolumes = np.empty(np.shape(meshCA))
        SecjetVolumes = np.empty(np.shape(meshCA))
        impactVolumes = np.empty(np.shape(meshCA))
        efficiency = np.empty(np.shape(meshCA))
        DispersalDists = np.empty(np.shape(meshCA))
        DispersalVars = np.empty(np.shape(meshCA))
        DispersalHeights = np.empty(np.shape(meshCA))
        jetNRJs = np.empty(np.shape(meshCA))
        jetVels = np.empty(np.shape(meshCA))

        

        # drop positions
        dropRs = np.sqrt(np.random.rand(ndrops)*(np.sqrt(ConeSurface*np.sin(np.max(coneAngles))/np.pi)+np.max(dropRadii)/2)**2)
        
    
        TotalVolume = np.sum(4/3*np.pi*dropRadii**3)/1000 # in [cm3] 
    
    
    
    
        ###### Simulations
    
        for a,ia in zip(coneAngles,range(len(coneAngles))):
            for dsc,ids in zip(dropScaling,range(len(dropScaling))):
    
                coneNum = ia*len(dropScaling) + ids + 1
                
                cr = np.sqrt(ConeSurface*np.sin(a)/np.pi)
    
                jetVolume = np.empty(np.shape(dropRadii))
                SecjetVolume = np.empty(np.shape(dropRadii))
                impactVolume = np.empty(np.shape(dropRadii))
                jetNRJ = np.empty(np.shape(dropRadii))
                jetVel = np.empty(np.shape(dropRadii))
                DispersalDist = np.empty(np.shape(dropRadii))
                DispersalHeight = np.empty(np.shape(dropRadii))
                
    
                for ds,dr,dv,di in zip(dropRadii,dropRs,dropVels,range(len(dropVels))):
    
                    print("Computing impacts for conditions " + str(coneNum) + "/" + str(len(dropScaling)*len(coneAngles)) 
                          + " : Impact n°"+str(di+1)+"/"+ str(len(dropRadii))+
                          # ". Params : DS = " + str(ds*dsc) + ", A = " + 
                          # str(a/np.pi*180) + ", DR =" + str(dr) + ", CR = " + str(cr) +
                          ".".ljust(10),end='\r')
    
                    if (dr<(cr+ds*dsc)*0.95) & (ds*dsc<(cr+dr)):
    
                        
                        I = dgc.Cone(cr,a).impact(dgc.Drop(ds*dsc,dr,72,dv),oriType,velIni,velType,meshType)
    
                        impactVolume[di] = I.VolFrac/100*4/3*np.pi*(ds/1000)**3
    
                        jetVolume[di] = impactVolume[di]*I.get_JetFrac()[0]/100
                                          
                        SecjetVolume[di] = impactVolume[di]*I.get_JetFrac()[1]/100
    
                        jetNRJ[di] = I.compute_JetNRJ()[0]
                        
                        jetVel[di] = I.compute_JetVel()
                        
                        DispersalDist[di] = I.compute_DispertionDist()[0]
                        
                        DispersalHeight[di] = I.compute_DispertionDist()[2]
                        
                        
                        del I
                        
                    else:
      
                        impactVolume[di] = 0
    
                        jetVolume[di] = 0
    
                        SecjetVolume[di] = 0
    
                        jetNRJ[di] = 0
                        
                        jetVel[di] = 0
                        
                        DispersalDist[di] = 0
    
                    
    
                impactVolumes[ids,ia] = np.sum(impactVolume) # in [m^3]
    
                jetVolumes[ids,ia] = np.sum(jetVolume) # in [m^3]
    
                SecjetVolumes[ids,ia] = np.sum(SecjetVolume) # in [m^3]
    
                efficiency[ids,ia] = np.sum(jetVolume)/np.sum(impactVolume)*100 # in %
    
                jetNRJs[ids,ia] = np.sum(jetNRJ) # [J]
    
                jetVels[ids,ia] = np.nanmedian(jetVel) # [J]
    
                DispersalDists[ids,ia] = np.median(DispersalDist[DispersalDist>0]) # [mm]
                
                DispersalHeights[ids,ia] = np.median(DispersalHeight[DispersalHeight>0]) # [mm]
                
                DispersalVars[ids,ia] =  np.diff(np.percentile(DispersalDist[DispersalDist>0],[25,75])) # [mm]
                
                
                ######### Plots of dispersal distance distributions
                
                # plt.hist(DispersalDist[DispersalDist>0],bins=30)
                # plt.plot(np.percentile(DispersalDist[DispersalDist>0],[25,25]),[0,30],'r--', lw = 1.5,label = '25-75% of data')
                # plt.plot(np.percentile(DispersalDist[DispersalDist>0],[75,75]),[0,30],'r--', lw =1.5,label = None)
                # plt.plot([np.median(DispersalDist[DispersalDist>0]),np.median(DispersalDist[DispersalDist>0])],[0,30],'r-',lw= 2,label='Median')
                # plt.title(str(coneNum) + "/" + str(len(dropScaling)*len(coneAngles)) 
                #       + " : Angle" + str(np.round(a/np.pi*180)) + '_ConeScale' + str(np.round(100/dsc**2)/100) + ' ' + str(sum(DispersalDist>0))+' impacts/'+str(len(DispersalDist)) + ' drops' )
                
                # plt.show()
                
                   
                
        np.save(savepath + '\Data_impactVolumes.npy',impactVolumes)
        np.save(savepath + '\Data_jetVolumes.npy',jetVolumes)
        np.save(savepath + '\Data_SecjetVolumes.npy',SecjetVolumes)
        np.save(savepath + '\Data_efficiency.npy',efficiency)
        np.save(savepath + '\Data_DispersalDist.npy',DispersalDists)
        np.save(savepath + '\Data_DispersalHeight.npy',DispersalHeights)
        np.save(savepath + '\Data_DispersalVar.npy',DispersalVars)
        np.save(savepath + '\Data_jetNRJs.npy',jetNRJs)
        np.save(savepath + '\Data_jetVels.npy',jetVels)
        np.save(savepath + '\Data_TotalVolume.npy',TotalVolume)

    # ############# Optimization diagrams

    pointSize = 250000/(npts**2) 

     # Impact volume
    f0,ax0 = plt.subplots(dpi=500,figsize = (7,6)) 
    ax0.set_title('Impact volume for random rain of ' + str(ndrops) + ' drops')
    ax0.set_xlabel('Cone angle (°)')
    ax0.set_ylabel('Cone area / MedianDrop area')

    sc0 = ax0.scatter(meshCA/(2*np.pi)*360,ConeSurface/(MedianDropR**2*np.pi*4)/meshDS**2,c=impactVolumes*1e6/TotalVolume*100,vmin=0,cmap='viridis',s=pointSize,marker='s')

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
    ax1.set_ylabel('Cone area / MedianDrop area')

    sc1 = ax1.scatter(meshCA/(2*np.pi)*360,ConeSurface/(MedianDropR**2*np.pi*4)/meshDS**2,c=jetVolumes*1e6/TotalVolume*100,vmin=0,cmap='PuOr',s=pointSize,marker='s')

    cbar1 = plt.colorbar(sc1)
    cbar1.set_label('Total volume ejected as jets (% of total rain volume)')
    f1.tight_layout()

    f1.savefig(savepath + '\OptiDgm_'
                + label + '_'+str(int(npts))+'npts_JetVolume.png')

    plt.close(f1)
    
     # Jet secondaire volume
    f1,ax1 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax1.set_title('Jet volume for random rain of ' + str(ndrops) + ' drops')
    ax1.set_xlabel('Cone angle (°)')
    ax1.set_ylabel('Cone area / MedianDrop area')

    sc1 = ax1.scatter(meshCA/(2*np.pi)*360,ConeSurface/(MedianDropR**2*np.pi*4)/meshDS**2,c=SecjetVolumes*1e6/TotalVolume*100,vmin=0,cmap='PuOr',s=pointSize,marker='s')

    cbar1 = plt.colorbar(sc1)
    cbar1.set_label('Total volume ejected as jets (% of total rain volume)')
    f1.tight_layout()

    f1.savefig(savepath + '\OptiDgm_'
                + label + '_'+str(int(npts))+'npts_SecJetVolume.png')

    plt.close(f1)

    # Jet/Impact volume
    f2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax2.set_title('Efficiency [jet/impact volumes] for random rain of ' + str(ndrops) + ' drops')
    ax2.set_xlabel('Cone angle (°)')
    ax2.set_ylabel('Cone area / MedianDrop area')

    sc2 = ax2.scatter(meshCA/(2*np.pi)*360,ConeSurface/(MedianDropR**2*np.pi*4)/meshDS**2,c=efficiency,vmin=0,cmap='cividis',s=pointSize,marker='s')
        

    cbar2 = plt.colorbar(sc2)
    cbar2.set_label('Efficiency [jet/impact volumes] (%)')
    f2.tight_layout()

    f2.savefig(savepath + '\OptiDgm_'
                + label + '_'+str(int(npts))+'npts_Efficiency.png')

    plt.close(f2)
    
    # Dispersal height
    f2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax2.set_title('Dispersal height for random rain of ' + str(ndrops) + ' drops')
    ax2.set_xlabel('Cone angle (°)')
    ax2.set_ylabel('Cone area / MedianDrop area')

    sc2 = ax2.scatter(meshCA/(2*np.pi)*360,ConeSurface/(MedianDropR**2*np.pi*4)/meshDS**2,c=DispersalHeights/100,vmin=0,cmap='plasma',s=pointSize,marker='s')
        

    cbar2 = plt.colorbar(sc2)
    cbar2.set_label('Average dispersal height [cm]')
    f2.tight_layout()

    f2.savefig(savepath + '\OptiDgm_'
                + label + '_'+str(int(npts))+'npts_DispersalHeight.png')

    plt.close(f2)
    
    # Dispersal dist
    f2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax2.set_title('Dispersal distance for random rain of ' + str(ndrops) + ' drops')
    ax2.set_xlabel('Cone angle (°)')
    ax2.set_ylabel('Cone area / MedianDrop area')

    sc2 = ax2.scatter(meshCA/(2*np.pi)*360,ConeSurface/(MedianDropR**2*np.pi*4)/meshDS**2,c=DispersalDists/1000,vmin=0,cmap='plasma',s=pointSize,marker='s')
        

    cbar2 = plt.colorbar(sc2)
    cbar2.set_label('Average dispersal distance [m]')
    f2.tight_layout()

    f2.savefig(savepath + '\OptiDgm_'
                + label + '_'+str(int(npts))+'npts_DispersalDist.png')

    plt.close(f2)
    
    # Dispersal dist variation
    f2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax2.set_title('Dispersal distance variability for random rain of ' + str(ndrops) + ' drops')
    ax2.set_xlabel('Cone angle (°)')
    ax2.set_ylabel('Cone area / MedianDrop area')

    sc2 = ax2.scatter(meshCA/(2*np.pi)*360,ConeSurface/(MedianDropR**2*np.pi*4)/meshDS**2,c=DispersalVars/1000,vmin=0,cmap='jet',s=pointSize,marker='s')
        

    cbar2 = plt.colorbar(sc2)
    cbar2.set_label('Variability between drops in dispersal distance [m]')
    f2.tight_layout()

    f2.savefig(savepath + '\OptiDgm_'
                + label + '_'+str(int(npts))+'npts_DispersalVar.png')

    plt.close(f2)

    DispersalVarRel = np.divide(DispersalVars,DispersalDists)*100
    
    f2,ax2 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax2.set_title('Dispersal distance relative variability for random rain of ' + str(ndrops) + ' drops')
    ax2.set_xlabel('Cone angle (°)')
    ax2.set_ylabel('Cone area / MedianDrop area')

    sc2 = ax2.scatter(meshCA/(2*np.pi)*360,ConeSurface/(MedianDropR**2*np.pi*4)/meshDS**2,c=DispersalVarRel,vmin=0,cmap='jet',s=pointSize,marker='s')
        

    cbar2 = plt.colorbar(sc2)
    cbar2.set_label('Variability between drops in dispersal distance [% of average]')
    f2.tight_layout()

    f2.savefig(savepath + '\OptiDgm_'
                + label + '_'+str(int(npts))+'npts_DispersalVarRel.png')

    plt.close(f2)

     # Kinetic energy
    f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax3.set_title('Kinetic energy in jets\n for random rain of ' + str(ndrops) + ' drops')
    ax3.set_xlabel('Cone angle (°)')
    ax3.set_ylabel('Cone area / MedianDrop area')

    sc3 = ax3.scatter(meshCA/(2*np.pi)*360,ConeSurface/(MedianDropR**2*np.pi*4)/meshDS**2,c=jetNRJs,cmap='rainbow',s=pointSize,marker='s')

    cbar3 = plt.colorbar(sc3)
    cbar3.set_label('Kinetic energy [J]')
    f3.tight_layout()

    f3.savefig(savepath + '\OptiDgm_'
                + label + '_'+str(int(npts))+'npts_KineticEnergy.png')

    plt.close(f3)

     # Jet velocity
    f3,ax3 = plt.subplots(dpi=150,figsize = (7,6)) 
    ax3.set_title('Jet velocity at ejection\n for random rain of ' + str(ndrops) + ' drops')
    ax3.set_xlabel('Cone angle (°)')
    ax3.set_ylabel('Cone area / MedianDrop area')

    sc3 = ax3.scatter(meshCA/(2*np.pi)*360,ConeSurface/(MedianDropR**2*np.pi*4)/meshDS**2,c=jetVels,cmap='rainbow',s=pointSize,marker='s')

    cbar3 = plt.colorbar(sc3)
    cbar3.set_label('Jet velocity [m/s]')
    f3.tight_layout()

    f3.savefig(savepath + '\OptiDgm_'
                + label + '_'+str(int(npts))+'npts_JetVelocity.png')

    plt.close(f3)
    
    
    
    
    print('\n\nFigures ploted and saved !')

