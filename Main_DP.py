# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 14:10:14 2023

@author: laplaud
"""

# plotting stuff
import matplotlib.pyplot as plt



# numbers handling
import numpy as np
import pandas as pd
    

# images handling
from skimage import io
from skimage.segmentation import active_contour

# to hide known warnings
import warnings
warnings.filterwarnings("ignore")

# General system functions
import os
import shutil
import sys

# my functions
import VallapFunc_DP as vf

print('Importations done.')


def DropCentering(P,StackList,Scale, **kwargs):     
    
    SD = pd.DataFrame(data=None,columns=['DropXc','DropYc','DropDiam','TargetXc','TargetYc','TargetDiam','OffCentering','OffCentering_pcDropDiam']) 
    
    #init and read kwargs
    DebugPlots = False
    
    for key, value in kwargs.items(): 
        if key == 'debug':
            DebugPlots = value  
        elif key == 'append':
            SD = value
        else:
            print('Unknown key : ' + key + '. Kwarg ignored.')
    
    for s in StackList:
        
        ImgD = io.imread(P + '\\Impacts\\' + s + '_Drop.tif') # get the  image from tiff stack
        ImgT = io.imread(P + '\\Impacts\\' + s + '_Target.tif') # get the  image from tiff stack
        
        TXs,TYs = vf.getContourPointsCoordinates(ImgT,'Select on the contour of the target (enter button to validate)')
        DXs,DYs = vf.getContourPointsCoordinates(ImgD,'Select on the contour of the drop (enter to validate)')
        
        dropCircle = vf.fitCircle(DXs,DYs)
        targetCircle = vf.fitCircle(TXs,TYs)
        
        dropXc = dropCircle[0]
        dropYc = dropCircle[1]
        dropR = dropCircle[2]
        
        targetXc = targetCircle[0]
        targetYc = targetCircle[1]
        targetR = targetCircle[2]
        
        OffC= np.sqrt(np.square(dropXc-targetXc)+np.square(dropYc-targetYc))
        OffCpc = OffC/dropR/2*100
        
        data = {'DropXc': dropXc,
            'DropYc': dropYc,
            'DropDiam': dropR*2/Scale,
            'TargetXc': targetXc,
            'TargetYc': targetYc,
            'TargetDiam': targetR*2/Scale,
            'OffCentering': OffC/Scale,
            'OffCentering_pcDropDiam': OffCpc} 

        SD = SD.append(pd.DataFrame(data=data,index = [s]))
        
        
            
        t = np.linspace(0,2*np.pi,1000)

        fig,ax = plt.subplots(dpi=300,facecolor='black')
        fig.suptitle('DropD : ' + str(round(dropR*2/Scale*10)/10) + ' mm. TargetD : ' +
                     str(round(targetR*2/Scale*10)/10) + ' mm. \nOffCent : ' + str(round(OffC/Scale*10)/10) + ' mm.' )
        plt.imshow(ImgD,cmap='gray')
        ax.plot(targetXc,targetYc,'.g',ms=2)
        ax.plot(targetXc+targetR*np.cos(t),targetYc+targetR*np.sin(t),'--g',lw=0.5)
        ax.plot(dropXc,dropYc,'.r',ms=2)
        ax.plot(dropXc+dropR*np.cos(t),dropYc+dropR*np.sin(t),'--r',lw=0.5)
        ax.plot([dropXc, targetXc],[dropYc,targetYc],'c',lw=0.2)
        ax.axis('off')      
        fig.savefig(P + '\\Impacts\\'+ s + '_Circles.png')
        
        if DebugPlots:
            plt.show()
        else:
            plt.close()
        
    return(SD)





def splashCarac(P,StackList,nimgs,SD,Scale, **kwargs):
    
    if not os.path.exists(P + '\\Splashes'):
            os.mkdir(P + '\\Splashes') # create global folder  
          
    
    #init and read kwargs
    DebugPlots = False
    
    for key, value in kwargs.items(): 
        if key == 'debug':
            DebugPlots = value                       
        else:
            print('Unknown key : ' + key + '. Kwarg ignored.')
    
    DD = pd.DataFrame(data=None)           
        
    for s in StackList:
        
        print('Analyzing ' + s + '...', end = '\r')
        
        if not os.path.exists(P + '\\Splashes\\' + s):
                os.mkdir(P + '\\Splashes\\' + s) # create global folder 
        
        SplashXcs = np.zeros(len(nimgs))
        SplashYcs = np.zeros(len(nimgs))        
        SplashRs = np.zeros(len(nimgs))
        
        for nimg,i in zip(nimgs,range(len(nimgs))): 
            
            print('Analyzing ' + s + '... Img ' + str(nimg), end = '\r')
        
            Img = io.imread(P + '\\' + s + '_Processed_Cut.tif', key=nimg) # get the  image from tiff stack

            DX = SD.loc[s,'DropXc']
            DY = SD.loc[s,'DropYc']
            DR = SD.loc[s,'DropDiam']/2*Scale

            TX = SD.loc[s,'TargetXc']
            TY = SD.loc[s,'TargetYc']
            TR = SD.loc[s,'TargetDiam']/2*Scale

            size = Img.shape

            ymax = size[1]
            xmax = size[0]

            npts = 50
            
            if (i == 0) | np.isnan(SplashRs[i-1]):
                side1 = np.array([np.linspace(0,xmax,npts), np.repeat([0],npts)])
                side2 = np.array([np.repeat([xmax],npts), np.linspace(0,ymax,npts)])
                side3 = np.array([np.linspace(xmax,0,npts), np.repeat([ymax],npts)])
                side4 = np.array([np.repeat([0],npts), np.linspace(ymax,0,npts)])

                init = np.concatenate((side1.T,side2.T, side3.T, side4.T))  
            else:
                t = np.linspace(0,2*np.pi,4*npts)
                Y = SplashXcs[i-1]
                X = SplashYcs[i-1]       
                R = SplashRs[i-1]*1.30*Scale
                init = np.array([X+R*np.cos(t), Y+R*np.sin(t)]).T
                
            
            snake = active_contour(Img,init, alpha=0.05, beta = 15, w_line=2, w_edge = 2)

            """
            alpha : float, optional
                Snake length shape parameter. Higher values makes snake contract faster.
    
            beta : float, optional
                Snake smoothness shape parameter. Higher values makes snake smoother.

            w_line : float, optional
                Controls attraction to brightness. Use negative values to attract toward dark regions.

            w_edge : float, optional
                Controls attraction to edges. Use negative values to repel snake from edges.
            """          
            SXs = snake[:,1]
            SYs = snake[:,0]
            
            
            ###### Faire de goodside un angle de 120+°
            GoodSide = np.sqrt(np.square(SXs-DX)+np.square(SYs-DY))>np.sqrt(np.square(SXs-TX)+np.square(SYs-TY))
            
            ##### Ajouté un tri pour ne prendre que les points assez proche du cercle précédents
            
            
            
            t = np.linspace(0,2*np.pi,200)

            fig, [ax1,ax2] = plt.subplots(ncols = 2,dpi=400)
            fig.suptitle('Circular fit on splash ' + s  +'_' + str(nimg))
            ax1.imshow(Img, cmap='gray')
            ax1.plot(init[:, 1], init[:, 0], '--r', lw=1)
            ax1.plot(snake[:, 1], snake[:, 0], '-b', lw=0.8)            
            ax1.axis([-5, Img.shape[1]+5, Img.shape[0]+5, -5])    
            ax1.set_xticks([])
            ax1.set_yticks([]) 
            ax1.set_title('Splash detection')
            ax2.imshow(Img,cmap='gray')
            ax2.plot(TX,TY,'.g',ms=2)
            ax2.plot(TX+TR*np.cos(t),TY+TR*np.sin(t),'--g',lw=0.5)
            ax2.plot(DX,DY,'.r',ms=2)
            ax2.plot(DX+DR*np.cos(t),DY+DR*np.sin(t),'--r',lw=0.5)
            ax2.plot(SXs[GoodSide],SYs[GoodSide],'*b',ms=0.1)
            ax2.axis([0, Img.shape[1], Img.shape[0], 0]) 
            ax2.set_xticks([])
            ax2.set_yticks([])
            
            if sum(GoodSide)>9:
                splashCircle = vf.fitCircle(SXs[GoodSide],SYs[GoodSide])  
                SplashXcs[i] = splashCircle[0]
                SplashYcs[i] = splashCircle[1]        
                SplashRs[i] = splashCircle[2]/Scale                
                
                ax2.set_title('Fit circle to selected points')
                ax2.plot(splashCircle[0],splashCircle[1],'*c',ms=2)
                ax2.plot(splashCircle[0]+splashCircle[2]*np.cos(t),splashCircle[1]+splashCircle[2]*np.sin(t),'-c',lw=0.8)
                

                if DebugPlots:
                    plt.show()
                else:
                    plt.close()
                    
            else:
                SplashXcs[i] = np.nan
                SplashYcs[i] = np.nan   
                SplashRs[i] = np.nan
                ax2.set_title('Not enough points to fit')
                
               
            fig.savefig(P + '\\Splashes\\'+ s + '\\' + str(nimg) + '.png')   
                
        data = {'SplashXs': SplashXcs,
                    'SplashYs': SplashYcs,
                    'SplashRs': SplashRs,
                    'SplashTargetDist': vf.dist(SplashXcs,SplashYcs,TX,TY)/Scale,
                    'SplashDropDist': vf.dist(SplashXcs,SplashYcs,DX,DY)/Scale,
                    'frames': nimgs
                   }

        DD = DD.append(pd.DataFrame(data=data,index = np.repeat(s,len(SplashXcs))))
        
        
        print('Analyzing ' + s + '... Done.'.ljust(10), end = '\n')
        
        OffC = SD.loc[s,'OffCentering_pcDropDiam']
        fig0, [ax01,ax02] = plt.subplots(ncols=2,dpi=300,facecolor='black')
        fig0.suptitle(s + ' - OffC = ' + str(round(OffC*10)/10) + '%')
        DD.loc[s].plot(x='frames',y='SplashRs',ax=ax01,style='bo-',legend=None)
        ax01.set_ylabel('Splash radius (mm)')    
        ax01.set_xlabel('Frames since impact')
        DD.loc[s].plot(x='frames',y='SplashTargetDist',ax=ax02,style='go-',legend=None)
        ax02.set_ylabel('Splash center distance from target (mm)')  
        ax02.set_xlabel('Frames since impact')     
        plt.tight_layout()
        fig0.savefig(P + '\\Splashes\\'+ s + '_SplashCarac.png')
        
    return(DD)
            