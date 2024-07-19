# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 09:11:18 2023

@author: laplaud
"""
import numpy as np
import numpy.matlib as ml
import matplotlib.pyplot as plt

from scipy.interpolate import LinearNDInterpolator, RegularGridInterpolator

import seaborn as sns

import DropGeometryFuncs as dgf

import VallapFunc_DP as vf

import time as time

from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'inline')

import os


##############################################################
# 1. A class for the cone where the drop is going to impact. #
##############################################################

class Cone:
    
    # Initialization. From user param (Radius, vertical slope angle).
    # Computes additional features using 'coneGeometry' 

    def __init__(self, radius, angle):

        self.Rcone = radius

        self.Alpha = angle
        
        self.Rcircle,self.Beta,self.Frac = dgf.coneGeometry(radius,angle)
        
        self.Xc = 0
        
        self.Yc = 0
        
        
        
    ###### Geometry related methods
    
    def Coordinates(self):
        
        return(self.Xc,self.Yc)
        
    def Parameters(self):
        
        return(self.Rcone,self.Alpha)
    
    def CircleParameters(self):
        
        return(self.Rcircle,self.Cone.Beta,self.Frac)
    
    def IsIn(self,R,Theta): 
        
        isinCone = R<self.Rcone
        
        return(isinCone)
    
    def IsInCircle(self,R,Theta): 
        
        isInCircle = R<self.Rcircle
        
        return(isInCircle)
    
    
    def Cone2Circle(self,X,Y):
        # (X,Y) points to transform, Alpha angle of the cone, Ad = 0 angle of removed sector bissecant (cone config)
    
        return(dgf.Cone2Circle(X,Y,self.Alpha,0))
    
    def Circle2Cone(self,X,Y): 
        # (X,Y) points to transform, Alpha angle of the cone, Ad = 0 angle of removed sector bissecant (circle config)
        
        return(dgf.Circle2Cone(X,Y,self.Alpha,0))
    
    
    
    ######  Impact with a drop
    
    def impact(self,drop,oriType,velIni,meshType):
        
        I = Impact(drop,self,oriType,velIni,meshType)
        
        return(I)
    
    
    ######  Representation method
    
    def draw(self, **kwargs):
        
        Drop = None
        
        Title     = 'kw: title= ''My title for the figure'''
        Xlabel_Ci = 'kw: xlabelCi= ''My xlabel for the circle config'''
        Ylabel_Ci = 'kw: ylabelCi= ''My ylabel for the circle config'''
        Xlabel_Co = 'kw: xlabelCo= ''My xlabel for the cone config'''
        Ylabel_Co = 'kw: ylabelCo= ''My ylabel for the cone config'''
        NoLabels  = False
        DropView = 'full'
        DropMesh = True
        
        ConeColor = 'g'
        ConeLW = 1
        
        for key, value in kwargs.items(): 
            if key == 'drop':
                Drop = value  
            elif key == 'title':
                Title = value
            elif key == 'xlabelCo':
                Xlabel_Co = value
            elif key == 'ylabelCo':
                Ylabel_Co = value
            elif key == 'xlabelCi':
                Xlabel_Ci = value
            elif key == 'ylabelCi':
                Ylabel_Ci = value
            elif key == 'nolabels':
                NoLabels = value
                if value == True:
                    Xlabel_Ci = ''
                    Ylabel_Ci = ''
                    Xlabel_Co = ''
                    Ylabel_Co = ''
            elif key == 'conecolor':
                ConeColor = value
            elif key == 'conelinewidth':
                ConeLW = value
            elif key == 'dropview':
                DropView = value
            elif key == 'dropmesh':
                DropMesh = value
                    
            else:
                print('Unknown key : ' + key + '. Kwarg ignored.')
        
        ### An abscisse vector
                
        tx = np.linspace(0,2*np.pi,100)
        
        ### Circle config : removed sector to form a cone 
                
        T1 = np.pi-self.Beta/2
        T2 = np.pi+self.Beta/2
        
        if T1>T2:
            Ts = np.linspace(np.mod(T1+np.pi,2*np.pi),np.mod(T2+np.pi,2*np.pi),20) - np.pi        
        else:                
            Ts = np.linspace(T1,T2,20)
        
        sectorT = np.append(np.append([T1],Ts),[T2])
        sectorR = np.append(np.append([0],self.Rcircle*np.ones(20)),[0])
        
        sectorX,sectorY = vf.ToCart(sectorT,sectorR,angle = 'rad')
        
        ### Adding the drop 
        
        if not Drop == None:
            
            Xdrop,Ydrop = Drop.Coordinates() # in cone config
                
            Rd = Drop.Parameters()[0]               
        
            ### Anaytical drop equation in circle config
            
            if Xdrop > Rd:
                ThetaBorder = np.sin(self.Alpha)*np.arcsin(Rd/Xdrop)
                Theta = np.linspace(-ThetaBorder,ThetaBorder,200)
                R1 = (Xdrop*np.cos(np.divide(Theta,np.sin(self.Alpha))) - np.sqrt(Rd**2-Xdrop**2*np.sin(np.divide(Theta,np.sin(self.Alpha)))**2))/(np.sin(self.Alpha))
                R2 = (Xdrop*np.cos(np.divide(Theta,np.sin(self.Alpha))) + np.sqrt(Rd**2-Xdrop**2*np.sin(np.divide(Theta,np.sin(self.Alpha)))**2))/(np.sin(self.Alpha))
                Theta = np.concatenate((Theta, np.flip(Theta)))
                R = np.concatenate((R1, R2))
            else:
                ThetaBorder = np.pi - self.Beta/2
                Theta = np.linspace(-ThetaBorder,ThetaBorder,200)
                R = (Xdrop*np.cos(np.divide(Theta,np.sin(self.Alpha))) + np.sqrt(Rd**2-Xdrop**2*np.sin(np.divide(Theta,np.sin(self.Alpha)))**2))/(np.sin(self.Alpha))
            
            
            X,Y = vf.ToCart(Theta,R, angle = 'rad')


        ### Plotting
        
        fig,ax = plt.subplots(dpi = 250, ncols = 2, figsize = (9,5))
        
        fig.suptitle(Title)
        ax[0].set_title('Circle config')
        ax[1].set_title('Cone config')
        
        ax[0].set_xlabel(Xlabel_Ci)
        ax[0].set_ylabel(Ylabel_Ci)
        
        ax[1].set_xlabel(Xlabel_Co)
        ax[1].set_ylabel(Ylabel_Co)
        
        if NoLabels:
            ax[0].set_xticks([])
            ax[0].set_yticks([])
            
            ax[1].set_xticks([])
            ax[1].set_yticks([])
        
        # (0,0) point
        ax[0].plot(0,0,'g*',ms = 3)
        ax[1].plot(0,0,'g*',ms = 3)
        

        ax[0].plot(self.Rcircle*np.cos(tx),self.Rcircle*np.sin(tx),color = ConeColor, lw=ConeLW,label = 'Cone',zorder=1);
        ax[0].plot(sectorX,sectorY,'--r',lw=1, label = 'Removed sector',zorder=2)

        ax[1].plot(self.Rcone*np.cos(tx),self.Rcone*np.sin(tx),color = ConeColor, lw=ConeLW, label = 'Cone',zorder=1);
        
        
        
        
        if not Drop == None:
            
            
            ax[0].plot(Drop.Xd/np.sin(self.Alpha),0,'b*',ms = 3)
            ax[1].plot(Drop.Xd,0,'b*',ms = 3)
            
            meshX,meshY,meshH = Drop.mesh()
            droplabel = 'Drop height'
            
            if DropView == 'impact':
                meshA,meshR = vf.ToCirc(meshX,meshY, angle='deg')
                inImpact = meshR<self.Rcone
                meshX,meshY,meshH = meshX[inImpact],meshY[inImpact],meshH[inImpact]
                droplabel = 'Drop height (impacting fraction)'
            
            cmap = plt.get_cmap('Blues')
            bluemap = vf.truncate_colormap(cmap,0.5,1,100)
            
            if DropMesh:
                ax[1].scatter(meshX,meshY,c=meshH,cmap=bluemap, s=15, zorder=3,label=droplabel)
            
            meshXci,meshYci = self.Cone2Circle(meshX, meshY)
            
            if DropMesh:
                ax[0].scatter(meshXci,meshYci,c=meshH,cmap=bluemap, s=15, zorder=3,label=droplabel)
            
            ax[0].plot(X,Y, 'b-', lw = 1,label='Deformed drop',zorder=4)      
            ax[1].plot(Xdrop + Rd*np.cos(tx), Ydrop + Rd*np.sin(tx), 'b-', lw = 1, label='Drop',zorder=4)
        
        
        # ax[0].set_xlim([0,500])
        # ax[0].set_ylim([-1.5,1.5])
        
        ax[0].set_aspect('equal')
        ax[1].set_aspect('equal')
        
        
        # ax[0].legend(fontsize='xx-small',loc='upper right')
        # ax[1].legend(fontsize='xx-small',loc='upper right')
        fig.tight_layout()
        
        return(fig,ax)
    
    
############################  
# 2. A class for the drop. #
############################


class Drop:
    
    
    def __init__(self, radius, offcent, npts, vel):

        self.Rdrop = radius
        
        self.Xd = offcent
        
        self.Yd = 0
        
        self.Vel = vel
        
        self.Npts = npts
        
        self.Volume = 4/3*np.pi*radius**3
        
        self.Mass = self.Volume*997e-9 # rho in [kg/mm3]
        
        ### Mesh of points in the drop
        
        
        ## Carthesian fixed mesh + borderpoints
        
        Xmin, Xmax, Ymin, Ymax = self.Xd-self.Rdrop,self.Xd+self.Rdrop,self.Yd-self.Rdrop,self.Yd+self.Rdrop
        
        Xs = np.linspace(Xmin, Xmax, npts)
        Ys = np.linspace(Ymin, Ymax, npts)
        
        Zs = np.linspace(-self.Rdrop,self.Rdrop,npts)

        
        Xgrid,Ygrid,Zgrid = np.meshgrid(Xs,Ys,Zs)
        

        XgridF = Xgrid.flatten()
        YgridF = Ygrid.flatten()
        ZgridF = Zgrid.flatten()
        
        XgridMid = XgridF[ZgridF ==0]
        YgridMid = YgridF[ZgridF ==0]

        Agrid,Rgrid = vf.ToCirc(XgridMid,YgridMid)

        isDrop = self.IsIn(Rgrid,Agrid)
        
        Rmax = self.Rdrop
        
        RgridBorder = np.ones((1,2*npts))*Rmax
        AgridBorder = np.linspace(0,2*np.pi,2*npts)
        
        XgridBorder,YgridBorder = vf.ToCart(AgridBorder,RgridBorder,angle='rad')
        
        XgridMidBord = np.append(XgridMid[isDrop],XgridBorder+self.Xd)
        YgridMidBord = np.append(YgridMid[isDrop],YgridBorder+self.Yd)
        

        self.meshXFull = Xgrid # in cone config
        self.meshYFull = Ygrid # in cone config
        self.meshZFull = Zgrid # in cone config
        
        self.meshX = XgridMidBord # in cone config
        self.meshY = YgridMidBord # in cone config
        self.meshH = dgf.SphereH(self.Rdrop,self.meshX,self.meshY,self.Xd)
        
        self.meshZmin = -self.meshH/2
        self.meshZmax = self.meshZmin + self.meshH
        
        meshH_z = dgf.SphereH(self.Rdrop,Xgrid,Ygrid,self.Xd)
        
        self.meshZmin3D = -meshH_z/2
        self.meshZmax3D = self.meshZmin3D + meshH_z
        
        
        
        inDrop3D = (ZgridF<self.meshZmax3D.flatten())&(ZgridF>self.meshZmin3D.flatten())
        
        self.meshX3D = XgridF[inDrop3D]
        self.meshY3D = YgridF[inDrop3D]
        self.meshZ3D = ZgridF[inDrop3D]
        
        ### H gradient 
        
        self.HgradInterp = []
        
        
    ###### Geometry related methods
    
    def Coordinates(self):
        
        return(self.Xd,self.Yd)
        
    def Parameters(self):
        
        return(self.Rdrop,self.Npts)
    
    def IsIn(self,R,Theta): 
        
        isInDrop = R**2 - 2*R*(self.Xd*np.cos(Theta)+self.Yd*np.sin(Theta)) + self.Xd**2 + self.Yd**2 <= self.Rdrop**2
        
        return(isInDrop)
        
    def mesh(self):
        
        return(self.meshX,self.meshY,self.meshH)    
    
    def mesh3D(self):
        
        return(self.meshX3D,self.meshY3D,self.meshZ3D)
    
    def Hgradient(self): 
        
        if self.HgradInterp == []:
        
            meshX,meshY,meshH = self.mesh() # Cone config
            
            Xmin, Xmax, Ymin, Ymax = self.Xd-self.Rdrop,self.Xd+self.Rdrop,self.Yd-self.Rdrop,self.Yd+self.Rdrop
            
            Xs = np.linspace(Xmin, Xmax, 200)
            Ys = np.linspace(Ymin, Ymax, 200)
            
            Xs = np.insert(Xs,0,Xmin*0.9)
            Xs = np.append(Xs,Xmax*1.1)
            Ys = np.insert(Ys,0,Ymin*0.9)
            Ys = np.append(Ys,Ymax*1.1)
            
            Xgrid,Ygrid = np.meshgrid(Xs,Ys)
            
            Hgrid = dgf.SphereH(self.Rdrop,Xgrid,Ygrid,self.Xd)
            
            Hgrad_y, Hgrad_x = np.gradient(Hgrid)
            
            HgradInterp_x = LinearNDInterpolator(list(zip(Xgrid.flatten(),Ygrid.flatten())),Hgrad_x.flatten())
            HgradInterp_y = LinearNDInterpolator(list(zip(Xgrid.flatten(),Ygrid.flatten())),Hgrad_y.flatten())
            
            self.HgradInterp = (HgradInterp_x,HgradInterp_y)
            
            return(HgradInterp_x,HgradInterp_y)
        
        else:
            
            return(self.HgradInterp[0],self.HgradInterp[1])
      

    ######  Impact with a cone
    
    def impact(self,cone,oriType,velIni,meshType):
        
        I = Impact(self,cone,oriType,velIni,meshType)
        
        return(I)

        
#################################################    
# 3. A class for an impact of a drop on a cone. #
################################################# 

class Impact:
    
    
    def __init__(self,drop,cone,oriType,velIni,meshType):
        
        self.Drop = drop
        
        self.meshType = meshType
        
        if meshType == 'point':
            if (oriType == 'Hgrad')|(oriType == 'Drop'):
                raise ValueError('Invalid case ! A single point mesh is incompatible with velocity origins "Drop" and "Hgrad"')
            if velIni == 'Radial':
                raise ValueError('Invalid case ! A single point mesh is incompatible with initial velocity type "Radial"')
                    

        
        self.Cone = cone
        
        self.velnorm = drop.Vel*np.sin(self.Cone.Alpha)
        
        self.veltan = drop.Vel*np.cos(self.Cone.Alpha)
        
        self.oriType = oriType
        self.velIni = velIni
        
        
        self.ori = []
        
        self.meshXci = []
        self.meshYci = []
        
        self.meshVXci_norm = []
        self.meshVYci_norm = []
        
        self.meshVXci_tan = []
        self.meshVYci_tan = []
        
        self.meshVXci_tan_div0 = []
        self.meshVYci_tan_div0 = []
        
        self.meshVXci = []
        self.meshVYci = []
        
        self.meshVXci_div0= []
        self.meshVYci_div0 = []

        self.SheetOpen = []
        self.wiXs = []
        self.wiYs = []

        self.trajX = []
        self.trajY = []
        self.TrajT = []
        
        self.JetFrac = []
        
        self.meshJFX = []
        self.meshJFY = []
        
        meshX,meshY,meshH = self.Drop.mesh() # Cone config
        
        inCone = np.sqrt(np.square(meshX) + np.square(meshY))<self.Cone.Rcone
        
        self.meshX = meshX[inCone]
        self.meshY = meshY[inCone]
        self.meshH = meshH[inCone]
        
        self.meshZmin = self.Drop.meshZmin[inCone]
        self.meshZmax = self.Drop.meshZmax[inCone]
        
        

        ##### Computation of contact height
        
        OffC = self.Drop.Xd 
        DropR = self.Drop.Rdrop
        ConeR = self.Cone.Rcone
        
        
        if (DropR < (ConeR - OffC)/np.cos(self.Cone.Alpha)):
            self.Zdrop = (DropR + OffC*np.cos(self.Cone.Alpha))/np.sin(self.Cone.Alpha)
        else :
            self.Zdrop = ConeR/np.tan(self.Cone.Alpha) + np.sqrt(DropR**2 - (OffC-ConeR)**2)

        
        
        
        ## Drop volume fraction in the impact
        
        self.VolFrac = dgf.volFrac([self.Drop.Xd],self.Drop.Rdrop,self.Cone.Rcone) 
                       
        
        
    
    ###############################################################################
    #                                                                             #
    #                          Simulation methods                                 #
    #                                                                             #
    ###############################################################################
    
    ## Trajectories orientation
    
    def compute_ori(self):
        
        if self.ori == []:
            
            
            
            meshXci,meshYci = self.Cone.Cone2Circle(self.meshX, self.meshY) 
            
            if self.meshType == 'point':
                nmesh = np.size(self.meshX)
                angles = np.linspace(0,np.pi*2,nmesh)
                R = self.Drop.Rdrop/100
                
                pointmeshX = R*np.cos(angles)
                pointmeshY = R*np.sin(angles)
            
            
            if self.oriType == 'Hgrad':
                
                           
                HgradInterp_x,HgradInterp_y = self.Drop.Hgradient()       
                
                meshOX = -HgradInterp_x(self.meshX,self.meshY) 
                meshOY = -HgradInterp_y(self.meshX,self.meshY) 
                
                # to circle config
                meshOXci,meshOYci = dgf.VelCone2Circle(meshOX, meshOY,self.meshX,self.meshY, self.Cone.Alpha) # Circle config, impacting fraction
    
                
                
                self.oriX = meshXci[np.argmax(self.meshH)]
                self.oriY = meshYci[np.argmax(self.meshH)]
                
                mesh1 = meshXci-meshXci[np.argmax(self.meshH)]
                mesh2 = meshYci-meshYci[np.argmax(self.meshH)]
                self.meshDist = np.sqrt(np.square(mesh1)+np.square(mesh2))
                
                self.impactR = np.max(self.meshDist)
                
                
            elif self.oriType == 'Drop':
                
                self.oriX,self.oriY = self.Cone.Cone2Circle(self.Drop.Xd, self.Drop.Yd)
                
                meshOXci = meshXci-self.oriX
                meshOYci = meshYci-self.oriY
                
                self.meshDist = np.sqrt(np.square(meshOXci)+np.square(meshOYci))
                self.impactR = np.max(self.meshDist)
                
                
                
                
            elif self.oriType == 'Hmax':
                
                if self.Drop.Xd>self.Cone.Rcone:
                    
                    self.oriX = self.Cone.Rcircle
                    
                else:
                    
                    self.oriX = self.Drop.Xd
                    
                self.oriY = 0
                
                if self.meshType == 'zone':
                    meshOXci = meshXci-self.oriX
                    meshOYci = meshYci-self.oriY
                
                    self.meshDist = np.sqrt(np.square(meshOXci)+np.square(meshOYci))
                    self.impactR = np.max(self.meshDist)
                
                elif self.meshType == 'point':
                    meshOXci = pointmeshX
                    meshOYci = pointmeshY
                    
                    meshXci = meshOXci + self.oriX
                    meshYci = meshOYci + self.oriY
            
                
            elif self.oriType == 'Central':
                
                self.oriX = meshXci.mean()
                self.oriY = meshYci.mean()
                
                if self.meshType == 'zone':
                    meshOXci = meshXci-self.oriX
                    meshOYci = meshYci-self.oriY
                
                    self.meshDist = np.sqrt(np.square(meshOXci)+np.square(meshOYci))
                    self.impactR = np.max(self.meshDist)
                
                elif self.meshType == 'point':
                    meshOXci = pointmeshX
                    meshOYci = pointmeshY
                    
                    meshXci = meshOXci + self.oriX
                    meshYci = meshOYci + self.oriY
                
            elif self.oriType == 'CentralHeight':
                
                self.oriX =  np.sum(np.multiply(meshXci,self.meshH))/np.sum(self.meshH)
                self.oriY = np.sum(np.multiply(meshYci,self.meshH))/np.sum(self.meshH)
                
                if self.meshType == 'zone':
                    meshOXci = meshXci-self.oriX
                    meshOYci = meshYci-self.oriY
                
                    self.meshDist = np.sqrt(np.square(meshOXci)+np.square(meshOYci))
                    self.impactR = np.max(self.meshDist)
                
                elif self.meshType == 'point':
                    meshOXci = pointmeshX
                    meshOYci = pointmeshY
                    
                    meshXci = meshOXci + self.oriX
                    meshYci = meshOYci + self.oriY
            
            self.ori = meshXci,meshYci,meshOXci,meshOYci
            
            self.meshXci = meshXci
            self.meshOXci = meshOXci
            self.meshYci = meshYci
            self.meshOYci = meshOYci
            
        return(self.ori)
    
        
    def orientation(self):
        
        ori = self.compute_ori()
        
        return(ori[0],ori[1],ori[2],ori[3])
    
    
    ## Trajectories speed vectors
    
    def compute_velocity_ini(self):
        
        if self.meshVXci_div0 == []:
    
            meshXci,meshYci,meshOXci,meshOYci = self.orientation() # Cone config, full drop       

            
            normCi = np.sqrt(np.square(self.meshOXci)+np.square(self.meshOYci))
                  
            if self.velIni == 'VelNorm':
                self.meshVXci_norm = np.divide(self.meshOXci,normCi)*self.velnorm
                self.meshVYci_norm = np.divide(self.meshOYci,normCi)*self.velnorm
                
            elif self.velIni == 'Radial':
                self.meshVXci_norm = np.multiply(np.divide(self.meshOXci,normCi),self.meshDist)/self.impactR*self.velnorm
                self.meshVYci_norm = np.multiply(np.divide(self.meshOYci,normCi),self.meshDist)/self.impactR*self.velnorm


            
            Tci,Rci = vf.ToCirc(meshXci,meshYci, angle = 'rad') 
            self.meshVXci_tan = -np.cos(Tci)*self.veltan
            self.meshVYci_tan = -np.sin(Tci)*self.veltan
            
            meshX,meshY = dgf.Circle2Cone(meshXci, meshYci, self.Cone.Alpha, 0)
            
            meshVX_tan,meshVY_tan = dgf.VelCircle2Cone(self.meshVXci_tan, self.meshVYci_tan, meshXci, meshYci, self.Cone.Alpha)
            
            meshVX_tan_div0 = np.ones(np.shape(meshVX_tan))*np.mean(meshVX_tan)
            meshVY_tan_div0 = np.ones(np.shape(meshVY_tan))*np.mean(meshVY_tan)
            
            
            self.meshVXci_tan_div0,self.meshVYci_tan_div0 = dgf.VelCone2Circle(meshVX_tan_div0,meshVY_tan_div0, meshX, meshY, self.Cone.Alpha)

            

            self.meshVXci = self.meshVXci_norm + self.meshVXci_tan
            self.meshVYci = self.meshVYci_norm + self.meshVYci_tan
            
            
            self.meshVXci_div0 = self.meshVXci_norm + self.meshVXci_tan_div0
            self.meshVYci_div0 = self.meshVYci_norm + self.meshVYci_tan_div0
            

    
    def velocity_ini(self,velType):
        
        self.compute_velocity_ini()
        
        if velType == 'norm':
            
            return(self.meshXci,self.meshYci,self.meshVXci_norm,self.meshVYci_norm)
        
        elif velType == 'tan':

            return(self.meshXci,self.meshYci,self.meshVXci_tan,self.meshVYci_tan)
        
        elif velType == 'tan_div0':

            return(self.meshXci,self.meshYci,self.meshVXci_tan_div0,self.meshVYci_tan_div0)
        
        elif velType == 'full':   

            return(self.meshXci,self.meshYci,self.meshVXci,self.meshVYci)
        
        elif velType == 'full_div0': 

            return(self.meshXci,self.meshYci,self.meshVXci_div0,self.meshVYci_div0)
        
        

    
    
    def compute_traj(self,Time):
       
        meshXci,meshYci,meshVXci,meshVYci = self.velocity_ini('full_div0')
        
        

        fieldnormCi = np.sqrt(np.square(meshVXci)+np.square(meshVYci))
            
        
               
        trajT = ml.repmat(Time,len(meshXci),1).T
        
        meshPts_X = ml.repmat(meshXci,len(Time),1)
        meshPts_Y = ml.repmat(meshYci,len(Time),1)
        
        meshVel_X = ml.repmat(meshVXci,len(Time),1)
        meshVel_Y = ml.repmat(meshVYci,len(Time),1)
        
                
        # meshPts_VX = ml.repmat(meshVXci,len(Time),1)
        # meshPts_VY = ml.repmat(meshVYci,len(Time),1)          
                
        # trajX = meshPts_X + np.multiply(meshPts_VX,trajT)
        # trajY = meshPts_Y + np.multiply(meshPts_VY,trajT)
                
        
        
        for it in range(1,len(Time)):
            
            meshPts_X[it,:] = meshPts_X[it-1,:] + meshVel_X[it-1,:]*(Time[it]-Time[it-1])
            meshPts_Y[it,:] = meshPts_Y[it-1,:] + meshVel_Y[it-1,:]*(Time[it]-Time[it-1])
            
            out_mask = np.sqrt(np.square(meshPts_X[it,:]) + np.square(meshPts_Y[it,:]))>self.Cone.Rcircle
            
            # out_mask = (np.sqrt(np.square(meshPts_X[it,:]) + np.square(meshPts_Y[it,:]))>self.Cone.Rcircle) & ()
            
            
            Xco,Yco = dgf.Circle2Cone(meshPts_X[it-1,:],meshPts_Y[it-1,:], self.Cone.Alpha, 0)
            Xco1,Yco1 = dgf.Circle2Cone(meshPts_X[it,:],meshPts_Y[it,:], self.Cone.Alpha, 0)
            
            
            tmp_velX,tmp_velY = dgf.VelCircle2Cone(meshVel_X[it-1,:],meshVel_Y[it-1,:],Xco,Yco, self.Cone.Alpha)
            
            
            velX_new,velY_new = dgf.VelCone2Circle(tmp_velX,tmp_velY,Xco1,Yco1, self.Cone.Alpha)
            
            chgmask = np.multiply(velX_new,meshVel_X[it-1,:]) < 0
            
            velX_new[chgmask] = meshVel_X[it-1,:][chgmask]
            velY_new[chgmask] = meshVel_Y[it-1,:][chgmask]
            
            
            meshVel_X[it,:][out_mask] = velX_new[out_mask]
            meshVel_Y[it,:][out_mask] = velY_new[out_mask]
            
            fieldnormCi = np.sqrt(np.square(meshVel_X[it,:])+np.square(meshVel_Y[it,:]))
            

            
        
        trajX = meshPts_X
        trajY = meshPts_Y

        self.trajVXci = meshVel_X
        self.trajVYci = meshVel_Y
       
        T,r =  vf.ToCirc(trajX,trajY,angle='rad')
        
        T = np.mod(T,2*np.pi)
        
        Beta = self.Cone.Beta # angle of sector to remove
        Xi,Yi = self.Cone.Cone2Circle(self.Drop.Xd, self.Drop.Yd) # Drop center in circle config
        OffA = vf.ToCirc(self.Cone.Xc-Xi,self.Cone.Yc-Yi,angle = 'rad')[0] #angle between impact and center
        T1 = np.mod(OffA - Beta/2,2*np.pi)
        XT1,YT1 = vf.ToCart(T1,self.Cone.Rcircle)
        T2 = np.mod(OffA + Beta/2,2*np.pi)
        XT2,YT2 = vf.ToCart(T2,self.Cone.Rcircle)
        if T1>T2:
            goodPts = ((T < T1) & (T > T2))            
        else:                
            goodPts = ((T < T1) | (T > T2))
        
        badPts = ~goodPts
        
        badX,badY = trajX[badPts],trajY[badPts]

        badT,badR = vf.ToCirc(badX-meshPts_X[badPts],badY-meshPts_Y[badPts],angle='rad')
        
        badIn = badR <= self.Cone.Rcircle
 
        
        # badOut = badR > self.Cone.Rcircle
        
        lx = -(meshPts_X[badPts]*np.tan(Beta/2)+np.abs(meshPts_Y[badPts]))/(np.abs(np.sin(badT))-np.abs(np.cos(badT))*np.tan(Beta/2))
        
        lz = (lx*np.abs(np.sin(badT))+np.abs(meshPts_Y[badPts]))/np.sin(Beta/2)

        # Sorting the bad point on both side of the removed sector
        closeT1 = np.abs((np.mod(badT,2*np.pi)-T1))<np.abs((np.mod(badT,2*np.pi)-T2))
        badT[closeT1] = T1
        badT[~closeT1] = T2
        
        badX[badIn],badY[badIn] = vf.ToCart(badT[badIn],badR[badIn]-lx[badIn]+lz[badIn])
        # badX[badOut],badY[badOut] = vf.ToCart(badT[badOut],badR[badOut]-lx[badOut]+lz[badOut])
        # 


        trajX[badPts] = badX
        trajY[badPts] = badY 
        
        self.trajX = trajX
        self.trajY = trajY
        self.trajT = trajT
                       
        
        return(self.trajX,self.trajY,self.trajT)
    
    
    
    
    ## Volume fraction in the jet (in % of impacting volume)
    
    def compute_JetFrac(self):
        
        meshXci,meshYci,meshVXci,meshVYci = self.velocity_ini('full_div0')
            
        # equations for the lines of the sector borders (y = c12*x)
        c1 = np.tan(np.pi-self.Cone.Beta/2)
        c2 = np.tan(-np.pi+self.Cone.Beta/2)
        
        # equation for the line along the trajectory (y = a*x + b)
        a = np.divide(meshVYci,meshVXci)
        b = np.divide((meshVXci*meshYci - meshVYci*meshXci),meshVXci)
        
        # intersection points xy coord
        xi1 = b/(c1-a)
        xi2 = b/(c2-a)
        
        yi1 = c1*xi1
        yi2 = c2*xi2
        
        dispX1 = xi1-meshXci
        dispX2 = xi2-meshXci
        
        velfactor1 = np.divide(dispX1,meshVXci)
        velfactor2 = np.divide(dispX2,meshVXci)
        
        inCone1 = np.sqrt(np.square(xi1)+np.square(yi1))<self.Cone.Rcircle
        inCone2 = np.sqrt(np.square(xi2)+np.square(yi2))<self.Cone.Rcircle
        
        
        # intersection 
        
        if self.Cone.Alpha>np.pi/6:     
            inter1 = ((velfactor1>0) & inCone1 & (xi1<0))
            inter2 = ((velfactor2>0) & inCone2 & (xi2<0))
            
            inter = inter1 | inter2
            
        elif  self.Cone.Alpha<np.pi/6:    
            inter1 = ((velfactor1>0)&(inCone1)&(xi1>0))
            inter2 = ((velfactor2>0)&(inCone2)&(xi2>0))
            
            inter = inter1 | inter2
            
        else:   
            inter1 = ((velfactor1>0)&(inCone1)) 
            inter2 = ((velfactor2>0)&(inCone2))
            
            inter = inter1 | inter2        
        
        
        cos1 = ( (np.square(meshXci-xi1) + np.square(meshYci-yi1)) + (np.square(xi1-meshXci) + np.square(yi1-c1*meshXci)) - (np.square(meshYci-c1*meshXci)))  / (2 * np.sqrt(np.square(meshXci-xi1) + np.square(meshYci-yi1)) * np.sqrt(np.square(xi1-meshXci) + np.square(yi1-c1*meshXci)))
        cos2 = ( (np.square(meshXci-xi2) + np.square(meshYci-yi2)) + (np.square(xi2-meshXci) + np.square(yi2-c2*meshXci)) - (np.square(meshYci-c2*meshXci)))  / (2 * np.sqrt(np.square(meshXci-xi2) + np.square(meshYci-yi2)) * np.sqrt(np.square(xi2-meshXci) + np.square(yi2-c2*meshXci)))
          
        
        self.meshJFxci = self.meshXci[inter]
        self.meshJFyci = self.meshYci[inter] 
        self.meshJFx = self.meshX[inter]
        self.meshJFy = self.meshY[inter]
        
        VXproj = self.meshVXci_div0.copy()
        VYproj = self.meshVYci_div0.copy()
        
        
        VXproj[inter1] = VXproj[inter1]*cos1[inter1]
        VYproj[inter1] = VYproj[inter1]*cos1[inter1]

        
        VXproj[inter2] = VXproj[inter2]*cos2[inter2]
        VYproj[inter2] = VYproj[inter2]*cos2[inter2]
        
        self.meshJFVx = VXproj[inter]
        self.meshJFVy = VYproj[inter]
        self.meshJH = self.meshH[inter] 
        
        
        JetFrac  = np.round(np.sum(inter)/np.size(inter)*1000)/10
        
        
        
        if (np.sqrt(np.square(self.Drop.Xd)+np.square(self.Drop.Yd))>(self.Drop.Rdrop+self.Cone.Rcone)) | ((self.Drop.Rdrop-np.sqrt(np.square(self.Drop.Xd)+np.square(self.Drop.Yd)))>(self.Cone.Rcone)):
            JetFrac = 0
            
        self.JetFrac= JetFrac
        
    
    def get_JetFrac(self):
        
        if self.JetFrac == []:
             self.compute_JetFrac()    
            
            
        
        return(self.JetFrac)
             
             
    
    # Sheet opening
    
    def compute_SheetOpening(self):
        
        if self.SheetOpen == []:
    
            trajX,trajY,trajT = self.compute_traj(np.linspace(0,5,50))
            
            tT,tR = vf.ToCirc(trajX,trajY, angle='rad')
            
            outCircle = tR>self.Cone.Rcircle
            
            trajX = trajX[outCircle]
            trajY = trajY[outCircle]
            
            trajXco,trajYco = self.Cone.Circle2Cone(trajX, trajY)
            
            
            tTco,tRco = vf.ToCirc(trajXco,trajYco, angle='rad')
            
            tTco[tTco<0] = tTco[tTco<0]+2*np.pi
            
            if len(tTco)>2:
                self.SheetOpen = (np.max(tTco)-np.min(tTco))
                self.wiXs = trajXco[(tTco==np.max(tTco))|(tTco==np.min(tTco))]
                self.wiYs = trajYco[(tTco==np.max(tTco))|(tTco==np.min(tTco))]
                
            else:
                self.SheetOpen = np.nan
                self.wiXs = np.nan
                self.wiYs = np.nan
                
    def SheetOpening(self):
        
        self.compute_SheetOpening()
        
        return(self.SheetOpen,self.wiXs,self.wiYs)
    
    
    def compute_ShapeFactor(self,Time,boxsize):
        
        if self.trajX == []:
            
            self.compute_traj(Time)
            
        trajXco, trajYco = self.Cone.Circle2Cone(self.trajX, self.trajY)
        
        trajVXco,trajVYco =dgf.VelCircle2Cone(self.trajVXci, self.trajVYci, trajXco, trajYco, self.Cone.Alpha)
        
                    
        jetmask = (trajXco < -self.Cone.Rcone) & (trajXco > -(1+boxsize)*self.Cone.Rcone) & (np.abs(trajYco) < self.Cone.Rcone*boxsize/2) 
        
        sidemask = (trajYco > self.Cone.Rcone) & (trajYco < (1+boxsize)*self.Cone.Rcone) & (np.abs(trajXco) < self.Cone.Rcone*boxsize/2) 
        
        
        # fig, ax = self.Cone.draw(conelinewidth=2, conecolor='g')
        
        # ax[1].plot(trajXco,trajYco,'ob',ms=1,zorder=0)
        # ax[1].plot(trajXco[jetmask],trajYco[jetmask],'or',ms=3,zorder=5)
        # ax[1].plot(trajXco[sidemask],trajYco[sidemask],'og',ms=3,zorder=5)
        
        # ax[1].set_aspect('equal')
        
        # plt.show()
        
        velsjetP = np.abs(trajVXco[jetmask])
        
        Tsjet = self.trajT[jetmask]
        
        Tsjoint = self.trajT[jetmask|sidemask]
        
        if np.sum(sidemask) == 0:
            
            velssideP = np.zeros(np.shape(trajVYco[jetmask]))
            Tsside = self.trajT[jetmask]
        else:
            velssideP = np.abs(trajVYco[sidemask])
            Tsside = self.trajT[sidemask]
            
                       
        velssidePmean = np.zeros(np.shape(np.unique(Tsjoint)))
        for i in np.unique(Tsside):
            tmp = np.where(Tsside == i)
            tmp2 = np.where(np.unique(Tsjoint) == i)
            velssidePmean[tmp2] = np.nanmean(velssideP[tmp])           
            
        velsjetPmean = np.zeros(np.shape(np.unique(Tsjoint)))    
        for i in np.unique(Tsjet):
            tmp = np.where(Tsjet == i)
            tmp2 = np.where(np.unique(Tsjoint) == i)
            velsjetPmean[tmp2] = np.nanmean(velsjetP[tmp])
            
            
        
        
        # f, ax = plt.subplots(dpi=200)
        # ax.plot(Tsside,velssideP,'or',ms=4,label='Side velocity // Y')
        # ax.plot(np.unique(Tsjoint),velssidePmean,'o-m',ms=4,label='Mean side velocity // Y')
        # ax.plot(Tsjet,velsjetP,'ob',ms=4,label = 'Jet velocity // X')
        # ax.plot(np.unique(Tsjoint),velsjetPmean,'o-c',ms=4,label = 'Mean jet velocity // X')
        # plt.legend()
        
        # ax.set_xlabel('Time (ms)')
        # ax.set_ylabel('Projected velocity (m/s)')
        
        # f.tight_layout()
        
        # f, ax = plt.subplots(dpi=200)
        # ax.plot(Tsside,velsside,'or',ms=4,label='Side velocity // Y')
        # ax.plot(Tsjet,velsjet,'ob',ms=4,label = 'Jet velocity // X')
        # ax.plot(Tsjet,velsjet,'ob',ms=4,label = 'Jet velocity // X')
        # plt.legend()
        
        # ax.set_xlabel('Time (ms)')
        # ax.set_ylabel('Velocity (m/s)')
        
        # f.tight_layout()
        
        OKmask = (velsjetPmean > 0) & (velssidePmean > 0)
        
        # Tok = np.unique(Tsjoint)[OKmask]
        
        # f, ax = plt.subplots(dpi=200)
        # ax.plot(Tok,np.divide(velssidePmean[OKmask],velsjetPmean[OKmask]),'og',ms=4,label='Side/Jet velocities')
        # plt.legend()
        
        # ax.set_xlabel('Time (ms)')
        # ax.set_ylabel('Relative projected velocity')
        
        # f.tight_layout()
        
         
        
        # veljet = np.mean(velsjet)
        
        # velside = np.mean(velsside)
        
        # ShapeFactor = velside/veljet

        if np.sum(OKmask) == 0:
            ShapeFactor = 0
        else:    
            ShapeFactor = np.nanmean(np.divide(velssidePmean[OKmask],velsjetPmean[OKmask]))
        
        return(ShapeFactor)
        
        
        
    
    def compute_JetNRJ(self):
        
        rho = 997e-9 # [kg/mm3]
        
        if np.size(self.meshJFVx)>0:
        
            Veq2 = np.average((np.square(self.meshJFVx)+np.square(self.meshJFVy)), weights=self.meshJH)
            # Weighted average of trajectories squared velocities by thickness
    
            self.JetNRJ = 4/3*np.pi*self.Drop.Rdrop**3*rho/2*Veq2*self.VolFrac/100*self.get_JetFrac()/100
            self.JetNRJ_Bal = self.JetNRJ*np.sin(2*(np.pi/2-self.Cone.Alpha))
        
        else:
            
            self.JetNRJ = 0
            self.JetNRJ_Bal = 0
        
        return(self.JetNRJ,self.JetNRJ_Bal)
    
    def compute_DispertionDist(self):
        
        g = 9.81*1e-3 # in [mm.ms-2]
        
        if np.size(self.meshJFVx)>0:
        
            Veq2 = np.average((np.square(self.meshJFVx)+np.square(self.meshJFVy)), weights=self.meshJH)
            
            self.DispertionDist = Veq2*np.sin(2*(np.pi/2-self.Cone.Alpha))/g
            
            self.DispertionDist_Var = np.sqrt(np.average(((np.square(self.meshJFVx)+np.square(self.meshJFVy))*np.sin(2*(np.pi/2-self.Cone.Alpha))/g-Veq2)**2, weights=self.meshJH))
      
        else:
            
            self.DispertionDist = 0
            
            self.DispertionDist_Var = 0
        
        
        return(self.DispertionDist,self.DispertionDist_Var)
    
    ###############################################################################
    #                                                                             #
    #                          Display methods                                    #
    #                                                                             #
    ###############################################################################
    
    def plot_splash_init(self,velType,**kwargs):
        
        Title     = 'kw: title= ''My title for the figure'''
        Xlabel_Ci = 'kw: xlabelCi= ''My xlabel for the circle config'''
        Ylabel_Ci = 'kw: ylabelCi= ''My ylabel for the circle config'''
        Xlabel_Co = 'kw: xlabelCo= ''My xlabel for the cone config'''
        Ylabel_Co = 'kw: ylabelCo= ''My ylabel for the cone config'''
        NoLabels  = False
        
        ConeColor = 'g'
        ConeLW = 1
        
        for key, value in kwargs.items(): 
            if key == 'title':
                Title = value
            elif key == 'xlabelCo':
                Xlabel_Co = value
            elif key == 'ylabelCo':
                Ylabel_Co = value
            elif key == 'xlabelCi':
                Xlabel_Ci = value
            elif key == 'ylabelCi':
                Ylabel_Ci = value
            elif key == 'nolabels':
                NoLabels = value
                if value == True:
                    Xlabel_Ci = ''
                    Ylabel_Ci = ''
                    Xlabel_Co = ''
                    Ylabel_Co = ''
            elif key == 'conecolor':
                ConeColor = value
            elif key == 'conelinewidth':
                ConeLW = value

                    
            else:
                print('Unknown key : ' + key + '. Kwarg ignored.')
        
       
        print('plotsplashinit')
        meshXci,meshYci,meshVXci,meshVYci = self.velocity_ini(velType)

        fieldnormCi = np.sqrt(np.square(meshVXci)+np.square(meshVYci))
 
        
        meshX, meshY = self.Cone.Circle2Cone(meshXci, meshYci)
  
        
        meshVX,meshVY = dgf.VelCircle2Cone(meshVXci,meshVYci, meshX, meshY, self.Cone.Alpha)
        
        fieldnorm = np.sqrt(np.square(meshVX)+np.square(meshVY))     
    
        
        fig, ax = self.Cone.draw(drop=self.Drop,nolabels=NoLabels,dropview = 'impact',conelinewidth=ConeLW,
                                 conecolor=ConeColor,title= Title + '_' + velType,xlabelCi=Xlabel_Ci,ylabelCi=Ylabel_Ci,xlabelCo=Xlabel_Co,ylabelCo=Ylabel_Co)
        
        
        
        q0 = ax[0].quiver(meshXci, meshYci, np.divide(meshVXci,fieldnormCi), np.divide(meshVYci,fieldnormCi),fieldnormCi,scale = 10,zorder=6,headlength=18,headaxislength=16)
        q1 = ax[1].quiver(meshX,   meshY,   np.divide(meshVX,fieldnorm),     np.divide(meshVY,fieldnorm),    fieldnorm,  scale = 15,zorder=6,headlength=18,headaxislength=16)

        ax[0].plot(self.oriX,self.oriY,'*m',ms = 10,zorder=7)
        ax[1].plot(self.oriX*np.sin(self.Cone.Alpha),self.oriY,'*m',ms = 10,zorder=7)

        if self.meshJFX == []:
            self.compute_JetFrac(velType)
        ax[0].scatter(self.meshJFxci,self.meshJFyci,c='r',zorder=5,s=7)
        ax[1].scatter(self.meshJFx,self.meshJFy,c='r',zorder=5,s=7)
        
        fig.colorbar(q0, ax = ax[0],orientation='horizontal',label = 'Velocity')
        fig.colorbar(q1, ax = ax[1],orientation='horizontal',label = 'Velocity')
    
    
    
    def plot_splash_traj(self,Time,**kwargs):
        
        Title     = 'kw: title= ''My title for the figure'''
        Xlabel_Ci = 'kw: xlabelCi= ''My xlabel for the circle config'''
        Ylabel_Ci = 'kw: ylabelCi= ''My ylabel for the circle config'''
        Xlabel_Co = 'kw: xlabelCo= ''My xlabel for the cone config'''
        Ylabel_Co = 'kw: ylabelCo= ''My ylabel for the cone config'''
        NoLabels  = False
        
        ConeColor = 'g'
        ConeLW = 1
        
        for key, value in kwargs.items(): 
            if key == 'title':
                Title = value
            elif key == 'xlabelCo':
                Xlabel_Co = value
            elif key == 'ylabelCo':
                Ylabel_Co = value
            elif key == 'xlabelCi':
                Xlabel_Ci = value
            elif key == 'ylabelCi':
                Ylabel_Ci = value
            elif key == 'nolabels':
                NoLabels = value
                if value == True:
                    Xlabel_Ci = ''
                    Ylabel_Ci = ''
                    Xlabel_Co = ''
                    Ylabel_Co = ''
            elif key == 'conecolor':
                ConeColor = value
            elif key == 'conelinewidth':
                ConeLW = value

                    
            else:
                print('Unknown key : ' + key + '. Kwarg ignored.')
                


        fig, ax = self.Cone.draw(drop=self.Drop,nolabels=NoLabels,dropview = 'impact',conelinewidth=ConeLW,dropmesh=False,
                                 conecolor=ConeColor,title=Title,xlabelCi=Xlabel_Ci,ylabelCi=Ylabel_Ci,xlabelCo=Xlabel_Co,ylabelCo=Ylabel_Co)
           

        trajX,trajY,trajT = self.compute_traj(Time)
            
        trajXco, trajYco = self.Cone.Circle2Cone(trajX, trajY)
    
        order = np.argsort(trajT.flatten())
        
        trajVci = np.sqrt(np.square(self.trajVXci)+ np.square(self.trajVYci))
        
        trajVXco,trajVYco =dgf.VelCircle2Cone(self.trajVXci, self.trajVYci, trajXco, trajYco, self.Cone.Alpha)
        
        
        trajVco = np.sqrt(np.square(trajVXco)+ np.square(trajVYco))
        

        sc0 = ax[0].scatter(trajX.flatten()[order],trajY.flatten()[order],c=trajVci.flatten()[order],cmap = 'Reds',s=1,zorder = -1)   
        ax[0].set_box_aspect(1)
        fig.colorbar(sc0, ax = ax[0],orientation='horizontal',label = 'Velocity')
        
        sc1 = ax[1].scatter(trajXco.flatten()[order],trajYco.flatten()[order],c=trajVco.flatten()[order],cmap = 'Reds',zorder = -1)
        
        
        ax[1].scatter(self.wiXs,self.wiYs,s=15,color='r',label='Sheet limits',zorder=4)

        fig.colorbar(sc1, ax = ax[1],orientation='horizontal',label = 'Velocity')
        
        
            
    
        return    
    
    
    def plot_3Dproj(self,title,resolution):
        
        R3D = np.sqrt(self.Drop.meshX3D**2 + self.Drop.meshY3D**2)
        inCone3D = R3D<self.Cone.Rcone
        
        MeshY3DinCone = self.Drop.meshY3D[inCone3D]
        MeshX3DinCone = self.Drop.meshX3D[inCone3D]
        MeshZ3DinCone = self.Drop.meshZ3D[inCone3D]+self.Zdrop
        
        MeshX2Dcircle,MeshY2Dcircle = self.Cone.Cone2Circle(MeshX3DinCone,MeshY3DinCone)
        
        MeshX3Dcircle,MeshY3Dcircle,MeshZ3Dcircle = dgf.Cone2CircleZ(MeshX3DinCone,MeshY3DinCone,MeshZ3DinCone,self.Cone.Alpha)
        
        
        ### Circle config : removed sector to form a cone 
                
        tx = np.linspace(0,2*np.pi,300)
        
        T1 = np.pi-self.Cone.Beta/2
        T2 = np.pi+self.Cone.Beta/2
        
        if T1>T2:
            Ts = np.linspace(np.mod(T1+np.pi,2*np.pi),np.mod(T2+np.pi,2*np.pi),20) - np.pi        
        else:                
            Ts = np.linspace(T1,T2,20)
        
        sectorT = np.append(np.append([T1],Ts),[T2])
        sectorR = np.append(np.append([0],self.Cone.Rcircle*np.ones(20)),[0])
        
        sectorX,sectorY = vf.ToCart(sectorT,sectorR,angle = 'rad')
        
        fig1,[[ax0,ax1],[ax2,ax3]] = plt.subplots(dpi=200,ncols=2,nrows=2)
        fig1.suptitle(title)
        ax0.set_title('Points')
        ax1.set_title('Hist')
        ax2.set_title('Interp')
        
        ax0.plot(self.Cone.Rcircle*np.cos(tx),self.Cone.Rcircle*np.sin(tx),color = 'g', lw=1,label = 'Circle')
        ax0.plot(sectorX,sectorY,'--m',lw=1, label = 'Removed sector')
        
        ax0.scatter(MeshX3Dcircle,MeshY3Dcircle,c='b',s=0.5,label='3Dproj',zorder=0)
        ax0.scatter(MeshX2Dcircle,MeshY2Dcircle,c='c',s=0.5,label='2Dproj',zorder=1)
        
        ax0.set_aspect('equal')
        xlim = ax0.get_xlim()
        ylim = ax0.get_ylim()
        
        
        
        ax1.plot(self.Cone.Rcircle*np.cos(tx),self.Cone.Rcircle*np.sin(tx),color = 'g', lw=1,label = 'Circle')
        ax1.plot(sectorX,sectorY,'--m',lw=1, label = 'Removed sector')

        ax2.plot(self.Cone.Rcircle*np.cos(tx),self.Cone.Rcircle*np.sin(tx),color = 'g', lw=1,label = 'Circle')
        ax2.plot(sectorX,sectorY,'--m',lw=1, label = 'Removed sector')
        
        ax3.plot(self.Cone.Rcircle*np.cos(tx),self.Cone.Rcircle*np.sin(tx),color = 'g', lw=1,label = 'Circle')
        ax3.plot(sectorX,sectorY,'--m',lw=1, label = 'Removed sector')
        
        
        
        I = ax1.hist2d(MeshX3Dcircle,MeshY3Dcircle,zorder=0,label='Drop height (density)',bins=resolution,cmap='plasma')

        ax1.set_xlim(xlim)
        ax1.set_ylim(ylim)
        
        ax1.set_aspect('equal')
        
        h,xedges,yedges = I[0],I[1],I[2]
        xvalues = (xedges[0:-1] + xedges[1:])/2
        yvalues = (yedges[0:-1] + yedges[1:])/2
        
        HeightInterp = RegularGridInterpolator((xvalues,yvalues),h,method='cubic')
        
        
        xx = np.linspace(np.min(xvalues),np.max(xvalues),50)
        yy = np.linspace(np.min(yvalues),np.max(yvalues),50)
        

        Xs, Ys = np.meshgrid(xx, yy, indexing='ij')

        # interpolator
        
        HeightGrid = HeightInterp((Xs, Ys))
        
        ax2.scatter(Xs, Ys,c= HeightGrid,cmap='plasma',marker='s')
        
        # sns.heatmap(HeightGrid,cmap = 'plasma',ax=ax2)


       
        ax2.set_xlim(xlim)
        ax2.set_ylim(ylim)
        
        ax2.set_aspect('equal')
        
        
        Hgrad_x, Hgrad_y = np.gradient(HeightGrid)
        
        

        
        xxx = np.linspace(np.min(xvalues),np.max(xvalues),15)
        yyy = np.linspace(np.min(yvalues),np.max(yvalues),15)
        
        GradXInterp = RegularGridInterpolator((xx,yy), Hgrad_x, method='cubic')
        GradYInterp = RegularGridInterpolator((xx,yy), Hgrad_y, method='cubic')
        
        

        XXs, YYs = np.meshgrid(xxx, yyy, indexing='ij')

        fieldnorm = np.sqrt(GradXInterp((XXs, YYs))**2+ GradYInterp((XXs, YYs))**2)
        
        ax3.quiver(XXs, YYs, np.divide(-GradXInterp((XXs, YYs)),fieldnorm), np.divide(-GradYInterp((XXs, YYs)),fieldnorm),fieldnorm,
                   scale=20,headlength=18,headaxislength=16,headwidth=8,cmap='plasma')

        
        ax3.set_xlim(xlim)
        ax3.set_ylim(ylim)
        
        ax3.set_aspect('equal')
        
        fig1.tight_layout()
        
        plt.show()
    
     
    def plot_3Dview(self,title):
        
        Hmax = 2*self.Drop.Rdrop + self.Cone.Rcone/np.tan(self.Cone.Alpha)
        
        Tmax = Hmax/self.Drop.Vel
        
        Times = np.linspace(0,4*Tmax,180)
        
        
        R3D = np.sqrt(self.Drop.meshX3D**2 + self.Drop.meshY3D**2)
        inCone3D = R3D<self.Cone.Rcone
        
        MeshY3DinCone = self.Drop.meshY3D[inCone3D]
        MeshX3DinCone = self.Drop.meshX3D[inCone3D]
         
        
        trajX,trajY,trajT = self.compute_traj(Times)
            
        trajXco, trajYco = self.Cone.Circle2Cone(trajX, trajY)
        
        trajZ = np.zeros(np.shape(trajXco))
        
        trajZco = np.sqrt(trajXco**2+trajYco**2)/np.tan(self.Cone.Alpha)
        
        
        timevect = np.unique(trajT)
        
        
        for tt,it in zip(Times,range(len(Times))):
            
            print(str(it+1) + '/' + str(len(Times)),end='\r')
        
            ZdropT = self.Zdrop - tt*self.Drop.Vel
            
            ZtotinCone = self.Drop.meshZ3D[inCone3D]+ZdropT
            
            
            ZtotUpinCone = self.Drop.meshZ3D[inCone3D]+ZdropT+2*self.Drop.Rdrop
            ZtotUpOutCone = self.Drop.meshZ3D[~inCone3D]+ZdropT+2*self.Drop.Rdrop
            
            
            
            Crossed = ZtotinCone < R3D[inCone3D]/np.tan(self.Cone.Alpha)
            
            ZtotinCone[Crossed] = np.nan
            
            
            
            CrossedUp = ZtotUpinCone < R3D[inCone3D]/np.tan(self.Cone.Alpha)
            
            ZtotUpinCone[CrossedUp] = np.nan
            
            
            
            CrossedUpOut = ZtotUpOutCone < 0
            
            ZtotUpOutCone[CrossedUpOut] = 0
            
            
            # get_ipython().run_line_magic('matplotlib', 'qt')
            fig = plt.figure(dpi=200)
            fig.suptitle(title)
            ax0 = fig.add_subplot(121, projection='3d')
            ax0.view_init(elev=25, azim=-100)
            ax0.set_title('Circle config')
            ax0.set_xlim(-4,4)
            ax0.set_ylim(-4,4)
            ax0.set_zlim(0,3)
            ax1 = fig.add_subplot(122, projection='3d')
            ax1.set_xlim(-2.5,2.5)
            ax1.set_ylim(-2.5,2.5)
            ax1.set_zlim(0,5)
            ax1.view_init(elev=30, azim=-100)
            ax1.set_title('Cone config')
            
            tx = np.linspace(0,2*np.pi,1000)
            
            
            ax1.plot(self.Cone.Rcone*np.cos(tx), self.Cone.Rcone*np.sin(tx),self.Cone.Rcone/np.tan(self.Cone.Alpha), 'g-', lw = 3, label='Cone',zorder=4)
            ax1.scatter(MeshX3DinCone,MeshY3DinCone,ZtotinCone,color='c',s=1)
            
            for t in tx:
                ax1.plot([0, self.Cone.Rcone*np.cos(t)],[0, self.Cone.Rcone*np.sin(t)],[0, self.Cone.Rcone/np.tan(self.Cone.Alpha)],'g-')
            
            
            ### Circle config : removed sector to form a cone 
                    
            T1 = np.pi-self.Cone.Beta/2
            T2 = np.pi+self.Cone.Beta/2
            
            if T1>T2:
                Ts = np.linspace(np.mod(T1+np.pi,2*np.pi),np.mod(T2+np.pi,2*np.pi),20) - np.pi        
            else:                
                Ts = np.linspace(T1,T2,20)
            
            sectorT = np.append(np.append([T1],Ts),[T2])
            sectorR = np.append(np.append([0],self.Cone.Rcircle*np.ones(20)),[0])
            
            sectorX,sectorY = vf.ToCart(sectorT,sectorR,angle = 'rad')
            
            ax0.plot(self.Cone.Rcircle*np.cos(tx),self.Cone.Rcircle*np.sin(tx),0,color = 'g', lw=3,label = 'Circle');
            ax0.plot(sectorX,sectorY,0,'--r',lw=3, label = 'Removed sector')
            Xi,Yi,Zi1 = dgf.Cone2CircleZ(MeshX3DinCone,MeshY3DinCone,ZtotinCone,self.Cone.Alpha)
            ax0.scatter(Xi,Yi,Zi1,color='c',s=1)
            
            
            H = tt*self.Drop.Vel
            
            points_mask = (trajT == 0)
            
            heightmask = np.zeros(np.shape(self.meshH), dtype=bool)
            heightmask = np.abs(self.Zdrop - self.meshH/2 - H) < np.sqrt(self.meshX**2+self.meshY**2)/np.tan(self.Cone.Alpha)
            
            ax0.scatter(trajX[points_mask][heightmask],trajY[points_mask][heightmask],trajZ[points_mask][heightmask],c='c',s = 1,vmin=0,vmax=self.Drop.Rdrop)
            
            ax1.scatter(trajXco[points_mask][heightmask],trajYco[points_mask][heightmask],trajZco[points_mask][heightmask],c='c',s = 1,vmin=0,vmax=self.Drop.Rdrop)
            
                
            for i in range(it):
                
                i= i+1
                
                Ti = timevect[i]
                

                points_mask = (trajT == timevect[it-i])
                
                
                
                Hi = Ti*self.Drop.Vel 
                
                
                
                heightmask = np.zeros(np.shape(self.meshH), dtype=bool)
                heightmask = np.abs(self.Zdrop - self.meshH/2 - Hi) < np.sqrt(self.meshX**2+self.meshY**2)/np.tan(self.Cone.Alpha)
                
                

                    
                
                
                ax0.scatter(trajX[points_mask][heightmask],trajY[points_mask][heightmask],trajZ[points_mask][heightmask],c='c',vmin=0,vmax=self.Drop.Rdrop,s = 1,zorder=10)
                ax1.scatter(trajXco[points_mask][heightmask],trajYco[points_mask][heightmask],trajZco[points_mask][heightmask],c='c',vmin=0,vmax=self.Drop.Rdrop,s = 1,zorder=10)
                
            
            ax0.set_aspect('equal')
            ax1.set_aspect('equal')
            
    
            
            
            savepath = r'd:\Users\laplaud\Desktop\PostDoc\Code\DropProject_WithAna\Figures\3D\Movies\\' + title
            
            os.makedirs(savepath,exist_ok = True)
            
            figname = savepath + '\Time_{:.2f}.png'
            
            
            fig.savefig(figname.format(tt))
            
            plt.close(fig)
            
            
            
    def rotating_3D(self,title):
        
        Hmax = 2*self.Drop.Rdrop + self.Cone.Rcone/np.tan(self.Cone.Alpha)
        
        Tmax = Hmax/self.Drop.Vel
        
        Times = np.linspace(0,4*Tmax,180)
        
        
        R3D = np.sqrt(self.Drop.meshX3D**2 + self.Drop.meshY3D**2)
        inCone3D = R3D<self.Cone.Rcone
        
        MeshY3DinCone = self.Drop.meshY3D[inCone3D]
        MeshX3DinCone = self.Drop.meshX3D[inCone3D]
        
        
        MeshY3DOutCone = self.Drop.meshY3D[~inCone3D]
        MeshX3DOutCone = self.Drop.meshX3D[~inCone3D]
        
        
        trajX,trajY,trajT = self.compute_traj(Times)
            
        trajXco, trajYco = self.Cone.Circle2Cone(trajX, trajY)
        
        trajZco = np.sqrt(trajXco**2+trajYco**2)/np.tan(self.Cone.Alpha)
        
        
        timevect = np.unique(trajT)
        
        
        for tt,it in zip(Times,range(len(Times))):
            
            print(str(it+1) + '/' + str(len(Times)),end='\r')
        
            ZdropT = self.Zdrop - tt*self.Drop.Vel
            
            ZtotinCone = self.Drop.meshZ3D[inCone3D]+ZdropT
            
            
            ZtotUpinCone = self.Drop.meshZ3D[inCone3D]+ZdropT+2*self.Drop.Rdrop
            ZtotUpOutCone = self.Drop.meshZ3D[~inCone3D]+ZdropT+2*self.Drop.Rdrop
            
            
            
            Crossed = ZtotinCone < R3D[inCone3D]/np.tan(self.Cone.Alpha)
            
            ZtotinCone[Crossed] = np.nan
            
            
            
            CrossedUp = ZtotUpinCone < R3D[inCone3D]/np.tan(self.Cone.Alpha)
            
            ZtotUpinCone[CrossedUp] = np.nan
            
            
            
            CrossedUpOut = ZtotUpOutCone < 0
            
            ZtotUpOutCone[CrossedUpOut] = 0
        
        ##### Rotating impact
        Az = 290-it*360/len(Times)
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.view_init(elev=18, azim=Az)
        
        # Make panes transparent
        ax.xaxis.pane.fill = False # Left pane
        ax.yaxis.pane.fill = False # Right pane
        
        # Remove grid lines
        ax.grid(False)
        
        # Remove tick labels
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        
        # Transparent spines
        ax.w_xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.w_yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        ax.w_zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
        
        # Transparent panes
        ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        
        # No ticks
        ax.set_xticks([]) 
        ax.set_yticks([]) 
        ax.set_zticks([])
        
        ax.set_xlim(-5,5)
        ax.set_ylim(-5,5)
        ax.set_zlim(0,5)
        
        
        tx1 = np.linspace(-1*np.pi/6,7*np.pi/6,1000)+np.pi/2+Az/180*np.pi
        tx2 = np.linspace(7*np.pi/6,11*np.pi/6,1000)+np.pi/2+Az/180*np.pi
        
        ax.plot(self.Cone.Rcone*np.cos(tx1), self.Cone.Rcone*np.sin(tx1),self.Cone.Rcone/np.tan(self.Cone.Alpha), 'g-', lw = 3, label='Cone',zorder=-2)
        
        ax.plot(self.Cone.Rcone*np.cos(tx2), self.Cone.Rcone*np.sin(tx2),self.Cone.Rcone/np.tan(self.Cone.Alpha), 'g-', lw = 3, label='Cone',zorder=201)
        
        ax.plot(self.Cone.Rcone*np.cos(0), self.Cone.Rcone*np.sin(0),self.Cone.Rcone/np.tan(self.Cone.Alpha), 'r.', ms = 3, label='Cone',zorder=202)
        
        ax.plot([0, self.Cone.Rcone*np.cos(0)],[0, self.Cone.Rcone*np.sin(0)],[0, self.Cone.Rcone/np.tan(self.Cone.Alpha)],
                '-r',zorder=199)
        
        for t in tx1[0:-1:2]:
            ax.plot([0, self.Cone.Rcone*np.cos(t)],[0, self.Cone.Rcone*np.sin(t)],[0, self.Cone.Rcone/np.tan(self.Cone.Alpha)],
                    '-',color = [0,0.7,0],zorder=-3)
        
        for t in tx2:
            ax.plot([0, self.Cone.Rcone*np.cos(t)],[0, self.Cone.Rcone*np.sin(t)],[0, self.Cone.Rcone/np.tan(self.Cone.Alpha)],
                    '-',color = [0,0.7,0],zorder=200)
          
        
        ax.plot([-3.5,3.5,3.5,-3.5],[3.5,3.5,-3.5,-3.5],[0,0,0,0],'.r',ms=5,zorder=-5)
        
        ax.scatter(MeshX3DinCone,MeshY3DinCone,ZtotUpinCone,color='b',s=1,zorder=0)
        
        ax.scatter(MeshX3DOutCone,MeshY3DOutCone,ZtotUpOutCone,color='b',s=1,zorder=0)
        
        
        H = tt*self.Drop.Vel
        
        points_mask = (trajT == 0)
        
        heightmask = np.zeros(np.shape(self.meshH), dtype=bool)
        heightmask = np.abs(self.Zdrop+2*self.Drop.Rdrop - self.meshH/2 - H) < np.sqrt(self.meshX**2+self.meshY**2)/np.tan(self.Cone.Alpha)
        
        ax.scatter(trajXco[points_mask][heightmask],trajYco[points_mask][heightmask],
                   trajZco[points_mask][heightmask],c=[0,0,0.9],s = 1,zorder=1)
        
            
        for i in range(it):
            
            i= i+1
            
            Ti = timevect[i]
            

            points_mask = (trajT == timevect[it-i])
            
            
            
            Hi = Ti*self.Drop.Vel 
            
            
            
            heightmask = np.zeros(np.shape(self.meshH), dtype=bool)
            heightmask = np.abs(self.Zdrop+2*self.Drop.Rdrop - self.meshH/2 - Hi) < np.sqrt(self.meshX**2+self.meshY**2)/np.tan(self.Cone.Alpha)
            

            ax.scatter(trajXco[points_mask][heightmask],trajYco[points_mask][heightmask],
                       trajZco[points_mask][heightmask],c=[0,0,0.9],s = 1,zorder=i)
            
        
        
        
        ax.set_aspect('equal')
        
        savepath = r'd:\Users\laplaud\Desktop\PostDoc\Code\DropProject_WithAna\Figures\3D\Movies\\RotatingView_' + title
        
        os.makedirs(savepath,exist_ok = True)
        
        figname = savepath + '\Right_Time_{:.2f}.png'
        
        
        fig.savefig(figname.format(tt))
        
        plt.close(fig)
        
        
    
    def plot_splash_movie(self,Time,label,xlims,ylims,**kwargs):
        
        Title     = 'kw: title= ''My title for the figure'''
        Xlabel_Ci = 'kw: xlabelCi= ''My xlabel for the circle config'''
        Ylabel_Ci = 'kw: ylabelCi= ''My ylabel for the circle config'''
        Xlabel_Co = 'kw: xlabelCo= ''My xlabel for the cone config'''
        Ylabel_Co = 'kw: ylabelCo= ''My ylabel for the cone config'''
        NoLabels  = False
        
        ConeColor = 'g'
        ConeLW = 1
        
        for key, value in kwargs.items(): 
            if key == 'title':
                Title = value
            elif key == 'xlabelCo':
                Xlabel_Co = value
            elif key == 'ylabelCo':
                Ylabel_Co = value
            elif key == 'xlabelCi':
                Xlabel_Ci = value
            elif key == 'ylabelCi':
                Ylabel_Ci = value
            elif key == 'nolabels':
                NoLabels = value
                if value == True:
                    Xlabel_Ci = ''
                    Ylabel_Ci = ''
                    Xlabel_Co = ''
                    Ylabel_Co = ''
            elif key == 'conecolor':
                ConeColor = value
            elif key == 'conelinewidth':
                ConeLW = value

                    
            else:
                print('Unknown key : ' + key + '. Kwarg ignored.')
                
        savepath = r'd:\Users\laplaud\Desktop\PostDoc\Code\DropProject_WithAna\Figures\Movies\\' + label 
         
        os.makedirs(savepath,exist_ok=True)
        
        trajX,trajY,trajT = self.compute_traj(Time)
            
        trajXco, trajYco = self.Cone.Circle2Cone(trajX, trajY)
        
        
        
        # movie
        
        timevect = np.unique(trajT)
        
        
        for T,iT in zip(timevect,range(len(timevect))):
            
            fig, ax = self.Cone.draw(drop=self.Drop,nolabels=NoLabels,dropview = 'impact',conelinewidth=ConeLW,dropmesh=False,
                                     conecolor=ConeColor,title=Title,xlabelCi=Xlabel_Ci,ylabelCi=Ylabel_Ci,xlabelCo=Xlabel_Co,ylabelCo=Ylabel_Co)
            
            
            
            H = T*self.Drop.Vel
            
            points_mask = (trajT == 0)
            
            contactDist = np.ones(np.shape(self.meshH))*self.Drop.Rdrop + ((self.Cone.Rcone - 
                                                                           np.sqrt(np.square(self.meshX)+np.square(self.meshY)))
                                                                           /np.tan(self.Cone.Alpha))
            
            heightmask = np.zeros(np.shape(self.meshH), dtype=bool)
            heightmask = self.meshH/2 > np.abs(H - contactDist)
            
            ax[0].scatter(trajX[points_mask][heightmask],trajY[points_mask][heightmask],c=self.meshH[heightmask]/2,cmap='Blues',s = 1,vmin=0,vmax=self.Drop.Rdrop)
            
            ax[1].scatter(trajXco[points_mask][heightmask],trajYco[points_mask][heightmask],c=self.meshH[heightmask]/2,cmap='Blues',s = 1,vmin=0,vmax=self.Drop.Rdrop)
            
                
            for i in range(iT):
                
                i= i+1
                
                Ti = timevect[i]
                

                points_mask = (trajT == timevect[iT-i])
                
                
                
                Hi = Ti*self.Drop.Vel 
                
                
                
                heightmask = np.zeros(np.shape(self.meshH), dtype=bool)
                heightmask = self.meshH/2 > np.abs(Hi - contactDist)
                

                    
                
                
                ax[0].scatter(trajX[points_mask][heightmask],trajY[points_mask][heightmask],c=self.meshH[heightmask]/2,cmap='Blues',vmin=0,vmax=self.Drop.Rdrop,s = 1,zorder=10)
                ax[1].scatter(trajXco[points_mask][heightmask],trajYco[points_mask][heightmask],c=self.meshH[heightmask]/2,cmap='Blues',vmin=0,vmax=self.Drop.Rdrop,s = 1,zorder=10)
                
            ax[1].set_xlim(xlims)
            ax[1].set_ylim(ylims)
                
            figname = savepath + '\Time_' + str(np.round(T*1000)) + '.png'
            
            plt.savefig(figname)
            
            plt.close(fig)
           
    
        return    
    
    
    
    
    