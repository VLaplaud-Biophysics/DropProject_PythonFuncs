# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 09:11:18 2023

@author: laplaud
"""
import numpy as np
import numpy.matlib as ml
import matplotlib.pyplot as plt

from scipy.interpolate import LinearNDInterpolator

# import seaborn as sns

import DropGeometryFuncs as dgf

import VallapFunc_Sim as vf

# import time as time

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
    
        return(dgf.Cone2Circle(X,Y,self.Alpha))
    
    def Circle2Cone(self,X,Y): 
        # (X,Y) points to transform, Alpha angle of the cone, Ad = 0 angle of removed sector bissecant (circle config)
        
        return(dgf.Circle2Cone(X,Y,self.Alpha))
    
    
    
    ######  Impact with a drop
    
    def impact(self,drop,oriType,velIni,velType,meshType):
        
        I = Impact(drop,self,oriType,velIni,velType,meshType)
        
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
        ax[0].set_title('Circle (unfolded) config')
        ax[1].set_title('Cone (original) config')
        
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
        
        
        ax[0].set_aspect('equal')
        ax[1].set_aspect('equal')
        
        fig.tight_layout()
        
        return(fig,ax)
    
    
############################  
# 2. A class for the drop. #
############################


class Drop:
    
    
    def __init__(self, radius, offcent, npts, vel):

        self.Rdrop = radius
        
        self.Xd = offcent
               
        self.Vel = vel
        
        self.Npts = npts
        
        self.Volume = 4/3*np.pi*radius**3
        
        self.Mass = self.Volume*997e-9 # rho in [kg/mm3]
        
        ### Mesh of points in the drop
        
        
        ## Carthesian fixed mesh + borderpoints
        
        Xmin, Xmax, Ymin, Ymax = self.Xd-self.Rdrop,self.Xd+self.Rdrop,-self.Rdrop,self.Rdrop
        
        Xs = np.linspace(Xmin, Xmax, npts)
        Ys = np.linspace(Ymin, Ymax, npts)
        
        Zs = np.linspace(-self.Rdrop,self.Rdrop,npts)

        
        Xgrid,Ygrid,Zgrid = np.meshgrid(Xs,Ys,Zs)
        
        Xgrid2D,Ygrid2D = np.meshgrid(Xs,Ys)
        

        XgridF = Xgrid.flatten()
        YgridF = Ygrid.flatten()
        ZgridF = Zgrid.flatten()
        
        
        XgridMid = Xgrid2D.flatten()
        YgridMid = Ygrid2D.flatten()

        Agrid,Rgrid = vf.ToCirc(XgridMid,YgridMid)

        isDrop = np.sqrt(np.square(XgridMid-self.Xd)+np.square(YgridMid)) < self.Rdrop
        
        Rmax = self.Rdrop*0.999
        
        RgridBorder = np.ones((1,2*npts))*Rmax
        AgridBorder = np.linspace(0,2*np.pi,2*npts)
        
        XgridBorder,YgridBorder = vf.ToCart(AgridBorder,RgridBorder,angle='rad')
        
        XgridMidBord = np.append(XgridMid[isDrop],XgridBorder+self.Xd)
        YgridMidBord = np.append(YgridMid[isDrop],YgridBorder)
        

        self.meshXFull = Xgrid # in cone config
        self.meshYFull = Ygrid # in cone config
        self.meshZFull = Zgrid # in cone config
        
        self.meshX =  XgridMidBord # in cone config
        self.meshY =  YgridMidBord # in cone config
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
        
        return(self.Xd,0)
        
    def Parameters(self):
        
        return(self.Rdrop,self.Npts)
    
    def IsIn(self,R,Theta): 
        
        isInDrop = R**2 - 2*R*self.Xd*np.cos(Theta)+ self.Xd**2  <= self.Rdrop**2
        
        return(isInDrop)
        
    def mesh(self):
        
        return(self.meshX,self.meshY,self.meshH)    
    
    def mesh3D(self):
        
        return(self.meshX3D,self.meshY3D,self.meshZ3D)
    
    def Hgradient(self): 
        
        if self.HgradInterp == []:
        
            meshX,meshY,meshH = self.mesh() # Cone config
            
            Xmin, Xmax, Ymin, Ymax = self.Xd-self.Rdrop,self.Xd+self.Rdrop,-self.Rdrop,self.Rdrop
            
            Xs = np.linspace(Xmin*0.9, Xmax*1.1, 200)
            Ys = np.linspace(Ymin*0.9, Ymax*1.1, 200)
        
            
            Xgrid,Ygrid = np.meshgrid(Xs,Ys)
            
            Hgrid = dgf.SphereH(self.Rdrop,Xgrid,Ygrid,self.Xd)
            
            Hgrad_y, Hgrad_x = np.gradient(Hgrid,np.abs(Xs[1]-Xs[0]),np.abs(Ys[1]-Ys[0]))
            
            HgradInterp_x = LinearNDInterpolator(list(zip(Xgrid.flatten(),Ygrid.flatten())),Hgrad_x.flatten())
            HgradInterp_y = LinearNDInterpolator(list(zip(Xgrid.flatten(),Ygrid.flatten())),Hgrad_y.flatten())
            
            self.HgradInterp = (HgradInterp_x,HgradInterp_y)
            
            return(HgradInterp_x,HgradInterp_y)
        
        else:
            
            return(self.HgradInterp[0],self.HgradInterp[1])
      

    ######  Impact with a cone
    
    def impact(self,cone,oriType,velIni,velType,meshType):
        
        I = Impact(self,cone,oriType,velIni,velType,meshType)
        
        return(I)

        
#################################################    
# 3. A class for an impact of a drop on a cone. #
################################################# 

class Impact:
    
    
    def __init__(self,drop,cone,oriType,velIni,velType,meshType):
        
        self.MISSED = False
        
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
        self.velType = velType
        
        
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
        self.trajT = []
        
        self.JetFrac = []
        
        self.meshJFX = []
        self.meshJFY = []
        
        meshX,meshY,meshH = self.Drop.mesh() # Cone config
        
        
        
        inCone = np.sqrt(np.square(meshX) + np.square(meshY))<=self.Cone.Rcone
        
        
        
        self.meshX = meshX[inCone]
        self.meshY = meshY[inCone]
        self.meshH = meshH[inCone]
        
        self.meshZmin = self.Drop.meshZmin[inCone]
        self.meshZmax = self.Drop.meshZmax[inCone]
        
        if self.meshH.size == 0 :
            self.MISSED = True

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
                       

        
        ## simulate trajectories
        if not self.MISSED:
            self.compute_velocity_ini()
            Time = 20 # ms
            self.compute_traj(np.linspace(0,Time,Time*50))
    
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
                self.oriY = 0
                
                mesh1 = meshXci-meshXci[np.argmax(self.meshH)]
                mesh2 = meshYci-0
                self.meshDist = np.sqrt(np.square(mesh1)+np.square(mesh2))
                
                self.impactR = np.max(self.meshDist)
                
                
            elif self.oriType == 'Drop':
                
                self.oriX,self.oriY = self.Cone.Cone2Circle(self.Drop.Xd, 0)
                
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
            
            meshX,meshY = dgf.Circle2Cone(meshXci, meshYci, self.Cone.Alpha)
            
            meshVX_tan,meshVY_tan = dgf.VelCircle2Cone(self.meshVXci_tan, self.meshVYci_tan, meshXci, meshYci, self.Cone.Alpha)
            
            meshVX_tan_div0 = np.ones(np.shape(meshVX_tan))*np.mean(meshVX_tan)
            meshVY_tan_div0 = np.ones(np.shape(meshVY_tan))*np.mean(meshVY_tan)
            
            
            self.meshVXci_tan_div0,self.meshVYci_tan_div0 = dgf.VelCone2Circle(meshVX_tan_div0,meshVY_tan_div0, meshX, meshY, self.Cone.Alpha)

            

            self.meshVXci = self.meshVXci_norm + self.meshVXci_tan
            self.meshVYci = self.meshVYci_norm + self.meshVYci_tan
            
            
            self.meshVXci_div0 = self.meshVXci_norm + self.meshVXci_tan_div0
            self.meshVYci_div0 = self.meshVYci_norm + self.meshVYci_tan_div0
            

    
    def velocity_ini(self,velType):
        
        
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
       
        meshXci,meshYci,meshVXci,meshVYci = self.velocity_ini(self.velType)
        
               
        trajT = ml.repmat(Time,len(meshXci),1).T
        
        meshPts_Xci = ml.repmat(meshXci,len(Time),1)
        meshPts_Yci = ml.repmat(meshYci,len(Time),1)
        
        meshVel_Xci = ml.repmat(meshVXci,len(Time),1)
        meshVel_Yci = ml.repmat(meshVYci,len(Time),1)
             
        
        for it in range(1,len(Time)):
            
            # position at it is position at it-1 + (velocity at it-1 * dt)
            
            meshPts_Xci[it,:] = meshPts_Xci[it-1,:] + meshVel_Xci[it-1,:]*(Time[it]-Time[it-1])
            meshPts_Yci[it,:] = meshPts_Yci[it-1,:] + meshVel_Yci[it-1,:]*(Time[it]-Time[it-1])
            
            

            #######################################################################
            ###### Finding points that crossed the removed sector 
            
            # trajectory points on  ith image
            
            
            trajXci1 = meshPts_Xci[it,:].copy()
            trajYci1 = meshPts_Yci[it,:].copy()

            # expression in radial coordinates
           
            T,r =  vf.ToCirc(trajXci1,trajYci1,angle='rad')
            
            T = np.mod(T,2*np.pi)
            
            # Finding points that crossed the sector border
            
            Beta = self.Cone.Beta # angle of sector to remove
            
            ### Circle config : removed sector to form a cone 
                    
            T1 = np.pi-Beta/2
            T2 = np.pi+Beta/2
            
            
            badPts1 = (T > T1) & (T<np.pi)
            badPts2 = (T < T2) & (T>np.pi)
            badPts = badPts1 | badPts2

            
            if not all(badPts==False):
                
                sectorX,sectorY = vf.ToCart(T1,self.Cone.Rcircle,angle = 'rad')
                sectorX2,sectorY2 = vf.ToCart(T2,self.Cone.Rcircle,angle = 'rad')
                
                # plt.figure(dpi=300)
                # plt.plot(0,0,'*y')
                # plt.quiver(0,0,sectorX,sectorY,color='gold',angles='xy')
                # plt.quiver(0,0,sectorX2,sectorY2,color='gold',angles='xy')
                # plt.plot([0,sectorX],[0,sectorY],'y-')
                # plt.plot([0,sectorX2],[0,sectorY2],'y-')
                # plt.scatter(meshPts_Xci[it,:][badPts], meshPts_Yci[it,:][badPts],s=17,zorder=1, color = 'r')
                # plt.scatter(meshPts_Xci[it,:][~badPts], meshPts_Yci[it,:][~badPts],s=10,zorder=1, color = 'g')
            
                # plt.quiver(meshPts_Xci[it,:][badPts], meshPts_Yci[it,:][badPts],meshVel_Xci[it-1,:][badPts],meshVel_Yci[it-1,:][badPts],color='r',angles='xy')
                
                
                 
                if not all(badPts1==False):
                    
                    # position correction in circle config            
                    
                    meshPts_Xci[it,:][badPts1], meshPts_Yci[it,:][badPts1] = vf.ToCart(T1, r[badPts1], angle='rad') 
                    
                    # plt.scatter(meshPts_Xci[it,:][badPts1], meshPts_Yci[it,:][badPts1],s=13,zorder=2,color = 'c')
                
                    
                    # Points that have crossed, corrected expressed in cone configuration
                    
                    trajXco1, trajYco1 = dgf.Circle2Cone(meshPts_Xci[it,:][badPts1], meshPts_Yci[it,:][badPts1], self.Cone.Alpha) 
                    
               
                    a = [[x,y] for x,y in zip(meshVel_Xci[it-1,:][badPts1],meshVel_Yci[it-1,:][badPts1])]
                    b = [sectorX,sectorY]
                    # velocity correction
                    
                    
                    
                    velX_new1 = np.sign(np.dot(a,b))*sectorX/np.sqrt(np.square(sectorX)+np.square(sectorY))*np.sqrt(np.square(meshVel_Xci[it-1,:][badPts1])+np.square(meshVel_Yci[it-1,:][badPts1]))
                    velY_new1 = np.sign(np.dot(a,b))*sectorY/np.sqrt(np.square(sectorX)+np.square(sectorY))*np.sqrt(np.square(meshVel_Xci[it-1,:][badPts1])+np.square(meshVel_Yci[it-1,:][badPts1]))
                    
                    # print(np.square(meshVel_Xci[it-1,:][badPts1])+np.square(meshVel_Yci[it-1,:][badPts1]))
                    # print(np.shape(meshVel_Yci[it-1,:][badPts1]))
                    # print(np.shape(velX_new1))
                    # print(np.shape(velY_new1))
                    # print(sectorX)
                    velX_new1[(np.sign(np.dot(a,b))<0)&(trajXco1<0)] = 0
                    velY_new1[(np.sign(np.dot(a,b))<0)&(trajXco1<0)] = 0
                    
                    # plt.quiver(meshPts_Xci[it,:][badPts1], meshPts_Yci[it,:][badPts1],velX_new1,velY_new1,color='c',angles='xy')
                    
                    
                    
                    meshVel_Xci[it-1,:][badPts1] = velX_new1.copy()
                    meshVel_Yci[it-1,:][badPts1] = velY_new1.copy()
                    
                    meshVel_Xci[it,:][badPts1] = velX_new1.copy()
                    meshVel_Yci[it,:][badPts1] = velY_new1.copy()
                    
                
                if not all(badPts2==False):
                    
                    
                    # position correction in circle config 
                    meshPts_Xci[it,:][badPts2], meshPts_Yci[it,:][badPts2] = vf.ToCart(T1, r[badPts2], angle='rad') 
                
                
                    
                    # plt.scatter(meshPts_Xci[it,:][badPts2], meshPts_Yci[it,:][badPts2],s=13,zorder=2,color = 'c')
                
                    
                    # Points that have crossed, corrected expressed in cone configuration
                    
                    trajXco12, trajYco12 = dgf.Circle2Cone(meshPts_Xci[it,:][badPts2], meshPts_Yci[it,:][badPts2], self.Cone.Alpha) 
                    
                    a2 = [[x,y] for x,y in zip(meshVel_Xci[it-1,:][badPts2],meshVel_Yci[it-1,:][badPts2])]
                    b2 = [sectorX2,sectorY2]
                    # velocity correction
                    
                    velX_new2 = np.sign(np.dot(a2,b2))*sectorX2/np.sqrt(np.square(sectorX2)+np.square(sectorY2))*np.sqrt(np.square(meshVel_Xci[it-1,:][badPts2])+np.square(meshVel_Yci[it-1,:][badPts2]))
                    velY_new2 = np.sign(np.dot(a2,b2))*sectorY2/np.sqrt(np.square(sectorX2)+np.square(sectorY2))*np.sqrt(np.square(meshVel_Xci[it-1,:][badPts2])+np.square(meshVel_Yci[it-1,:][badPts2]))
                    
                    
                    
                    velX_new2[(np.sign(np.dot(a2,b2))<0)&(trajXco12<0)] = 0
                    velY_new2[(np.sign(np.dot(a2,b2))<0)&(trajXco12<0)] = 0
                    
                    # plt.quiver(meshPts_Xci[it,:][badPts2], meshPts_Yci[it,:][badPts2],velX_new2,velY_new2,color='c',angles='xy')
                    
                    # plt.show()
                    
                    meshVel_Xci[it-1,:][badPts2] = velX_new2.copy()
                    meshVel_Yci[it-1,:][badPts2] = velY_new2.copy()
                    
                    meshVel_Xci[it,:][badPts2] = velX_new2.copy()
                    meshVel_Yci[it,:][badPts2] = velY_new2.copy()
                    
                ############################################################
                ##### Correcting points trajectory outside of the cone
                
                # Point whose ith position is outside of the cone
                
                out_mask = np.sqrt(np.square(meshPts_Xci[it,:]) + np.square(meshPts_Yci[it,:]))>self.Cone.Rcircle
                
                           
                # Point position in cone configuration for it-1 and it
                
                Xco,Yco = dgf.Circle2Cone(meshPts_Xci[it-1,:],meshPts_Yci[it-1,:], self.Cone.Alpha) 
                Xco1,Yco1 = dgf.Circle2Cone(meshPts_Xci[it,:],meshPts_Yci[it,:], self.Cone.Alpha)
                
                # velocity of the points in the cone configuration
             
                tmp_velX,tmp_velY = dgf.VelCircle2Cone(meshVel_Xci[it-1,:],meshVel_Yci[it-1,:],Xco,Yco, self.Cone.Alpha)
                
                
                # transfert in circle config
                
                velX_new,velY_new = dgf.VelCone2Circle(tmp_velX,tmp_velY,Xco1,Yco1, self.Cone.Alpha)
                
                
                meshVel_Xci[it,:][out_mask] = velX_new[out_mask]
                meshVel_Yci[it,:][out_mask] = velY_new[out_mask]

            
            
            
            
            
            

        
        self.trajVXci = meshVel_Xci
        self.trajVYci = meshVel_Yci
        
        self.trajX = meshPts_Xci
        self.trajY = meshPts_Yci
        self.trajT = trajT
                       
        return
    
    def get_traj(self):
        
        return(self.trajX,self.trajY,self.trajT)
    
    
    
    ###############################################################################
    #                                                                             #
    #                          Quantification methods                             #
    #                                                                             #
    ###############################################################################
    
    

    ## Volume fraction in the jet (in % of impacting volume)
    
    def compute_JetFrac(self):
        
        meshXci,meshYci,meshVXci,meshVYci = self.velocity_ini('full_div0')

        
        ############### Intersection from trajectories (inside and outside)
        
        trajX,trajY,trajT = self.get_traj()
        
        trajR = np.sqrt(np.square(trajX)+np.square(trajY))
        
        OUT = np.argwhere(np.transpose(trajR > self.Cone.Rcircle))
        
        TimeOut = np.empty(np.shape(meshXci))
        
        
        for i in range(len(meshXci)):
            
            ParticleMask = OUT[:,0] == i
            
            if all(ParticleMask == 0):
            
                TimeOut[i] = len(trajT)-1 # trajectory never got out
                
            else:
                
                TimeOut[i] = np.min(OUT[ParticleMask,1]) # time index of first time out 
            
        
        OutX = np.array([trajX[int(TimeOut[i]),i] for i in range(len(TimeOut))])
        OutY =  np.array([trajY[int(TimeOut[i]),i] for i in range(len(TimeOut))])
        
        OutXco,OutYco = dgf.Circle2Cone(OutX, OutY, self.Cone.Alpha)
        OutRco = np.sqrt(np.square(OutXco)+np.square(OutYco))
        
        InJetMask = (np.abs(OutYco)<0.01) & (OutXco<-self.Cone.Rcone)
        InSecondJetMask = (np.abs(OutYco)<0.01) & (OutXco>-self.Cone.Rcone) & (OutXco<0)
        
        
        OutVelX = np.array([self.trajVXci[int(TimeOut[i]),i] for i in range(len(TimeOut))])
        OutVelY = np.array([self.trajVYci[int(TimeOut[i]),i] for i in range(len(TimeOut))])
        
        
        OutVelXco,OutVelYco = dgf.VelCircle2Cone(OutVelX, OutVelY, OutXco, OutYco, self.Cone.Alpha)

        
        OutJetMask =( np.abs(OutYco)>=0.01) & (OutYco*OutVelYco<0) & (OutRco>self.Cone.Rcone)
        
        
        
        
        ########################### DISPLAY
               
        
        # tx = np.linspace(0,2*np.pi,100)
        
        
        # f = plt.figure(dpi=200)
        
        # plt.plot(self.Cone.Rcone*np.cos(tx),self.Cone.Rcone*np.sin(tx),'g')
        
        # plt.quiver(OutXco[~(InJetMask|OutJetMask|InSecondJetMask)],OutYco[~(InJetMask|OutJetMask|InSecondJetMask)],OutVelXco[~(InJetMask|OutJetMask|InSecondJetMask)],OutVelYco[~(InJetMask|OutJetMask|InSecondJetMask)],color='b')
        
        # plt.quiver(OutXco[InJetMask|OutJetMask],OutYco[InJetMask|OutJetMask],OutVelXco[InJetMask|OutJetMask],OutVelYco[InJetMask|OutJetMask],color='r')
        
        # plt.scatter(OutXco,OutYco,color = 'c',s=4)
        # plt.scatter(self.meshX,self.meshY,color = 'c',label='Non jet trajectories',s=4)
        
        # plt.scatter(OutXco[InJetMask],OutYco[InJetMask],color='crimson',s=4.1,zorder=10)
        # plt.scatter(self.meshX[InJetMask],self.meshY[InJetMask],color = 'crimson',label='Jet formed inside the cone',s=4.1,zorder=10)
        
        
        # plt.scatter(OutXco[OutJetMask],OutYco[OutJetMask],color='blueviolet',s=4.1,zorder=10)
        # plt.scatter(self.meshX[OutJetMask],self.meshY[OutJetMask],color = 'blueviolet',s=4.1,label='Jet formed outside the cone',zorder=10)
        
        # plt.scatter(OutXco[InSecondJetMask],OutYco[InSecondJetMask],color='gold',s=4.1,zorder=10)
        # plt.scatter(self.meshX[InSecondJetMask],self.meshY[InSecondJetMask],color = 'yellow',label='Secondary jet',s=4.1)
        
        
        
        # ax = plt.gca()
        # ax.set_xlim([-1.5*self.Cone.Rcone,1.2*self.Cone.Rcone])
        # ax.set_ylim([-1.2*self.Cone.Rcone,1.2*self.Cone.Rcone])
        # ax.set_aspect('equal')
        
        # plt.legend(fontsize='xx-small')
        
        # f.tight_layout()
        # plt.show()
        
        
        # f = plt.figure(dpi=200)
        
        # plt.plot(self.Cone.Rcircle*np.cos(tx),self.Cone.Rcircle*np.sin(tx),'g')
        
        # # plt.quiver(OutX[~(InJetMask|OutJetMask)],OutY[~(InJetMask|OutJetMask)],OutVelXnorm[~(InJetMask|OutJetMask)],OutVelYnorm[~(InJetMask|OutJetMask)],color='w')
        
        # # plt.quiver(OutX[InJetMask|OutJetMask],OutY[InJetMask|OutJetMask],OutVelXnorm[InJetMask|OutJetMask],OutVelYnorm[InJetMask|OutJetMask],color='c',scale=50)
        
        # plt.scatter(OutX,OutY,color = 'b',label='Non jet trajectories',s=4)
        # plt.scatter(self.meshXci,self.meshYci,color = 'b',s=4)
        
        # plt.scatter(OutX[InJetMask],OutY[InJetMask],color='r',label='Jet formed inside the cone',s=4.1)
        # plt.scatter(self.meshXci[InJetMask],self.meshYci[InJetMask],color = 'r',s=4.1)
        
        
        # plt.scatter(OutX[OutJetMask],OutY[OutJetMask],color='m',label='Jet formed outside the cone',s=4.1)
        # plt.scatter(self.meshXci[OutJetMask],self.meshYci[OutJetMask],color = 'm',s=4.1)
        # # plt.scatter(trajX[:,inter],trajY[:,inter],color = 'w',s=1.1)
        # # plt.scatter(trajX[:,OutJetMask],trajY[:,OutJetMask],color = 'c',s=1.1)
        
        
        
        # ax = plt.gca()
        # ax.set_xlim([-1.5*self.Cone.Rcircle,1.2*self.Cone.Rcircle])
        # ax.set_ylim([-1.2*self.Cone.Rcircle,1.2*self.Cone.Rcircle])
        # ax.set_aspect('equal')
        
        # f.tight_layout()
        # plt.show()
        
    
        ########################### DISPLAY END
        
        
        inter = InJetMask|OutJetMask
        
        
        fracmap = np.zeros(np.shape(self.meshXci))
        
        fracmap[inter] = 1
        fracmap[InSecondJetMask] = 2
        
        self.fracmap = fracmap
        
        self.meshJFxci = self.meshXci[inter]
        self.meshJFyci = self.meshYci[inter] 
        self.meshJFx = self.meshX[inter]
        self.meshJFy = self.meshY[inter]
        
        self.meshJFVx = OutVelXco[inter]
        self.meshJFVy = OutVelYco[inter]
        self.meshJ2FVx = OutVelXco[InSecondJetMask]
        self.meshJ2FVy = OutVelYco[InSecondJetMask]
        self.meshJH = self.meshH[inter] 
        self.meshJ2H = self.meshH[InSecondJetMask] 
        
        
        JetFrac  = np.round(np.sum(inter)/np.size(inter)*1000)/10
        SecondaryJetFrac = np.round(np.sum(InSecondJetMask)/np.size(InSecondJetMask)*1000)/10
        
        
        
        if (self.Drop.Xd>(self.Drop.Rdrop+self.Cone.Rcone)) | ((self.Drop.Rdrop-self.Drop.Xd)>(self.Cone.Rcone)):
            JetFrac = 0
            SecondaryJetFrac = 0
            
        self.JetFrac= JetFrac
        self.SecJetFrac= SecondaryJetFrac
        
    
    def get_JetFrac(self):
        
        if self.MISSED:
            self.JetFrac = 0
            self.SecJetFrac = 0
            
        else:
            
        
            if self.JetFrac == []:
                 self.compute_JetFrac()    
            
            
        
        return(self.JetFrac,self.SecJetFrac)
             
             
    
    # Sheet opening
    
    def compute_SheetOpening(self):
        
        if self.SheetOpen == []:
    
            trajX,trajY,trajT = self.get_traj()
            
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
        
        if self.MISSED:
            self.SheetOpen,self.wiXs,self.wiYs = [0,0,0]
            
        else:
        
            self.compute_SheetOpening()
        
        return(self.SheetOpen,self.wiXs,self.wiYs)
    
    
    def compute_ShapeFactor(self):
        
        ######################## TO CHANGE TO MATCH ANA
        
        if self.MISSED:
            
            self.ShapeFactor = np.nan
        
        else:
            
            trajX_save, trajY_save, trajT_save = self.get_traj()
            
            self.compute_traj(np.linspace(0,15,50))
                
            trajXco, trajYco = self.Cone.Circle2Cone(self.trajX, self.trajY)
            
            time = np.unique(self.trajT)
            
            ShapeFactorTime = np.empty(np.shape(time))
            
            # f,ax = plt.subplots(dpi=200)

            for t,it in zip(time,range(len(time))) :
                
                tmask = self.trajT == t
                
                
                # ax.scatter(trajXco[tmask],trajYco[tmask],color = np.random.rand(1,3))                
                
                Xsize = np.abs(np.nanmin(trajXco[tmask])) - self.Cone.Rcone
                
                Ysize = (np.nanmax(trajYco[tmask]) - np.nanmin(trajYco[tmask]))/2 - self.Cone.Rcone

                Xsize = np.max([Xsize,0])

                Ysize = np.max([Ysize,0])
                
                ShapeFactorTime[it] = Ysize/Xsize
                
            self.ShapeFactor = ShapeFactorTime[-1]

            # f,ax = plt.subplots()
            # ax.set_xlabel('Time [ms]')
            # ax.set_ylabel('ShapeFactor')
            # ax.plot(time,ShapeFactorTime,'-.',label = 'ShapeFactor in time')
            # plt.plot([time[0],time[-1]],[self.ShapeFactor, self.ShapeFactor], label = 'ShapeFactor median')
            # plt.plot([time[0],time[-1]],[ShapeFactorTime[-1], ShapeFactorTime[-1]], label = 'ShapeFactor last value')
            # plt.show()
            
            self.trajX, self.trajY, self.trajT = trajX_save, trajY_save, trajT_save

        
        return(self.ShapeFactor)
        
    def compute_JetVel(self):
        
        # Jet velocity at the edge of the cone at the time it first gets out
        
        if self.MISSED:
            
            self.jetVel = np.nan
        
        else:
            
                
            # Trajectories in cone configuration
            
            trajXco, trajYco = self.Cone.Circle2Cone(self.trajX, self.trajY)
            
            trajVXco,trajVYco =dgf.VelCircle2Cone(self.trajVXci, self.trajVYci, trajXco, trajYco, self.Cone.Alpha)
            
            # Jet trajectories and time abscisse
                        
            mask = (trajXco < -self.Cone.Rcone) & (np.abs(trajYco) < self.Cone.Rcone*0.01)
            
            
            times = self.trajT[mask]
            
            OutTimes = np.sort(np.unique(times))[0:2] # three first time point of a trajectory being out of the cone
            
            timeMask = np.isin(times,OutTimes) 
            
            self.jetVel = np.abs(np.mean(trajVXco[mask][timeMask]))


        
        return(self.jetVel)
        
    
    def compute_JetNRJ(self):
        
        rho = 997e-9 # [kg/mm3]
        
        if (hasattr(self, 'meshJH')) :
            
            if np.sum(self.meshJH)>0:
            
                V2 = np.square(self.meshJFVx)+np.square(self.meshJFVy)

                Veq2 = np.nansum( np.multiply( V2, self.meshJH) ) / np.nansum(self.meshJH)

                
                # Weighted average of trajectories squared velocities by thickness
        
                self.JetNRJ = 4/3*np.pi*self.Drop.Rdrop**3*rho/2*Veq2*self.VolFrac/100*self.get_JetFrac()[0]/100
                self.JetNRJ_Bal = self.JetNRJ*np.sin(2*(np.pi/2-self.Cone.Alpha))
            
            else:
                
                self.JetNRJ = 0
                
                self.JetNRJ_Bal = 0
        
        else:
            
            self.JetNRJ = 0
            
            self.JetNRJ_Bal = 0
        
        return(self.JetNRJ,self.JetNRJ_Bal)
    
    def compute_DispertionDist(self):
        
        g = 9.81*1e-3 # in [mm.ms-2] 
        
        if (hasattr(self, 'meshJH')) :
            
            
            if np.sum(self.meshJH)>0:
                
                V2 = np.square(self.meshJFVx)+np.square(self.meshJFVy)

                Veq2 = np.nansum( np.multiply( V2, self.meshJH) ) / np.nansum(self.meshJH)


                
                self.DispertionDist = Veq2*np.sin(2*(np.pi/2-self.Cone.Alpha))/g
                
                
                self.DispertionHeight = Veq2*np.sin((np.pi/2-self.Cone.Alpha))**2/(2*g)
                
                self.DispertionDistVar = np.sqrt( np.nansum( np.multiply((( np.square(self.meshJFVx)+np.square(self.meshJFVy)) *np.sin(2*(np.pi/2-self.Cone.Alpha)) /g-Veq2)**2,self.meshJH)) / np.nansum(self.meshJH))
                
            else:
                
                self.DispertionDist = 0
                
                self.DispertionHeight = 0
                
                self.DispertionDistVar = 0
        
        else:
            
            self.DispertionDist = 0
            
            self.DispertionHeight = 0
            
            self.DispertionDistVar = 0
        
        
        return(self.DispertionDist,self.DispertionDistVar,self.DispertionHeight)
    
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
            self.compute_JetFrac()
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
           
        self.compute_traj(Time)
        trajX,trajY,trajT = self.get_traj()
            
        trajXco, trajYco = self.Cone.Circle2Cone(trajX, trajY)
    
        order = np.argsort(trajT.flatten())
        
        trajVci = np.sqrt(np.square(self.trajVXci)+ np.square(self.trajVYci))
        
        trajVXco,trajVYco =dgf.VelCircle2Cone(self.trajVXci, self.trajVYci, trajXco, trajYco, self.Cone.Alpha)
        
        
        trajVco = np.sqrt(np.square(trajVXco)+ np.square(trajVYco))
        

        sc0 = ax[0].scatter(trajX.flatten()[order],trajY.flatten()[order],c=trajVci.flatten()[order],cmap = 'Reds',s=0.5,zorder = -1)   
        ax[0].set_box_aspect(1)
        fig.colorbar(sc0, ax = ax[0],orientation='horizontal',label = 'Velocity')
        
        sc1 = ax[1].scatter(trajXco.flatten()[order],trajYco.flatten()[order],c=trajVco.flatten()[order],s=0.5,cmap = 'Reds',zorder = -1)
        
        
        ax[1].scatter(self.wiXs,self.wiYs,s=15,color='r',label='Sheet limits',zorder=4)

        fig.colorbar(sc1, ax = ax[1],orientation='horizontal',label = 'Velocity')
        
        
            
    
        return    
    
    
    # def plot_3Dproj(self,title,resolution):
        
    #     R3D = np.sqrt(self.Drop.meshX3D**2 + self.Drop.meshY3D**2)
    #     inCone3D = R3D<self.Cone.Rcone
        
    #     MeshY3DinCone = self.Drop.meshY3D[inCone3D]
    #     MeshX3DinCone = self.Drop.meshX3D[inCone3D]
    #     MeshZ3DinCone = self.Drop.meshZ3D[inCone3D]+self.Zdrop
        
    #     MeshX2Dcircle,MeshY2Dcircle = self.Cone.Cone2Circle(MeshX3DinCone,MeshY3DinCone)
        
    #     MeshX3Dcircle,MeshY3Dcircle,MeshZ3Dcircle = dgf.Cone2CircleZ(MeshX3DinCone,MeshY3DinCone,MeshZ3DinCone,self.Cone.Alpha)
        
        
    #     ### Circle config : removed sector to form a cone 
                
    #     tx = np.linspace(0,2*np.pi,300)
        
    #     T1 = np.pi-self.Cone.Beta/2
    #     T2 = np.pi+self.Cone.Beta/2
        
    #     if T1>T2:
    #         Ts = np.linspace(np.mod(T1+np.pi,2*np.pi),np.mod(T2+np.pi,2*np.pi),20) - np.pi        
    #     else:                
    #         Ts = np.linspace(T1,T2,20)
        
    #     sectorT = np.append(np.append([T1],Ts),[T2])
    #     sectorR = np.append(np.append([0],self.Cone.Rcircle*np.ones(20)),[0])
        
    #     sectorX,sectorY = vf.ToCart(sectorT,sectorR,angle = 'rad')
        
    #     fig1,[[ax0,ax1],[ax2,ax3]] = plt.subplots(dpi=200,ncols=2,nrows=2)
    #     fig1.suptitle(title)
    #     ax0.set_title('Points')
    #     ax1.set_title('Hist')
    #     ax2.set_title('Interp')
        
    #     ax0.plot(self.Cone.Rcircle*np.cos(tx),self.Cone.Rcircle*np.sin(tx),color = 'g', lw=1,label = 'Circle')
    #     ax0.plot(sectorX,sectorY,'--m',lw=1, label = 'Removed sector')
        
    #     ax0.scatter(MeshX3Dcircle,MeshY3Dcircle,c='b',s=0.5,label='3Dproj',zorder=0)
    #     ax0.scatter(MeshX2Dcircle,MeshY2Dcircle,c='c',s=0.5,label='2Dproj',zorder=1)
        
    #     ax0.set_aspect('equal')
    #     xlim = ax0.get_xlim()
    #     ylim = ax0.get_ylim()
        
        
        
    #     ax1.plot(self.Cone.Rcircle*np.cos(tx),self.Cone.Rcircle*np.sin(tx),color = 'g', lw=1,label = 'Circle')
    #     ax1.plot(sectorX,sectorY,'--m',lw=1, label = 'Removed sector')

    #     ax2.plot(self.Cone.Rcircle*np.cos(tx),self.Cone.Rcircle*np.sin(tx),color = 'g', lw=1,label = 'Circle')
    #     ax2.plot(sectorX,sectorY,'--m',lw=1, label = 'Removed sector')
        
    #     ax3.plot(self.Cone.Rcircle*np.cos(tx),self.Cone.Rcircle*np.sin(tx),color = 'g', lw=1,label = 'Circle')
    #     ax3.plot(sectorX,sectorY,'--m',lw=1, label = 'Removed sector')
        
        
        
    #     I = ax1.hist2d(MeshX3Dcircle,MeshY3Dcircle,zorder=0,label='Drop height (density)',bins=resolution,cmap='plasma')

    #     ax1.set_xlim(xlim)
    #     ax1.set_ylim(ylim)
        
    #     ax1.set_aspect('equal')
        
    #     h,xedges,yedges = I[0],I[1],I[2]
    #     xvalues = (xedges[0:-1] + xedges[1:])/2
    #     yvalues = (yedges[0:-1] + yedges[1:])/2
        
    #     HeightInterp = RegularGridInterpolator((xvalues,yvalues),h,method='cubic')
        
        
    #     xx = np.linspace(np.min(xvalues),np.max(xvalues),50)
    #     yy = np.linspace(np.min(yvalues),np.max(yvalues),50)
        

    #     Xs, Ys = np.meshgrid(xx, yy, indexing='ij')

    #     # interpolator
        
    #     HeightGrid = HeightInterp((Xs, Ys))
        
    #     ax2.scatter(Xs, Ys,c= HeightGrid,cmap='plasma',marker='s')
        
    #     # sns.heatmap(HeightGrid,cmap = 'plasma',ax=ax2)


       
    #     ax2.set_xlim(xlim)
    #     ax2.set_ylim(ylim)
        
    #     ax2.set_aspect('equal')
        
        
    #     Hgrad_x, Hgrad_y = np.gradient(HeightGrid)
        
        

        
    #     xxx = np.linspace(np.min(xvalues),np.max(xvalues),15)
    #     yyy = np.linspace(np.min(yvalues),np.max(yvalues),15)
        
    #     GradXInterp = RegularGridInterpolator((xx,yy), Hgrad_x, method='cubic')
    #     GradYInterp = RegularGridInterpolator((xx,yy), Hgrad_y, method='cubic')
        
        

    #     XXs, YYs = np.meshgrid(xxx, yyy, indexing='ij')

    #     fieldnorm = np.sqrt(GradXInterp((XXs, YYs))**2+ GradYInterp((XXs, YYs))**2)
        
    #     ax3.quiver(XXs, YYs, np.divide(-GradXInterp((XXs, YYs)),fieldnorm), np.divide(-GradYInterp((XXs, YYs)),fieldnorm),fieldnorm,
    #                scale=20,headlength=18,headaxislength=16,headwidth=8,cmap='plasma')

        
    #     ax3.set_xlim(xlim)
    #     ax3.set_ylim(ylim)
        
    #     ax3.set_aspect('equal')
        
    #     fig1.tight_layout()
        
    #     plt.show()
    
     
    def plot_3Dview(self,path,label):
        
        Hmax = 2*self.Drop.Rdrop + self.Cone.Rcone/np.tan(self.Cone.Alpha)
        
        Tmax = Hmax/self.Drop.Vel
        
        Times = np.linspace(0,4*Tmax,np.round(4*Tmax*50).astype(int))
        
        
        R3D = np.sqrt(self.Drop.meshX3D**2 + self.Drop.meshY3D**2)
        inCone3D = R3D<self.Cone.Rcone
        
        MeshY3DinCone = self.Drop.meshY3D[inCone3D]
        MeshX3DinCone = self.Drop.meshX3D[inCone3D]
         
        
        self.compute_traj(Times)
        trajX,trajY,trajT = self.get_traj()
            
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
            fig.suptitle(label)
            ax0 = fig.add_subplot(121, projection='3d')
            ax0.view_init(elev=25, azim=-100)
            ax0.set_title('Circle config')
            # ax0.set_xlim(-4,4)
            # ax0.set_ylim(-4,4)
            # ax0.set_zlim(0,3)
            ax1 = fig.add_subplot(122, projection='3d')
            # ax1.set_xlim(-2.5,2.5)
            # ax1.set_ylim(-2.5,2.5)
            # ax1.set_zlim(0,5)
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
            
    
            
            
            savepath = path + label
            
            os.makedirs(savepath,exist_ok = True)
            
            figname = savepath + '\Time_{:.2f}.png'
            
            
            fig.savefig(figname.format(tt))
            
            plt.close(fig)
            
            
            
    def rotating_3D(self,path,label):
        
        Hmax = 2*self.Drop.Rdrop + self.Cone.Rcone/np.tan(self.Cone.Alpha)
        
        Tmax = Hmax/self.Drop.Vel
        
        Times = np.linspace(0,4*Tmax,np.round(4*Tmax*50).astype(int))
        
        
        R3D = np.sqrt(self.Drop.meshX3D**2 + self.Drop.meshY3D**2)
        inCone3D = R3D<self.Cone.Rcone
        
        MeshY3DinCone = self.Drop.meshY3D[inCone3D]
        MeshX3DinCone = self.Drop.meshX3D[inCone3D]
        
        
        MeshY3DOutCone = self.Drop.meshY3D[~inCone3D]
        MeshX3DOutCone = self.Drop.meshX3D[~inCone3D]
        
        
        self.compute_traj(Times)
        trajX,trajY,trajT = self.get_traj()
            
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
            
            savepath = path + label
            
            os.makedirs(savepath,exist_ok = True)
            
            figname = savepath + '\Right_Time_{:.2f}.png'
            
            
            fig.savefig(figname.format(tt))
            
            plt.close(fig)
        
        
    
    def plot_splash_movie(self,Time,path,label,xlims,ylims,**kwargs):
        
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
                
        savepath = path + label 
         
        os.makedirs(savepath,exist_ok=True)
        
        self.compute_traj(Time)
        trajX,trajY,trajT = self.get_traj()
        
        self.compute_JetFrac()
            
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
            
            ax[0].scatter(trajX[points_mask][heightmask],trajY[points_mask][heightmask],c=self.fracmap[heightmask],cmap='rainbow',s = 0.5,vmin=0,vmax=self.Drop.Rdrop)
            
            ax[1].scatter(trajXco[points_mask][heightmask],trajYco[points_mask][heightmask],c=self.fracmap[heightmask],cmap='rainbow',s = 0.5,vmin=0,vmax=self.Drop.Rdrop)
            
                
            for i in range(iT):
                
                i= i+1
                
                Ti = timevect[i]
                

                points_mask = (trajT == timevect[iT-i])
                
                
                
                Hi = Ti*self.Drop.Vel 
                
                
                
                heightmask = np.zeros(np.shape(self.meshH), dtype=bool)
                heightmask = self.meshH/2 > np.abs(Hi - contactDist)
                

                    
                
                
                ax[0].scatter(trajX[points_mask][heightmask],trajY[points_mask][heightmask],c=self.fracmap[heightmask],cmap='rainbow',vmin=0,vmax=self.Drop.Rdrop,s = 1,zorder=10)
                ax[1].scatter(trajXco[points_mask][heightmask],trajYco[points_mask][heightmask],c=self.fracmap[heightmask],cmap='rainbow',vmin=0,vmax=self.Drop.Rdrop,s = 1,zorder=10)
                
            ax[1].set_xlim(xlims)
            ax[1].set_ylim(ylims)
                
            figname = savepath + '\Time_' + str(np.round(T*1000)) + '.png'
            
            plt.savefig(figname)
            
            plt.close(fig)
           
    
        return    
    
    
    
    
    