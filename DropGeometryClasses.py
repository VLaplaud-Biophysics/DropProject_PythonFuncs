# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 09:11:18 2023

@author: laplaud
"""
import numpy as np
import numpy.matlib as ml
import matplotlib.pyplot as plt

from scipy.interpolate import LinearNDInterpolator

import DropGeometryFuncs as dgf

import VallapFunc_DP as vf

import time as time


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
        
        isInCone = R<self.Rcone
        
        return(isInCone)
    
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
    
    def impact(self,drop,oriType,velIni):
        
        I = Impact(drop,self,oriType,velIni)
        
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
        
        ax[0].set_aspect('equal')
        ax[1].set_aspect('equal')
        
        # (0,0) point
        ax[0].plot(0,0,'g*',ms = 3)
        ax[1].plot(0,0,'g*',ms = 3)
        

        ax[0].plot(self.Rcircle*np.cos(tx),self.Rcircle*np.sin(tx),color = ConeColor, lw=ConeLW,label = 'Cone',zorder=1);
        ax[0].plot(sectorX,sectorY,'--r',lw=1, label = 'Removed sector',zorder=2)

        ax[1].plot(self.Rcone*np.cos(tx),self.Rcone*np.sin(tx),color = ConeColor, lw=ConeLW, label = 'Cone',zorder=1);
        
        
        
        if not Drop == None:
            
            
            ax[0].plot(Drop.Xd,Drop.Yd,'b*',ms = 3)
            ax[1].plot(Drop.Xd,Drop.Yd,'b*',ms = 3)
            
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
        
        ### Mesh of points in the drop
        
        
        ## Carthesian fixed mesh + borderpoints
        
        Xmin, Xmax, Ymin, Ymax = self.Xd-self.Rdrop,self.Xd+self.Rdrop,self.Yd-self.Rdrop,self.Yd+self.Rdrop
        
        Xs = np.linspace(Xmin, Xmax, npts)
        Ys = np.linspace(Ymin, Ymax, npts)
        
        Xgrid,Ygrid = np.meshgrid(Xs,Ys)
        

        XgridF = Xgrid.flatten()
        YgridF = Ygrid.flatten()

        Agrid,Rgrid = vf.ToCirc(XgridF,YgridF)

        isDrop = self.IsIn(Rgrid,Agrid)
        
        Rmax = self.Rdrop
        
        RgridBorder = np.ones((1,2*npts))*Rmax
        AgridBorder = np.linspace(0,2*np.pi,2*npts)
        
        XgridBorder,YgridBorder = vf.ToCart(AgridBorder,RgridBorder,angle='rad')
        
        XgridF = np.append(XgridF[isDrop],XgridBorder+self.Xd)
        YgridF = np.append(YgridF[isDrop],YgridBorder+self.Yd)
        

        self.meshXFull = Xgrid # in cone config
        self.meshYFull = Ygrid # in cone config
        
        self.meshX = XgridF # in cone config
        self.meshY = YgridF # in cone config
        self.meshH = dgf.SphereH(self.Rdrop,self.meshX,self.meshY,self.Xd)
        
      
        
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
    
    def impact(self,cone,oriType,velIni):
        
        I = Impact(self,cone,oriType,velIni)
        
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
            if velType == 'Radial':
                raise ValueError('Invalid case ! A single point mesh is incompatible with velocity type "Radial"')
                    

        
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
        
        self.meshVXci = []
        self.meshVYci = []
        
        
        
        self.SheetOpen = []
        self.wiXs = []
        self.wiYs = []
        
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
            
            meshX,meshY,meshH = self.Drop.mesh() # Cone config
            
            InCone = np.sqrt(np.square(meshX) + np.square(meshY))<self.Cone.Rcone
            
            self.meshX = meshX[InCone]
            self.meshY = meshY[InCone]
            self.meshH = meshH[InCone]
            
            meshXci,meshYci = self.Cone.Cone2Circle(self.meshX, self.meshY) 
            
            if self.meshType == 'point':
                nmesh = np.numel(meshX)
                angles = np.linspace(0,np.pi*2,nmesh)
                R = self.Drop.Rdrop/100
                
                pointmeshX = R*np.cos(angles)
                pointmeshY = R*np.sin(angles)
            
            
            if self.oriType == 'Hgrad':
            
                HgradInterp_x,HgradInterp_y = self.Drop.Hgradient()       
                
                meshOX = -HgradInterp_x(self.meshX,self.meshY) 
                meshOY = -HgradInterp_y(self.meshX,self.meshY) 
                
                
                mesh1 = meshXci-meshXci[np.argmax(self.meshH)]
                mesh2 = meshYci-meshYci[np.argmax(self.meshH)]
                self.meshDist = np.sqrt(np.square(mesh1)+np.square(mesh2))
                
                self.impactR = np.max(self.meshDist)
                
                # to circle config
                meshOXci,meshOYci = dgf.VelCone2Circle(meshOX, meshOY,self.meshX,self.meshY, self.Cone.Alpha) # Circle config, impacting fraction
    
            elif self.oriType == 'Drop':
                
                Xdci,Ydci = self.Cone.Cone2Circle(self.Drop.Xd, self.Drop.Yd)
                
                meshOXci = meshXci-Xdci
                meshOYci = meshYci-Ydci
                
                self.meshDist = np.sqrt(np.square(meshOXci)+np.square(meshOYci))
                self.impactR = np.max(self.meshDist)
                
                
                
                
            elif self.oriType == 'Hmax':
                
                oriX = meshXci[np.argmax(self.meshH)]
                oriY = meshYci[np.argmax(self.meshH)]
                
                if self.meshType == 'zone':
                    meshOXci = meshXci-oriX
                    meshOYci = meshYci-oriY
                
                    self.meshDist = np.sqrt(np.square(meshOXci)+np.square(meshOYci))
                    self.impactR = np.max(self.meshDist)
                
                elif self.meshType == 'point':
                    meshOXci = pointmeshX
                    meshOYci = pointmeshY
                    
                    meshXci = meshOXci + oriX
                    meshYci = meshOYci + oriY
            
                
            elif self.oriType == 'Central':
                
                oriX = meshXci.mean()
                oriY = meshYci.mean()
                
                if self.meshType == 'zone':
                    meshOXci = meshXci-oriX
                    meshOYci = meshYci-oriY
                
                    self.meshDist = np.sqrt(np.square(meshOXci)+np.square(meshOYci))
                    self.impactR = np.max(self.meshDist)
                
                elif self.meshType == 'point':
                    meshOXci = pointmeshX
                    meshOYci = pointmeshY
                    
                    meshXci = meshOXci + oriX
                    meshYci = meshOYci + oriY
                
            elif self.oriType == 'CentralHeight':
                
                oriX =  np.sum(np.multiply(meshXci,self.meshH))/np.sum(self.meshH)
                oriY = np.sum(np.multiply(meshYci,self.meshH))/np.sum(self.meshH)
                
                if self.meshType == 'zone':
                    meshOXci = meshXci-oriX
                    meshOYci = meshYci-oriY
                
                    self.meshDist = np.sqrt(np.square(meshOXci)+np.square(meshOYci))
                    self.impactR = np.max(self.meshDist)
                
                elif self.meshType == 'point':
                    meshOXci = pointmeshX
                    meshOYci = pointmeshY
                    
                    meshXci = meshOXci + oriX
                    meshYci = meshOYci + oriY
            
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
        
        if self.meshVXci == []:
    
            meshXci,meshYci,meshOXci,meshOYci = self.orientation() # Cone config, full drop       

            
            normCi = np.sqrt(np.square(self.meshOXci)+np.square(self.meshOYci))/np.sqrt(2)
                  
            if self.velIni == 'VelNorm':
                self.meshVXci_norm = np.divide(self.meshOXci,normCi)*self.velnorm/np.sqrt(2)
                self.meshVYci_norm = np.divide(self.meshOYci,normCi)*self.velnorm/np.sqrt(2)
                
            elif self.velIni == 'Radial':
                self.meshVXci_norm = np.multiply(np.divide(self.meshOXci,normCi),self.meshDist)/self.impactR*self.velnorm/np.sqrt(2) 
                self.meshVYci_norm = np.multiply(np.divide(self.meshOYci,normCi),self.meshDist)/self.impactR*self.velnorm/np.sqrt(2) 
######################################################################################################################################
#######################################################################################################################################
            
            meshX,meshY = self.Cone.Circle2Cone(meshXci, meshYci)
            T,R = vf.ToCirc(meshX,meshY, angle = 'rad') 
            
            
            meshVX_tan = -np.cos(T)*self.veltan*np.sin(self.Cone.Alpha)
            meshVY_tan = -np.sin(T)*self.veltan*np.sin(self.Cone.Alpha)
            
            Tci,Rci = vf.ToCirc(meshXci,meshYci, angle = 'rad') 
            self.meshVXci_tan = -np.cos(Tci)*self.veltan
            self.meshVYci_tan = -np.sin(Tci)*self.veltan
            
            meshVX_tan_div0 = np.ones(np.shape(meshVX_tan))*np.mean(meshVX_tan)
            meshVY_tan_div0 = np.ones(np.shape(meshVY_tan))*np.mean(meshVY_tan)
            
            
            self.meshVXci_tan_div0,self.meshVYci_tan_div0 = dgf.VelCone2Circle(meshVX_tan_div0,meshVY_tan_div0, meshX, meshY, self.Cone.Alpha)

            

            self.meshVXci = self.meshVXci_norm + self.meshVXci_tan
            self.meshVYci = self.meshVYci_norm + self.meshVYci_tan
            
            self.meshVXci_div0 = self.meshVXci_norm + self.meshVXci_tan_div0
            self.meshVYci_div0 = self.meshVYci_norm + self.meshVYci_tan_div0
        
        
    
    def velocity_ini(self,veltype):
        
        self.compute_velocity_ini()
        
        if veltype == 'norm':
            
            self.meshVX,self.meshVY = dgf.VelCone2Circle(self.meshVXci_norm,self.meshVYci_norm, self.meshX, self.meshY, self.Cone.Alpha)
            
            return(self.meshXci,self.meshYci,self.meshVXci_norm,self.meshVYci_norm)
        
        elif veltype == 'tan':
            
            self.meshVX,self.meshVY = dgf.VelCone2Circle(self.meshVXci_tan,self.meshVYci_tan, self.meshX, self.meshY, self.Cone.Alpha)
            
            return(self.meshXci,self.meshYci,self.meshVXci_tan,self.meshVYci_tan)
        
        elif veltype == 'tan_div0':
            
            self.meshVX,self.meshVY = dgf.VelCone2Circle(self.meshVXci_tan_div0,self.meshVYci_tan_div0, self.meshX, self.meshY, self.Cone.Alpha)
            
            return(self.meshXci,self.meshYci,self.meshVXci_tan_div0,self.meshVYci_tan_div0)
        
        elif veltype == 'full':   
            
            self.meshVX,self.meshVY = dgf.VelCone2Circle(self.meshVXci,self.meshVYci, self.meshX, self.meshY, self.Cone.Alpha)         
            
            return(self.meshXci,self.meshYci,self.meshVXci,self.meshVYci)
        
        elif veltype == 'full_div0': 
            
            self.meshVX,self.meshVY = dgf.VelCone2Circle(self.meshVXci_div0,self.meshVYci_div0, self.meshX, self.meshY, self.Cone.Alpha)           
            
            return(self.meshXci,self.meshYci,self.meshVXci_div0,self.meshVYci_div0)
        
        

    
    
    def compute_traj(self,Time,velType):
       
        meshXci,meshYci,meshVXci,meshVYci = self.velocity_ini(velType)
               
        trajT = ml.repmat(Time,len(meshXci),1).T
        
        meshPts_X = ml.repmat(meshXci,len(Time),1)
        meshPts_Y = ml.repmat(meshYci,len(Time),1)
        
        meshPts_VX = ml.repmat(meshVXci,len(Time),1)
        meshPts_VY = ml.repmat(meshVYci,len(Time),1)     
        
        trajX = meshPts_X + np.multiply(meshPts_VX,trajT)
        trajY = meshPts_Y + np.multiply(meshPts_VY,trajT)
        
       
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
        
        lx = -(meshPts_X[badPts]*np.tan(Beta/2)+np.abs(meshPts_Y[badPts]))/(np.abs(np.sin(badT))-np.abs(np.cos(badT))*np.tan(Beta/2))
        
        lz = (lx*np.abs(np.sin(badT))+np.abs(meshPts_Y[badPts]))/np.sin(Beta/2)

        # Sorting the bad point on both side of the removed sector
        closeT1 = np.abs((np.mod(badT,2*np.pi)-T1))<np.abs((np.mod(badT,2*np.pi)-T2))
        badT[closeT1] = T1
        badT[~closeT1] = T2

        trajX[badPts],trajY[badPts] = vf.ToCart(badT,badR-lx+lz)
        
        trajX[badPts],trajY[badPts] = trajX[badPts]+self.Cone.Xc,trajY[badPts]+self.Cone.Yc    
        
        return(trajX,trajY,trajT)
    
    
    
    
    ## Volume fraction in the jet (in % of impacting volume)
    
    def compute_JetFrac(self,veltype):
        
        meshXci,meshYci,meshVXci,meshVYci = self.velocity_ini(veltype)
            
        # equations for the lines of the sector borders (y = c12*x)
        c1 = np.tan(np.pi-self.Cone.Beta/2)
        c2 = np.tan(-np.pi+self.Cone.Beta/2)
        
        # equation for the line along the trajectory (y = a*x + b)
        a = np.divide(meshVYci,meshVXci)
        b = np.divide((meshVXci*meshYci - meshVYci*meshXci),meshVXci)
        
        # intersection points x coord
        xi1 = b/(c1-a)
        xi2 = b/(c2-a)
        
        yi1 = c1*xi1
        yi2 = c2*xi2
        
        
        # intersection 
        
        if self.Cone.Alpha>np.pi/6:
            inter = (meshVXci <=0)&(((xi1<=0)&(xi1>=-1*self.Cone.Rcircle*np.cos(self.Cone.Beta/2)))|((xi2<=0)&(xi2>=-1*self.Cone.Rcircle*np.cos(self.Cone.Beta/2))))
        elif  self.Cone.Alpha<np.pi/6: 
            inter =  ((np.divide(meshVXci,xi1-meshXci)>=0)|(np.divide(meshVXci,xi2-meshXci)>=0))&(((xi1>=0)&(xi1<=-1*self.Cone.Rcircle*np.cos(self.Cone.Beta/2)))|((xi2>=0)&(xi2<=-1*self.Cone.Rcircle*np.cos(self.Cone.Beta/2))))
        else:
            inter = (meshVXci <=0)&((np.abs(yi1)<self.Cone.Rcircle)|(np.abs(yi2)<self.Cone.Rcircle))
        
        
        
        self.meshJFxci = self.meshXci[inter]
        self.meshJFyci = self.meshYci[inter] 
        self.meshJFx = self.meshX[inter]
        self.meshJFy = self.meshY[inter]
        self.meshJFVx = self.meshVX[inter]
        self.meshJFVy = self.meshVY[inter]
        self.meshJH = self.meshH[inter] 
        JetFrac  = np.round(np.sum(inter)/np.size(inter)*1000)/10
        
        if (np.sqrt(np.square(self.Drop.Xd)+np.square(self.Drop.Yd))>(self.Drop.Rdrop+self.Cone.Rcone)) | ((self.Drop.Rdrop-np.sqrt(np.square(self.Drop.Xd)+np.square(self.Drop.Yd)))>(self.Cone.Rcone)):
            JetFrac = 0
        
        return(JetFrac)
    
            
            
    
    # Sheet opening
    
    def compute_SheetOpening(self,velType):
        
        if self.SheetOpen == []:
    
            trajX,trajY,trajT = self.compute_traj(np.linspace(0,50,500),velType)
            
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
                
    def SheetOpening(self,velType):
        
        self.compute_SheetOpening(velType)
        
        return(self.SheetOpen,self.wiXs,self.wiYs)
    
    def compute_JetNRJ(self,velType):
        
        rho = 997e-9 # [kg/mm3]
        
        Veq2 = np.sum(np.multiply(self.meshJH,np.sqrt(np.square(self.meshJFVx)+np.square(self.meshJFVy))))/np.sum(self.meshJH)
        
        # print(np.sqrt(Veq2))
        self.JetNRJ = 4/3*np.pi*self.Drop.Rdrop**3*rho/2*self.VolFrac/100*self.compute_JetFrac(velType)/100*Veq2
        
        return(self.JetNRJ)
    
    ###############################################################################
    #                                                                             #
    #                          Display methods                                    #
    #                                                                             #
    ###############################################################################
    
    def plot_splash_init(self,**kwargs):
        
        Title     = 'kw: title= ''My title for the figure'''
        Xlabel_Ci = 'kw: xlabelCi= ''My xlabel for the circle config'''
        Ylabel_Ci = 'kw: ylabelCi= ''My ylabel for the circle config'''
        Xlabel_Co = 'kw: xlabelCo= ''My xlabel for the cone config'''
        Ylabel_Co = 'kw: ylabelCo= ''My ylabel for the cone config'''
        NoLabels  = False
        VelType = 'full'
        
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
            elif key == 'veltype':
                VelType = value

                    
            else:
                print('Unknown key : ' + key + '. Kwarg ignored.')
        
       
        meshXci,meshYci,meshVXci,meshVYci = self.velocity_ini(VelType)

        fieldnormCi = np.sqrt(np.square(meshVXci)+np.square(meshVYci))  
        # fieldnormCi = 4
        
        meshX, meshY = self.Cone.Circle2Cone(meshXci, meshYci)
        
        # speedres = 10000

        # meshVX, meshVY = self.Cone.Circle2Cone(meshVXci/speedres+meshXci, meshVYci/speedres+meshYci)
        
        # meshVX = (meshVX - meshX)*speedres
        # meshVY = (meshVY - meshY)*speedres
        
        
        meshVX,meshVY = dgf.VelCircle2Cone(meshVXci,meshVYci, meshX, meshY, self.Cone.Alpha)
        
        fieldnorm = np.sqrt(np.square(meshVX)+np.square(meshVY))     
        # fieldnorm = 4

    
        
        fig, ax = self.Cone.draw(drop=self.Drop,nolabels=NoLabels,dropview = 'impact',conelinewidth=ConeLW,
                                 conecolor=ConeColor,title= Title + '_' + VelType,xlabelCi=Xlabel_Ci,ylabelCi=Ylabel_Ci,xlabelCo=Xlabel_Co,ylabelCo=Ylabel_Co)
        
        self.compute_JetFrac(VelType)
        ax[0].scatter(self.meshJFxci,self.meshJFyci,c='r',zorder=5,s=7)
        ax[1].scatter(self.meshJFx,self.meshJFy,c='r',zorder=5,s=7)
        
        q0 = ax[0].quiver(meshXci, meshYci, np.divide(meshVXci,fieldnormCi), np.divide(meshVYci,fieldnormCi),fieldnormCi,scale = 20,zorder=5,headlength=18,headaxislength=16)
        q1 = ax[1].quiver(meshX,   meshY,   np.divide(meshVX,fieldnorm),     np.divide(meshVY,fieldnorm),    fieldnorm,  scale = 15,zorder=20,headlength=18,headaxislength=16)

        fig.colorbar(q0, ax = ax[0],orientation='horizontal',label = 'Velocity',ticks=np.linspace(0,self.Drop.Vel,6))
        fig.colorbar(q1, ax = ax[1],orientation='horizontal',label = 'Velocity',ticks=np.linspace(0,self.Drop.Vel,6))
    
    def plot_div(self,**kwargs):
        
        Title     = 'kw: title= ''My title for the figure'''
        VelType = 'full'
        
        ConeColor = 'g'
        ConeLW = 1
        
        
        
        for key, value in kwargs.items(): 
            if key == 'title':
                Title = value
            elif key == 'conecolor':
                ConeColor = value
            elif key == 'conelinewidth':
                ConeLW = value
            elif key == 'veltype':
                VelType = value

                    
            else:
                print('Unknown key : ' + key + '. Kwarg ignored.')
        
       
        meshXci,meshYci,meshVXci,meshVYci = self.velocity_ini(VelType)
        
        meshX,meshY = self.Drop.meshX,self.Drop.meshY
        Xgrid,Ygrid = self.Drop.meshXFull,self.Drop.meshYFull
        
        Xgridci,Ygridci = self.Cone.Cone2Circle(Xgrid,Ygrid)
        
        VXgridci = np.zeros(np.shape(Xgrid))
        # VXgridci[:] = np.nan
        VYgridci = np.zeros(np.shape(Xgrid))
        # VYgridci[:] = np.nan
        
        
        Maskgridci = np.zeros(np.shape(Xgrid))
        
        
        for pos,vx,vy in zip([np.argwhere((Xgrid == x)&(Ygrid == y)) for x,y in zip(meshX,meshY)],meshVXci,meshVYci):
    
            VXgridci[pos[0][0],pos[0][1]] = vx
            VYgridci[pos[0][0],pos[0][1]] = vy
            Maskgridci[pos[0][0],pos[0][1]] = 1
         
            
         # HgradInterp_x = LinearNDInterpolator(list(zip(Xgrid.flatten(),Ygrid.flatten())),Hgrad_x.flatten())
        
        dVXci = np.gradient(VXgridci, axis=1)[Maskgridci.astype(bool)]
        dVYci = np.gradient(VXgridci, axis=0)[Maskgridci.astype(bool)]
        divci = dVXci+dVYci
        
        fieldnormCi = np.sqrt(np.square(meshVXci)+np.square(meshVYci))        
        
        meshX, meshY = self.Cone.Circle2Cone(meshXci, meshYci)
        
        meshVX, meshVY = self.Cone.Circle2Cone(meshVXci+meshXci, meshVYci+meshYci)
        
        meshVX = meshVX - meshX
        meshVY = meshVY - meshY

        
        ### An abscisse vector
                
        tx = np.linspace(0,2*np.pi,100)
        
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
        
        ### Adding the drop 

        Xdrop,Ydrop = self.Drop.Coordinates() # in cone config
            
        Rd = self.Drop.Parameters()[0]               
    
        ### Anaytical drop equation in circle config
        
        if Xdrop > Rd:
            ThetaBorder = np.sin(self.Cone.Alpha)*np.arcsin(Rd/Xdrop)
            Theta = np.linspace(-ThetaBorder,ThetaBorder,200)
            R1 = (Xdrop*np.cos(np.divide(Theta,np.sin(self.Cone.Alpha))) - np.sqrt(Rd**2-Xdrop**2*np.sin(np.divide(Theta,np.sin(self.Cone.Alpha)))**2))/(np.sin(self.Cone.Alpha))
            R2 = (Xdrop*np.cos(np.divide(Theta,np.sin(self.Cone.Alpha))) + np.sqrt(Rd**2-Xdrop**2*np.sin(np.divide(Theta,np.sin(self.Cone.Alpha)))**2))/(np.sin(self.Cone.Alpha))
            Theta = np.concatenate((Theta, np.flip(Theta)))
            R = np.concatenate((R1, R2))
        else:
            ThetaBorder = np.pi - self.Cone.Beta/2
            Theta = np.linspace(-ThetaBorder,ThetaBorder,200)
            R = (Xdrop*np.cos(np.divide(Theta,np.sin(self.Cone.Alpha))) + np.sqrt(Rd**2-Xdrop**2*np.sin(np.divide(Theta,np.sin(self.Cone.Alpha)))**2))/(np.sin(self.Cone.Alpha))
        
        
        X,Y = vf.ToCart(Theta,R, angle = 'rad')

        
        ### Plotting
        
        fig,ax = plt.subplots(dpi = 250, ncols = 3, figsize = (9,5))
        
        fig.suptitle(Title)
        ax[0].set_title('dVx/dx')
        ax[1].set_title('dVy/dy')
        ax[2].set_title('div(V)')

        
        for axx in ax:
            
            axx.set_aspect('equal')
            
            axx.set_xticks([])
            axx.set_yticks([])
            
            axx.plot(0,0,'g*',ms = 3)
            axx.plot(self.Cone.Rcircle*np.cos(tx),self.Cone.Rcircle*np.sin(tx),color = ConeColor, lw=ConeLW,label = 'Cone',zorder=1);
            axx.plot(sectorX,sectorY,'--r',lw=1, label = 'Removed sector',zorder=2)

        
        
        q0 = ax[0].quiver(meshXci, meshYci, np.divide(meshVXci,fieldnormCi), np.divide(meshVYci,fieldnormCi),dVXci,scale = 20,zorder=5,headlength=18,headaxislength=16)

        q1 = ax[1].quiver(meshXci, meshYci, np.divide(meshVXci,fieldnormCi), np.divide(meshVYci,fieldnormCi),dVYci,scale = 20,zorder=5,headlength=18,headaxislength=16)

        q2 = ax[2].quiver(meshXci, meshYci, np.divide(meshVXci,fieldnormCi), np.divide(meshVYci,fieldnormCi),divci,scale = 20,zorder=5,headlength=18,headaxislength=16)

        fig.colorbar(q0, ax = ax[0],orientation='horizontal',label = 'dVx/dx')
        fig.colorbar(q1, ax = ax[1],orientation='horizontal',label = 'dVy/dy')
        fig.colorbar(q2, ax = ax[2],orientation='horizontal',label = 'Div(V)')
        
        fig.tight_layout()
        
    
    def plot_splash_traj(self,Time,velType,**kwargs):
        
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
                                 conecolor=ConeColor,title='JetFrac = ' + str(self.compute_JetFrac(velType)) + '. '+ Title,xlabelCi=Xlabel_Ci,ylabelCi=Ylabel_Ci,xlabelCo=Xlabel_Co,ylabelCo=Ylabel_Co)
           
                    
        trajX,trajY,trajT = self.compute_traj(Time,velType)
            
        trajXco, trajYco = self.Cone.Circle2Cone(trajX, trajY)
    
        order = np.argsort(trajT.flatten())

        sc0 = ax[0].scatter(trajX.flatten()[order],trajY.flatten()[order],c=trajT.flatten()[order],s=2,zorder = -1)   
        ax[0].set_box_aspect(1)
        fig.colorbar(sc0, ax = ax[0],orientation='horizontal',label = 'Time')
        
        sc1 = ax[1].scatter(trajXco.flatten()[order],trajYco.flatten()[order],c=trajT.flatten()[order],s=2,zorder = -1)
        
        
        ax[1].scatter(self.wiXs,self.wiYs,s=15,color='r',label='Sheet limits')
        
        # ax[1].set_box_aspect(1)
        
        
        ax[1].set_xlim(np.array([-10,7]))
        ax[1].set_ylim(np.array([-7,7]))
        
        fig.colorbar(sc1, ax = ax[1],orientation='horizontal',label = 'Time')
    
        return(fig)
    
    
    
    
    
    