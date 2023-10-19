# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 09:11:18 2023

@author: laplaud
"""
import numpy as np
import matplotlib.pyplot as plt

import DropGeometryFuncs as dgf

import VallapFunc_DP as vf


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
        
        return(self.Rcircle,self.Beta,self.Frac)
    
    def IsIn(self,R,Theta): 
        
        isInCone = R<self.Rcone
        
        return(isInCone)
    
    
    def Cone2Circle(self,X,Y):
        # (X,Y) points to transform, Alpha angle of the cone, Ad = 0 angle of removed sector bissecant (cone config)
    
        return(dgf.Cone2Circle(X,Y,self.Alpha,0))
    
    def Circle2Cone(self,X,Y): 
        # (X,Y) points to transform, Alpha angle of the cone, Ad = 0 angle of removed sector bissecant (circle config)
        
        return(dgf.Circle2Cone(X,Y,self.Alpha,0))
    
    
    
    ######  Impact with a drop
    
    def impact(self,drop):
        
        I = Impact(drop,self)
        
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
        
        fig,ax = plt.subplots(dpi = 250, ncols = 2, figsize = (9,4))
        
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
            
            meshX,meshY,meshH = Drop.mesh()
            droplabel = 'Drop height'
            
            if DropView == 'impact':
                meshA,meshR = vf.ToCirc(meshX,meshY, angle='deg')
                inImpact = meshR<self.Rcone
                meshX,meshY,meshH = meshX[inImpact],meshY[inImpact],meshH[inImpact]
                droplabel = 'Drop height (impacting fraction)'
            
            cmap = plt.get_cmap('Blues')
            bluemap = vf.truncate_colormap(cmap,0.5,1,100)
            
            ax[1].scatter(meshX,meshY,c=meshH,cmap=bluemap, s=2, zorder=3,label=droplabel)
            
            meshXci,meshYci = self.Cone2Circle(meshX, meshY)
            
            ax[0].scatter(meshXci,meshYci,c=meshH,cmap=bluemap, s=2, zorder=3,label=droplabel)
            
            ax[0].plot(X,Y, 'b-', lw = 1,label='Deformed drop',zorder=4)      
            ax[1].plot(Xdrop + Rd*np.cos(tx), Ydrop + Rd*np.sin(tx), 'b-', lw = 1, label='Drop',zorder=4)
        
        ax[0].legend(fontsize='xx-small')
        ax[1].legend(fontsize='xx-small')
        fig.tight_layout()
        
        return(fig,ax)
    
    
############################  
# 2. A class for the drop. #
############################


class Drop:
    
    
    def __init__(self, radius, offcent, npts, speed):

        self.Rdrop = radius
        
        self.Xd = offcent
        
        self.Yd = 0
        
        self.Speed = speed
        
        self.Npts = npts
        
        ### Mesh of points in the drop
                
        Xmin, Xmax, Ymin, Ymax = self.Xd-self.Rdrop,self.Xd+self.Rdrop,self.Yd-self.Rdrop,self.Yd+self.Rdrop
        
        Xs = np.linspace(Xmin, Xmax, npts)
        Ys = np.linspace(Ymin, Ymax, npts)
        
        Xgrid,Ygrid = np.meshgrid(Xs,Ys)

        Xgrid = Xgrid.flatten()
        Ygrid = Ygrid.flatten()

        Agrid,Rgrid = vf.ToCirc(Xgrid,Ygrid)

        isDrop = self.IsIn(Rgrid,Agrid)

        self.meshX = Xgrid[isDrop] # in cone config
        self.meshY = Ygrid[isDrop] # in cone config
        self.meshH = dgf.SphereH(self.Rdrop,self.meshX,self.meshY,self.Xd)
        
        
        
    ###### Geometry related methods
    
    def Coordinates(self):
        
        return(self.Xd,self.Yd)
        
    def Parameters(self):
        
        return(self.Rdrop,self.Npts)
    
    def IsIn(self,R,Theta): 
        
        isInDrop = R**2 - 2*R*(self.Xd*np.cos(Theta)+self.Yd*np.sin(Theta)) + self.Xd**2 + self.Yd**2 < self.Rdrop**2
        
        return(isInDrop)
        
    def mesh(self):
        
        return(self.meshX,self.meshY,self.meshH)
      

    ######  Impact with a cone
    
    def impact(self,cone):
        
        I = Impact(self,cone)
        
        return(I)

        
#################################################    
# 3. A class for an impact of a drop on a cone. #
################################################# 

class Impact:
    
    
    def __init__(self,drop,cone):
        
        self.Drop = drop
        
        self.Cone = cone
        
        self.speednorm = drop.Speed*np.cos(self.Cone.Alpha)
        
        self.speedtan = drop.Speed*np.sin(self.Cone.Alpha)
    
    
    def orientation(self):
        
        meshX,meshY,meshH = self.Drop.mesh() # Cone config
        
        Xmin, Xmax, Ymin, Ymax = self.Drop.Xd-self.Drop.Rdrop,self.Drop.Xd+self.Drop.Rdrop,self.Drop.Yd-self.Drop.Rdrop,self.Drop.Yd+self.Drop.Rdrop
        
        Xs = np.linspace(Xmin, Xmax, self.Drop.Npts)
        Ys = np.linspace(Ymin, Ymax, self.Drop.Npts)
        
        Xs = np.insert(Xs,0,Xmin*0.9)
        Xs = np.append(Xs,Xmin*1.1)
        Ys = np.insert(Ys,0,Ymin*0.9)
        Ys = np.append(Ys,Ymin*1.1)
        
        Xgrid,Ygrid = np.meshgrid(Xs,Ys)
        
        Hgrid = np.ones((len(Xs),len(Ys)))*np.min(meshH)
        
        xlist = Xgrid[0,:]
        ylist = Ygrid[:,0]
        
        xcoords = []
        ycoords = []
        
        for x,y,h in zip(meshX,meshY,meshH):
            
            xcoord = np.argwhere(vf.ismember(x,xlist))
            ycoord = np.argwhere(vf.ismember(y,ylist))
            
            xcoords=np.append(xcoords,xcoord)
            ycoords=np.append(ycoords,ycoord)
            
            Hgrid[xcoord,ycoord] = h
        
        o_field_y, o_field_x = np.gradient(Hgrid)
        
        meshOX = o_field_x[xcoords.astype(int),ycoords.astype(int)]
        meshOY = o_field_y[xcoords.astype(int),ycoords.astype(int)]
        
    
        return(Xgrid[xcoords.astype(int),ycoords.astype(int)],Ygrid[xcoords.astype(int),ycoords.astype(int)],meshOX,meshOY)
    
    def plot_splash_init(self,norm,**kwargs):
        
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
        
        meshX,meshY,meshOX,meshOY = self.orientation() # Cone config
        
        
        Theta,R = vf.ToCirc(meshX,meshY)
        
        inCone = self.Cone.IsIn(R, Theta)
        
        
        meshX = meshX[inCone]
        meshOX = meshOX[inCone]
        meshY = meshY[inCone]
        meshOY = meshOY[inCone]
        
        meshXci,meshYci = self.Cone.Cone2Circle(meshX, meshY) # Circle config
        meshOXci,meshOYci = self.Cone.Cone2Circle(meshOX, meshOY) # Circle config

         
        normCi = np.sqrt(np.square(meshOXci)+np.square(meshOYci))
        
        meshVXci = meshOXci/normCi*self.speednorm + self.speedtan
        meshVYci = meshOYci/normCi*self.speednorm

        fieldnormCi = np.sqrt(np.square(meshVXci)+np.square(meshVYci))
        
        meshVX, meshVY = self.Cone.Circle2Cone(meshVXci, meshVYci)

        fieldnorm = np.sqrt(np.square(meshVX)+np.square(meshVY))
        
        
        
        
        fig, ax = self.Cone.draw(drop=self.Drop,nolabels=NoLabels,dropview = 'full',conelinewidth=ConeLW,
                                 conecolor=ConeColor,title=Title,xlabelCi=Xlabel_Ci,ylabelCi=Ylabel_Ci,xlabelCo=Xlabel_Co,ylabelCo=Ylabel_Co)
        
        ax[0].quiver(meshXci, meshYci, -1*np.divide(meshVXci,fieldnormCi), -1*np.divide(meshVYci,fieldnormCi),fieldnormCi,scale = 15,zorder=5)
        ax[1].quiver(meshX, meshY, -1*np.divide(meshVX,fieldnorm), -1*np.divide(meshVY,fieldnorm),fieldnorm,scale = 15,zorder=5)

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    