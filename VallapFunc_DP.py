# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 11:22:26 2021

@author: Valentin Laplaud

Useful generic functions 

"""


import numpy as np
import pandas as pd
import mpmath as mpm

import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.colors as colors

from scipy.spatial.distance import directed_hausdorff 
from scipy import optimize

from IPython import get_ipython

import tkinter as tk
from tkinter import ttk

# Intersection volume between sphere and cylinder

def interVolSC(Rs,Rc,Dsc):
    
                
    """ ref : Boersma and Kamminga, 1961. (https://core.ac.uk/download/pdf/82412251.pdf) """
    """ Equation (5) and (8) are implemented with the use of mpmath for computation of the elliptic integrals. """
    
    # Rs sphere radius, Rc cylinder radius, Dsc distance between the two centers
    
    # normalization to sphere radius
    rho = Rc/Rs
    eta = Dsc/Rs
    
    """ Formula is valid if sphere and cylinder intersects """
    if not (eta-rho) < 1:
        raise ValueError('Invalid case ! This code''s formula is only valid for intersecting sphere and cylinders \n -> (Dsc/Rs) - (Rc/Rs) < 1 !!')
        

    """ Heuman's lambda function """
    def Lambda0(beta,m):      
        """ From (https://link.springer.com/content/pdf/10.1007%2F978-3-642-65138-0.pdf) form 150.3 page 36 """
        L = 2/np.pi*(mpm.ellipe(m)*mpm.ellipf(beta,(1-m))+mpm.ellipk(m)*mpm.ellipe(beta,(1-m))-mpm.ellipk(m)*mpm.ellipf(beta,(1-m)))
        return(L)   
    
    
    if (eta+rho)>1 :
        """ Formula (5) is valid for eta+rho > 1 , this means that the part of the cylinder is outside the drop  """ 
        
        m = (1-(eta-rho)**2)/(4*rho*eta) # parameter for eliptic functions of mpmath (=k² in the paper)    
        theta = np.arcsin(eta-rho)  
        
        V = (2/3*np.pi*(1-Lambda0(theta,m)) 
             -8/9*np.sqrt(rho*eta)*(6*rho**2+2*rho*eta-3)*(1-m)*mpm.ellipk(m) 
             +8/9*np.sqrt(rho*eta)*(7*rho**2+eta**2-4)*mpm.ellipe(m))
    
    else:
        """ Formula (8) is valid for eta+rho <= 1 , this means that the cylinder is completely inside the drop  """ 
        
        m = (4*rho*eta)/(1-(eta-rho)**2) # parameter for eliptic functions of mpmath (=k² in the paper)    
        theta = np.arcsin((eta-rho)/(eta+rho))  
        
        V = (2/3*np.pi*(1-Lambda0(theta,m))
             -(4*np.sqrt(1-(eta-rho)**2))/(9*(eta+rho))*(2*rho-4*eta+(eta+rho)*(eta-rho)**2)*(1-m)*mpm.ellipk(m)
             +4/9*np.sqrt(1-(eta-rho)**2)*(7*rho**2+eta**2-4)*mpm.ellipe(m))

    
    return(float(V*Rs**3))



#  Compute normal to a vector

# Normal vector is normalized, and by default rotated counter clockwise from given vector
# the input vector is from (x1,y1) to (x2,y2)


def getNormal(x1,y1,x2,y2, **kwargs):
    
    rotation = 'CCW'
    
    for key, value in kwargs.items(): 
        if key == 'rotation':
            rotation = value
        else:            
            print('Unknown key : ' + key+ '. Kwarg ignored.')
            
    dx = x2 - x1
    dy = y2 - y1
    
    Norm = np.sqrt(np.square(dx) + np.square(dy))
    
    dxN = np.divide(dx,Norm)
    dyN = np.divide(dy,Norm)
    
    if rotation == 'CCW':
        
        x = -dyN
        y = dxN
        
    elif rotation == 'CW':
        
        x = dyN
        y = -dxN

    else:
        print('Wrong rotation parameter !! Should be ''CW'' or ''CCW''. Default is ''CCW'' with no input.')
            
    return(x,y) 
    
#  Plotting boxplots with data points on top 

# A function combining boxplot with seaborn's swarmplot for a better display of data

def boxswarmplot(Title,Ylabel,Data,facecolors,Labels,**kwargs):

    FS = (5,3)    

    for key, value in kwargs.items(): 
        if key == 'figsize':
            FS = value

    fig,ax = plt.subplots(dpi = 250,facecolor='white',figsize = FS)
    fig.suptitle(Title)
 
    
    cap= [None]*len(Data)
    med= [None]*len(Data)
    
    grouping = []
    
    for dat,col,lab,i in zip(Data,facecolors,Labels,range(len(Data))):
    
        # plots properties
        plotprops = {'color':'black'}
        boxprops = {'color':'black','facecolor':col}
        
        lab = lab + '\nn = ' + str(len(dat))
        
        Labels[i] = lab

        bp = ax.boxplot(dat, positions = [i], labels = [lab],patch_artist =True, boxprops=boxprops, capprops =plotprops,
                    showfliers=False,whiskerprops=plotprops,medianprops =plotprops)
        
        grouping = np.append(grouping,np.ones(len(dat))*i)
    
        cap[i] = bp['caps'][1].get_ydata(orig=True)[0]
        med[i] = bp['medians'][0].get_ydata(orig=True)[0]
    
    sns.swarmplot(x=grouping,y=pd.concat(Data),color = 'gray', size=2, ax = ax)
    
    ax.set_ylabel(Ylabel)
    
    ax.set_xticklabels(Labels)
    
    return(fig,ax,cap,med)
    

# Coordinate conversion from cartesian to circular (in deg) an vice versa

def ToCirc(X,Y, **kwargs):
    
    Angle = 'rad'
    
    for key, value in kwargs.items(): 
        if key == 'angle':
            Angle = value
        else:            
            print('Unknown key : ' + key+ '. Kwarg ignored.')
    
    
    if Angle == 'deg':
        Alpha = np.rad2deg(np.arctan2(Y,X))
    elif Angle == 'rad':
        Alpha = np.arctan2(Y,X)
    else:
        print('Wrong angle unit : ' + Angle + '. Default to radians.')         
        Alpha = np.arctan2(Y,X)
        
    Radius = np.sqrt(np.square(X)+np.square(Y))
    
    return(Alpha,Radius)



def ToCart(Alpha,Radius, **kwargs):
    
    Angle = 'rad'
    
    for key, value in kwargs.items(): 
        if key == 'angle':
            Angle = value
        else:            
            print('Unknown key : ' + key + '. Kwarg ignored.')
    
    if Angle == 'deg':
        Alpharad = np.deg2rad(Alpha)
    elif Angle == 'rad':
        Alpharad = Alpha
    else:
        print('Wrong angle unit : ' + Angle + '. Default to radians.') 
        Alpharad = Alpha
    
    
    X = Radius*np.cos(Alpharad)
    Y = Radius*np.sin(Alpharad)
    
    return(X,Y)

# Euclidian distance between two arrays of points in carthesian coordinates 
def dist(x1,y1,x2,y2):
    
    d = np.sqrt(np.square(x1-x2)+np.square(y1-y2))
    
    return(d)

# Computation of Hausdorff distance (https://en.wikipedia.org/wiki/Hausdorff_distance) between two contours
    
    
def HausdorffDist(x1,y1,x2,y2, **kwargs):
    
    DebugPlots = False

    for key, value in kwargs.items():
        if key == 'debug':
            DebugPlots = value
    
    c1 = [[x,y] for x,y in zip(x1,y1)]
    c2 = [[x,y] for x,y in zip(x2,y2)]
    
    d1, i11, i12 = directed_hausdorff(c1, c2)
    d2, i21, i22 = directed_hausdorff(c2, c1)
    D = max(d1,d2)
    
    
    if DebugPlots:   
        f0, ax0 = plt.subplots(dpi=200,facecolor='white')
        ax0.set_aspect('equal', adjustable='box')
        ax0.plot(x1,y1,'r.')
        ax0.plot(x2,y2,'b.')
        ax0.plot(x1[i11],y1[i11],'*g')
        ax0.plot(x2[i12],y2[i12],'*g')
        ax0.plot([x1[i11],x2[i12]],[y1[i11],y2[i12]],'g')
        ax0.plot(x1[i21],y1[i21],'*c')
        ax0.plot(x2[i22],y2[i22],'*c')
        ax0.plot([x1[i21],x2[i22]],[y1[i21],y2[i22]],'c')
         
   
    return(D)


# Simple ismember function, checks if A is within B
def ismember(A, B):
    return [ np.sum(b == A) for b in B ]


# R2 computation for a fit
def computeR2(Ydata,Yfit):
    # Ydata are the fitted data, Yfit the comuted value from the fit
    
    SumResidues = np.sum(np.square(np.subtract(Ydata,Yfit)))
    TotalVariance = np.sum(np.square(np.subtract(Ydata,np.mean(Ydata))))
    
    R2 = 1 - SumResidues/TotalVariance
    
    return R2
    


def dataSummary(GDs,Ns,labels,Mult,col,name,unit,method):
    DataPooled = np.empty(0)
    nPooled = np.sum(Ns)
    
    
    
    if method == 'mean' :
        print(name + ' (mean + STE) : ')
        
        for GD,n,lab in zip(GDs,Ns,labels):
            DataMean =  np.round(GD.loc[GD['Img'] == 0,col].mean()*Mult*100)/100
            Var =  np.round(GD.loc[GD['Img'] == 0,col].std()/np.sqrt(n)*Mult*100)/100
            DataPooled = np.append(DataPooled,GD.loc[GD['Img'] == 0,col].to_numpy())
            print(lab + ' -> ' + str(DataMean) + ' ' + unit + ' ' + u"\u00B1" + ' ' + str(Var)  + unit + ' (n = ' + str(n) + ')')
    
        PooledMean = np.round(np.mean(DataPooled)*Mult*100)/100
        PooledVar = np.round(np.std(DataPooled)/np.sqrt(np.sum(Ns))*Mult*100)/100
        print('Pooled -> ' + str(PooledMean) + ' ' + unit + ' ' + u"\u00B1" + ' ' + str(PooledVar)  + unit +  '  (n = ' + str(nPooled) + ')' )
        
    else:
        print(name + ' (median + AAD/med) : ')
        
        for GD,n,lab in zip(GDs,Ns,labels):
            DataMedian =  np.round(GD.loc[GD['Img'] == 0,col].median()*Mult*100)/100
            Var =  np.round(   np.mean(np.abs(GD.loc[GD['Img'] == 0,col]*Mult-DataMedian))*10000/DataMedian)/100
            DataPooled = np.append(DataPooled,GD.loc[GD['Img'] == 0,col].to_numpy())
            print(lab + ' -> ' + str(DataMedian) + ' ' + unit + ' ' + u"\u00B1" + ' ' + str(Var)  + ' % (n = ' + str(n) + ')')
    
        PooledMedian = np.round(np.median(DataPooled)*Mult*100)/100
        PooledVar = np.round(np.mean(np.abs(DataPooled*Mult - PooledMedian))*10000/PooledMedian)/100
        print('Pooled -> ' + str(PooledMedian) + ' ' + unit + ' ' + u"\u00B1" + ' ' + str(PooledVar)  + ' % (n = ' + str(nPooled) + ')' )
    
#  Fitting a circle to a collection of points  
def fitCircle(X,Y):
    
    x_m = np.mean(X)
    y_m = np.mean(Y)    
    
    def calc_R(xc, yc):

        return np.sqrt((X - xc) ** 2 + (Y - yc) ** 2)
    
    def f_2(c):
        
        Ri = calc_R(*c)
        return Ri - Ri.mean()
    
    center_estimate = x_m, y_m
    center_fit, _ = optimize.leastsq(f_2, center_estimate)

    xc_fit, yc_fit = center_fit
    Ri_fit       = calc_R(xc_fit, yc_fit)
    
     #Fitting the radius of the circle
    R_fit        = Ri_fit.mean()
   
    
    if False:
        fig,ax = plt.subplots(dpi=250, facecolor = 'black')
        ax.plot(X,Y,'ro')

        xcircle = xc_fit + R_fit*np.cos(np.linspace(-np.pi,np.pi,100))
        ycircle = yc_fit + R_fit*np.sin(np.linspace(-np.pi,np.pi,100))
        ax.plot(xcircle,ycircle,'k--')

    return([xc_fit, yc_fit,R_fit])
    
# Getting coordinates from user clicks on an image

def getContourPointsCoordinates(Img,Title):

    # ask user to click
    get_ipython().run_line_magic('matplotlib', 'qt')
    fig,ax = plt.subplots(dpi=250)
    ax.imshow(Img,cmap='gray')
    ax.set_xticks([])
    ax.set_yticks([]) 
    fig.suptitle(Title)    
    fig.tight_layout()
    pts = np.asarray(plt.ginput(n=-1, timeout=-1))
    plt.close()
    get_ipython().run_line_magic('matplotlib', 'inline')

    Xs = pts[:,0]
    Ys = pts[:,1]

    return(Xs,Ys)

# Creating a multiple choice button for the user


def button_choice(choices,questiontext):

    
    # selection window
    root = tk.Tk()
    root.geometry('300x220')
    root.resizable(False, False)
    root.title('User selection')
    
    # label
    label = ttk.Label(text=questiontext)
    label.pack(fill='x', padx=5, pady=5)
    
    # Variable for click result
    result = tk.StringVar()

    # what happens when the confirmation button is clicked
    def buttonAction():
        global EventChoice
        EventChoice = result.get()
        root.destroy()

   # choice buttons
    for choice in choices:
        r = tk.Radiobutton(
            root,
            text=choice[0],
            value=choice[1],
            variable=result,
        )
        r.pack(fill='x', padx=5, pady=5)

    # confirmation button
    button = tk.Button(
        root,
        text="Validate",
        command=buttonAction)
    
    button.pack(fill='x')

    root.mainloop()
    
    
    return(EventChoice)

# Truncate a colormap
def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap