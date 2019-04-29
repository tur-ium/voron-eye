# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 18:43:57 2019

@author: admin
"""

#%% VORONOI GRID ERRORS
# =============================================================================
# Estimate_area_errors
# =============================================================================
from scipy.ndimage.morphology import binary_erosion
import pickle as pkl
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
import numpy as np
from numpy import pi,cos,sin,sqrt
from voron_eye import optimum_bin_number
from setup import regions_directory, areas_directory

#%%
# =============================================================================
# ESTIMATE ERRORS HEALTHY
# =============================================================================
def estimate_errs(names,areas_directory,database):
    """Returns norm_areas_list,abs_area_errs_list"""
    phi_N = 500
    lat_N = 500
    epsilon = 0.0001 #Epsilon for floating point comparison
    
    save_data = False

    # =============================================================================
    # FORMATTING PARAMETERS
    # =============================================================================
    plt.rcParams['font.size'] = 14
    
    # =============================================================================
    # END OF PARAMETERS
    # =============================================================================
    
    if database == "HRF":
        HRF = True
        DRIVE = False
        
    elif database=="DRIVE":
        print("LOADING DRIVE")
        HRF = False
        DRIVE = True
    else:
        print("Database parameter can either be HRF or DRIVE")
    if HRF:
        R_im = int(.5*3190) #Radius of region of retina visible on image
        FOV = pi/3
        R_eye = R_im/cos(.5*pi-.5*FOV) #Actual Radius of eye
        im_directory = "C:\\Users\\admin\\Documents\\Physics\\year 4\\Vascular Networks\\iBranch\\"
#        dr_names = ["{:02d}_dr".format(n) for n in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]]
#        g_names = ["{:02d}_g".format(n) for n in range(1,16)]
#        h_names = ["{:02d}_h".format(n) for n in range(1,16)]
#        names = h_names+dr_names+g_names
#        
    if DRIVE:
        R_im = int(.5*536) #Radius of region of retina visible on image
        FOV = pi/4
        R_eye = R_im/cos(.5*pi-.5*FOV) #Actual Radius of eye
        im_directory = "C:\\Users\\admin\\Documents\\Physics\\year 4\\Vascular Networks\\data\\DRIVE\\"
#        dr_names = [3,8,14,17,25,26,32] #Subjects with DR
#        #according to researchers https://drive.grand-challenge.org/
#        names = list(range(1,41))
#        h_names = names
#        for n in dr_names:
#            h_names.remove(n)
    
    #%%
    lat_min = .5*pi-.5*FOV
    lat_max = .5*pi
    d_phi = 2*pi/phi_N
    d_lat = (lat_max-lat_min)/lat_N
    phis,lats = np.ogrid[-pi:pi:phi_N*1j,lat_min:lat_max:lat_N*1j]
    def area_on_eye(lam,R,d_phi,d_angle):
            return R**2*cos(lam)*d_phi*d_lat
    rel_errs = list()
    norm_areas_list = list()
    abs_area_errs_list=list()
    
    for i in range(len(names)):
        name = names[i]
        print(name)
        f = open(regions_directory+"{}_regions.pkl".format(name),"rb")
        regions = pkl.load(f)
        area_errs = list()
        f.close()
        
        #Read areas
        f = open(areas_directory+"curved_core_cell_areas_{}_phiN={}latN={}.txt".format(name,phi_N,lat_N),"r")
                
        i = 0
        for line in f.readlines():
            if i > 1:
                areas = np.array([float(x) for x in line.rsplit(",")[0:-1]])
            i+=1
        f.close()
        r_max = int(regions.max())
        
        print("Number of areas saved to file: {}".format(len(areas)))
        print("Maximum region index:          {}".format(regions.max()))
        prop_boundary = list()
    
        if r_max != areas.size:
            print("{} diff # of regions and areas".format(name))
        for r in range(r_max+1):
            region = (regions==r)
            region_area = areas[r]
            if region.sum()>0 and region_area>0:
    
                boundary = binary_erosion(region,iterations=1)^region
                bindices = boundary.nonzero()
                bcells = bindices[0].size
    
#                if bcells>.5*region.sum():
#                    print("Large boundary")
    #                plt.figure()
    #                plt.imshow(binary_erosion(region,iterations=1)^region)    
                
                barea = 0 #Area of gridsquares on boundary
                for z in range(bcells):
                    i1 = bindices[0][z]
                    i2 = bindices[1][z]
                
                    barea += area_on_eye(lats[0,i2],R_eye,d_phi,d_lat)
                
                
                prop_boundary.append(bcells/region.sum()) #Proportion of gridsquares in boundary
                rel_err = .5*barea/region_area
            else:
                rel_err = 0
            area_errs.append(rel_err)
        area_errs = np.array(area_errs)
        print("{} proportion of gridsquares on boundary: mean: {}, max: {}".format(name,np.mean(prop_boundary),max(prop_boundary)))
    
        if save_data:
            fwrite = open("voronoi_cell_areas\\{}\\curved_core_cell_area_errors_{}_phiN={}latN={}.txt".format(database,name,phi_N,lat_N),"w")
            fwrite.write("{},{}\n".format(name,len(regions)))
            fwrite.write(",".join(area_errs.astype(str)))
            fwrite.close()
        
        #Remove the errors corresponding to cells with zero area
        definedErrors = (area_errs != 0)
        rel_errs.extend(area_errs[definedErrors])
        
        #Normalize the areas
        nonzero_indices = areas.nonzero()
        areas_mean = areas[nonzero_indices].mean()
        norm_areas = areas/areas_mean
        norm_areas_list.extend(norm_areas)
        
        #Remove zero area errors
        abs_area_errs = area_errs*norm_areas[:r_max+1]
        abs_area_errs_list.extend(abs_area_errs)
        padding = len(norm_areas)-(r_max+1) #Padding due to difference between maximum region index and length of array
        abs_area_errs_list.extend([np.inf for x in range(padding)])
    return norm_areas_list,abs_area_errs_list
#%%
#mean_err=np.mean(rel_errs)
#print("Mean err {:.3f}".format(mean_err))
#plt.figure()
#plt.title("{}".format(database))
#vals,bins,patches = plt.hist(rel_errs,density=True,bins="auto")
#plt.ylabel("Probability Density")
#plt.xlabel("Relative error in area") #Number of gridpoints in boundary / number of gridpoints in each Voronoi region

# =============================================================================
# END OF VORONOI GRID ERRORS
# =============================================================================
mean_errs = {500:0.193,800:0.139} #x by x Gridsize, mean error
def pbin(l,u,s,m,N=1000):
    '''
    Probability Gaussian distributed random variable is in histogram bin
    l: lower bound
    u: upper bound
    s: standard deviation of Gaussian
    m: mean of Gaussian
    N: Number of points to calculate between l and u
    '''
    x = np.linspace(l,u,N)
    w = abs(u-l)/N
    pd = (s*sqrt(2*pi))**-1 * np.exp(-(x-m)**2/(2*s**2))
    return pd.sum()*w


# =============================================================================
# TEST HISTOGRAMMING
# =============================================================================
#from scipy.stats import gamma
#Ndata=100
#Nbins = 30
#x = gamma.rvs(2.5,0,0.5,Ndata)# [.5,1,2,3,3,3.5,3,4]
#err = np.random.random(size=Ndata)#[.15,1,1,1,1,1,1,1]
#
#bin_edges= np.linspace(0,10,Nbins+1)#[0,1,2,3,4,5]
#bin_centres = [sum(bin_edges[i:i+2])*.5 for i in range(len(bin_edges)-1)]
#masses = np.zeros(shape=Nbins)
#varns =  np.zeros(shape=Nbins) #Variances
#std_errs = np.zeros(shape=Nbins)
#for j in range(len(bin_centres)):
#    mass=0
#    varn=0
#    for i in range(len(x)):
#        l = bin_edges[j]
#        u = bin_edges[j+1]
#        p_i = pbin(l,u,err[i],x[i],N=1000)
#        mass += p_i
#        varn += p_i*(1-p_i)
#    masses[j] = mass
#    varns[j] = varn
#    std_errs[j] = varn**.5
#totmass = masses.sum()
#pdf = masses/totmass
#plt.figure()
#
#upper_ci = masses+std_errs
#lower_ci = masses-std_errs
##plt.errorbar(bin_centres,masses,yerr=std_errs.tolist())
#
##PLOT CONFIDENCE INTERVAL
#coords=np.zeros((Nbins*2,2))
#coords[:Nbins,0]=bin_centres
#coords[:Nbins,1]=lower_ci
#coords[Nbins:,0]=bin_centres[::-1]
#coords[Nbins:,1]=upper_ci[::-1]
#
#p = Polygon(coords,closed=True,alpha=.3,fill=(.1,.1,.1,.4),hatch=False)
#ax=plt.gca()
#ax.add_patch(p)
#
#plt.plot(bin_centres,upper_ci,'b--',label="$1 \sigma$")
#plt.plot(bin_centres,masses,label="Mean")
#plt.plot(bin_centres,lower_ci,'b--')
#plt.xlabel("x")
#plt.ylabel("Probability mass")
#plt.legend()

def plot_ci(x,err,density=True,fmt={'linecolor':'b-','fill':(.1,.1,.1,.4),'hatch':None},label="Mean",ax=None):
    """
    Plot distribution of a dataset x with errors err showing Confidence Intervals
    PARAMS
    ------
    x: np.array
        Data
    err: np.array
        Errors on data
    fmt: dict
        Formatting ('line':linestyle for plt.plot(), 'fill':fill style for matplotlib patch, 'hatch':hatches if True)
    density: boolean
        Plots pdf if True else plots mass in each bin
    ax: matplotlib axis object
        Plots on this axis, else if None creates an axis to plot on
    RETURNS
    -------
    fig,ax: matplotlib fig and ax
        matplotlib figure and axes of histogram
    """
    xmin = min(x)
    xmax = max(x)

    nonzero_indices=x.nonzero()
    Nbins=optimum_bin_number(x)
    bin_edges= np.linspace(xmin,xmax,Nbins+1)#[0,1,2,3,4,5]
    bin_centres = [sum(bin_edges[i:i+2])*.5 for i in range(len(bin_edges)-1)]
    bin_width = (xmax-xmin)/Nbins
    
    masses = np.zeros(shape=Nbins)
    varns =  np.zeros(shape=Nbins) #Variances
    std_errs = np.zeros(shape=Nbins)
    for j in range(len(bin_centres)):
        mass=0
        varn=0
        for i in nonzero_indices[0]:
            l = bin_edges[j]
            u = bin_edges[j+1]
            p_i = pbin(l,u,err[i],x[i],N=1000)
            mass += p_i
            varn += p_i*(1-p_i)
        masses[j] = mass
        varns[j] = varn
        std_errs[j] = varn**.5
    totmass = sum(masses)*bin_width
    pdf = masses/totmass
    std_errs_pdf = std_errs/totmass
    
# =============================================================================
#     PLOT
# =============================================================================
    if ax is None:
        fig=plt.figure()
        ax = fig.add_subplot(111)

    if density:
        #Plot PDF
        upper_ci = pdf+std_errs_pdf
        lower_ci = pdf-std_errs_pdf
        
        #PLOT CONFIDENCE INTERVAL
        coords=np.zeros((Nbins*2,2))
        coords[:Nbins,0]=bin_centres
        coords[:Nbins,1]=lower_ci
        coords[Nbins:,0]=bin_centres[::-1]
        coords[Nbins:,1]=upper_ci[::-1]
        
        p = Polygon(coords,closed=True,fill=fmt["fill"],hatch=fmt["hatch"])
        ax.add_patch(p)
        
        
        ax.plot(bin_centres,upper_ci,fmt["linecolor"]+"--",label="$1 \sigma$")
        ax.plot(bin_centres,pdf,fmt["linecolor"]+"-",label=label)
        ax.plot(bin_centres,lower_ci,fmt["linecolor"]+"--")
    else:
        upper_ci = masses+std_errs
        lower_ci = masses-std_errs
        
        #PLOT CONFIDENCE INTERVAL
        coords=np.zeros((Nbins*2,2))
        coords[:Nbins,0]=bin_centres
        coords[:Nbins,1]=lower_ci
        coords[Nbins:,0]=bin_centres[::-1]
        coords[Nbins:,1]=upper_ci[::-1]
        p = Polygon(coords,closed=True,fill=fmt["fill"],hatch=fmt["hatch"])
        ax.add_patch(p)
        
        ax.plot(bin_centres,upper_ci,fmt["linecolor"]+"--",label="$1 \sigma$")
        ax.plot(bin_centres,masses,fmt["linecolor"]+"-",label=label)
        
        ax.plot(bin_centres,lower_ci,fmt["linecolor"]+"--")
        ax.xlabel("x")
        ax.ylabel("Probability mass")
        ax.legend()
    fig=plt.gcf()
    return fig,ax

def plot_errorbar(x,err,density=True,color='black',fmt="+",label=None,ax=None,linestyle=None):
    """
    Plot distribution of a dataset x with errors err showing errorbars
    PARAMS
    ------
    x: np.array
        Data
    err: np.array
        Errors on data
    linestyle: str
        Line style for line between points
    color: str
        Color of error bars
    fmt: string
        Marker format
    density: boolean
        Plots pdf if True else plots mass in each bin
    ax: matplotlib axis object
        Plots on this axis, else if None creates an axis to plot on
    RETURNS
    -------
    fig,ax: matplotlib fig and ax
        matplotlib figure and axes of histogram
    """
    xmin = min(x)
    xmax = max(x)

    nonzero_indices=x.nonzero()
    Nbins=optimum_bin_number(x)
    bin_edges= np.linspace(xmin,xmax,Nbins+1)#[0,1,2,3,4,5]
    bin_centres = [sum(bin_edges[i:i+2])*.5 for i in range(len(bin_edges)-1)]
    bin_width = (xmax-xmin)/Nbins
    
    masses = np.zeros(shape=Nbins)
    varns =  np.zeros(shape=Nbins) #Variances
    std_errs = np.zeros(shape=Nbins)
    for j in range(len(bin_centres)):
        mass=0
        varn=0
        for i in nonzero_indices[0]:
            l = bin_edges[j]
            u = bin_edges[j+1]
            p_i = pbin(l,u,err[i],x[i],N=1000)
            mass += p_i
            varn += p_i*(1-p_i)
        masses[j] = mass
        varns[j] = varn
        std_errs[j] = varn**.5
    totmass = sum(masses)*bin_width
    pdf = masses/totmass
    std_errs_pdf = std_errs/totmass
    
# =============================================================================
#     PLOT
# =============================================================================
    if ax is None:
        fig=plt.figure()
        ax = fig.add_subplot(111)

    if density:
        #Plot PDF

        ax.errorbar(bin_centres,pdf,std_errs_pdf,linestyle=linestyle,color=color,fmt=fmt,capsize=5,label=label)

    fig=plt.gcf()
    return fig,ax