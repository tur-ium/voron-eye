# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 12:22:19 2019

@author: admin
"""
from scipy.stats import ks_2samp
import numpy as np
from estimate_area_errors import estimate_errs,plot_errorbar
from setup import areas_directory
# =============================================================================
# PARAMETERS
# =============================================================================

#Limits on axes
xmin = 0
xmax = 4
ymin = 0
ymax = 1

x=np.linspace(xmin,xmax,1000)
def load_data(database,areas_directory,names,phi_N=500,lat_N=500):
    """Load Voronoi cell area data. Returns two dictionaries containing normalized areas on each individual image.
The first dictionary contains cell areas loaded from curved data file the second loaded from the flat.
    PARAMS
    ----
    Database: 'HRF' or 'DRIVE'
    areas_directory: str
        Directory containing both the flat and curved core cell areas
    names: list
        Names of images to load
    phiN,latN: float,float
        Resolution of data to load.
    RETURNS
    -------
    ca_sep, fa_sep: dict, dict
        Curved areas, Flat areas on each individual image as described in the first line
"""
    #Directory containing branchpoint data and general info
    areas_directory = areas_directory+database+"\\"
    fa_sep = dict()
    ca_sep = dict()
    
    for name in names:
        print(name)
        f = open(areas_directory+"curved_core_cell_areas_{}_phiN={}latN={}.txt".format(name,phi_N,lat_N),"r")
        areas = list()
        i = 0
        for line in f.readlines():
            if i > 1:
                areas = [float(x) for x in line.rsplit(",")[0:-1]]
            i+=1
        areas = np.array(areas)
        #Remove zero areas
        nonzero_areas = areas[areas.nonzero()]
        #Normalize
        ca_sep[name]=nonzero_areas/nonzero_areas.mean()
    
        f = open(areas_directory+"flat_core_cell_areas_{}.txt".format(name),"r")
        
        areas = list()
        i = 0
        for line in f.readlines():
            if i > 1:
                areas = np.array([float(x) for x in line.rstrip().rsplit(",")[0:-1]])
            i+=1
        fa_sep[name]=areas/areas.mean()
    return ca_sep,fa_sep
#%%
# =============================================================================
# LOAD HRF
# =============================================================================
R_im = int(.5*3190) #Radius of region of retina visible on image
im_directory = "im\\HRF"
dr_names =["{:02d}_dr".format(n) for n in range(1,16)]
npdr_names =["{:02d}_dr".format(n) for n in [1,2,3,4,5,7,8,9,10,11,13,15]] 
pdr_names = ["{:02d}_dr".format(n) for n in [6,12,14]]
g_names = ["{:02d}_g".format(n) for n in range(1,16)]
h_names = ["{:02d}_h".format(n) for n in range(1,16)]
names = h_names+dr_names+g_names
caHRF,faHRF = load_data("HRF",areas_directory,names)


ca_hHRF=list()
ca_drHRF=list()
ca_hDRIVE=list()
ca_drDRIVE=list()

for n in h_names:
    ca_hHRF.extend(caHRF[n])
for n in npdr_names:
    ca_drHRF.extend(caHRF[n])
# =============================================================================
# PLOT HRF
# =============================================================================
ca_sep = caHRF
c_params = {
        "Healthy":[3.504,0,0.285,0.012,0.85],  #H
        "NPDR":[3.235,0,0.309,0.024,0.26],  #NPDR

}
marker="o"
fig=plt.figure()
ax=plt.subplot(111)
for condition in ["Healthy","NPDR"]:
    if condition=="Healthy":
        indices = h_names
        fmt={"color":"g"}
    elif condition=="NPDR":
        indices = npdr_names
        fmt={"color":"r"}
        
    else:
        print("DONT KNOW THIS CONDITION")
    
    params = c_params[condition][:3]#2.334,0,0.428
    
    
#    
#    for name in indices:
#        print(name)
#        if len(ca_sep[name])>20:
#            fig,ax = plot_linear_distribution(ca_sep[name],name="{}".format(name),n="auto",density=True,fig=fig,ax=ax)
#        else:
#            print("Can't plot",name,"has",len(ca_sep[name]))
    #plot fit
    fit=gamma.pdf(x,*params)
    ax.plot(x,fit,fmt["color"]+"-",label=condition+" HRF fit")

    
    #Plot histogram with confidence interval
    print("#########################")
    a,err=estimate_errs(indices,areas_directory+"HRF\\","HRF")
    fig,ax=plot_errorbar(np.array(a),np.array(err),label=condition+" HRF",ax=ax,color=fmt["color"],fmt=marker)
#    a,err=estimate_errs(indices,areas_directory+"HRF\\","HRF")
#    fig,ax=plot_errorbar(np.array(a),np.array(err),fmt={'linecolor':fmt[0],'fill':(.1,.1,.1,0),'hatch':None},label=condition,ax=ax)
    #ax.legend()
# =============================================================================
# LOAD DRIVE
# =============================================================================
c_params = {
        "Healthy":[3.151,0,0.317,0.012,0.85],  #H
        "NPDR":[2.692,0,0.371,0.024,0.26],  #NPDR

}
npdr_names = [3,8,14,17,25,26,32] #Subjects with DR
#according to researchers https://drive.grand-challenge.org/
names = list(range(1,41))
h_names = list(names)
for n in npdr_names:
    h_names.remove(n)
caDRIVE,faDRIVE = load_data("DRIVE",areas_directory,names)
for n in h_names:
    ca_hDRIVE.extend(caDRIVE[n])
for n in npdr_names:
    ca_drDRIVE.extend(caDRIVE[n])
    
print("DRIVE vs HRF HEALTHY KS TEST")
print(ks_2samp(ca_hDRIVE,ca_hHRF))
print("DRIVE vs HRF KS NPDR TEST")
print(ks_2samp(ca_drDRIVE,ca_drHRF))
#generalized_ks_test(ca_drDRIVE,ca_drHRF,floc=True,label="NPDR DRIVE HRF comparison")

# =============================================================================
# PLOT
# =============================================================================
ca_sep = caDRIVE
marker="+"
for condition in ["Healthy","NPDR"]:
    if condition=="Healthy":
        indices = h_names
        fmt={"color":"g"}
    elif condition=="NPDR":
        indices = npdr_names
        fmt={"color":"r"}
        
    else:
        print("DONT KNOW THIS CONDITION")
    
    params = c_params[condition][:3]#2.334,0,0.428
    
    
    for name in indices:
        print(name)
#        if len(ca_sep[name])>20:
#            fig,ax = plot_linear_distribution(ca_sep[name],name="{}".format(name),n="auto",density=True,fig=fig,ax=ax)
##            #fig,ax=plot_linear_distribution(ca_sep[name],"{}".format(name),n="auto",density=True,fig=fig,ax=ax,formatting=formatting)
#        else:
#            print("Can't plot",name,"has",len(ca_sep[name]))
    #plot fit
    fit=gamma.pdf(x,*params)
    ax.plot(x,fit,fmt["color"]+"--",label=condition+" DRIVE fit")
#    fit2=gamma.pdf(x,*DRIVEparams)
#    ax.plot(x,fit2,label="Fit DRIVE")
    #Plot histogram with confidence interval
    a,err=estimate_errs(indices,areas_directory+"DRIVE\\","DRIVE")
    fig,ax=plot_errorbar(np.array(a),np.array(err),label=condition+" DRIVE",ax=ax,color=fmt["color"],fmt=marker)
    #fig,ax=plot_errorbar(np.array(a),np.array(err),fmt={'linecolor':fmt[0],'fill':(.1,.1,.1,0),'hatch':None},label=condition,ax=ax)
    ax.legend()
# =============================================================================
#  FORMATTING    
# =============================================================================
fig.suptitle("{} Core Cell Areas".format(condition))
ax.set_ylabel(r'Probability density')
ax.set_xlabel(r"$ A/\langle A \rangle$")
ax.set_xlim(xmin,xmax)
ax.set_ylim(ymin,ymax)

