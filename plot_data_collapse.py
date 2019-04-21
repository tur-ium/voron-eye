# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 13:55:21 2019

@author: Artur Donaldson
"""
from setup import im_directory,bp_directory,areas_directory
import numpy as np
import matplotlib.pyplot as plt
from voron_eye import plot_linear_distribution
from generalized_analysis import generalized_ad_test,generalized_ks_test
from scipy.stats import gamma,lognorm
from scipy.stats import ks_2samp
from estimate_area_errors import estimate_errs,plot_ci,plot_errorbar
# =============================================================================
# FORMATTING PARAMETERS
# =============================================================================
plt.rcParams['font.size'] = 14
formatting={"marker":"+","s":30} #Formatting for scatter plot for retina
xmin = 0
xmax = 5
ymin = 0
ymax = 1
n_formatting = dict(formatting)
n_formatting["color"]="#411900"
dr_formatting = dict(formatting)
dr_formatting["color"]="#ff0000"
pdr_formatting = dict(formatting)
pdr_formatting["marker"]="o"
h_formatting = dict(formatting)
h_formatting["color"]="#00ff00"
g_formatting = dict(formatting)
g_formatting["color"]="#0000ff"
# =============================================================================
# PARAMETERS
# =============================================================================
database = "HRF"
phi_N = 500
lat_N = 500
#Calculated parameters from author's thesis
c_params = {
        "Healthy":[3.504,0,0.285,0.012,0.85],  #H
        "Glaucomatous":[3.321,0,0.301,0.018,0.23],  #G
        "NPDR":[3.235,0,0.309,0.024,0.26],  #NPDR
        "PDR":[1.363,0,0.734,0.036,0.36]   #PDR

}
# =============================================================================
# LOAD DATA
# =============================================================================
if database == "HRF":
    HRF = True
    DRIVE = False
elif database=="DRIVE":
    HRF = False
    DRIVE = True

if HRF:
    R_im = int(.5*3190) #Radius of region of retina visible on image
    im_directory = im_directory+"HRF\\"
    dr_names =["{:02d}_dr".format(n) for n in range(1,16)]
    npdr_names =["{:02d}_dr".format(n) for n in [1,2,3,4,5,7,8,9,10,11,13,15]] 
    pdr_names = ["{:02d}_dr".format(n) for n in [6,12,14]]
    g_names = ["{:02d}_g".format(n) for n in range(1,16)]
    h_names = ["{:02d}_h".format(n) for n in range(1,16)]
    names = h_names+dr_names+g_names
if DRIVE:
    im_directory = im_directory+"DRIVE\\"
    dr_names = [3,8,14,17,25,26,32] #Subjects with DR
    #according to researchers https://drive.grand-challenge.org/
    names = list(range(1,41))
    h_names = list(names)
    for n in dr_names:
        h_names.remove(n)
    
#Directory containing branchpoint data and general info
areas_directory = areas_directory+database+"\\"
    

fa_sep = dict()
ca_sep = dict()

for name in names:
    print(name)
    if HRF:
        f = open(areas_directory+"curved_core_cell_areas_"+name+"_phiN={}latN={}.txt".format(phi_N,lat_N),"r")
    if DRIVE:
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

    if HRF:
        f = open(areas_directory+"flat_core_cell_areas_{}.txt".format(name),"r")
    if DRIVE:
        f = open(areas_directory+"flat_core_cell_areas_{}.txt".format(name),"r")
    areas = list()
    i = 0
    for line in f.readlines():
        if i > 1:
            areas = np.array([float(x) for x in line.rstrip().rsplit(",")[0:-1]])
        i+=1
    fa_sep[name]=areas/areas.mean()

#%%
# =============================================================================
# JOIN TOGETHER DISTRIBUTIONS
# =============================================================================
fa_h=list()
fa_dr=list()
fa_npdr=list()
fa_pdr=list()
fa_g=list()

ca_h=list()
ca_dr=list()
ca_npdr=list()
ca_pdr=list()
ca_g=list()

for n in h_names:
    fa_h.extend(fa_sep[n])
    ca_h.extend(ca_sep[n])
for n in dr_names:
    fa_dr.extend(fa_sep[n])
    ca_dr.extend(ca_sep[n])

if not DRIVE:
    for n in npdr_names:
        fa_npdr.extend(fa_sep[n])
        ca_npdr.extend(ca_sep[n])
    for n in pdr_names:
        fa_pdr.extend(fa_sep[n])
        ca_pdr.extend(ca_sep[n])
    for n in g_names:
        fa_g.extend(fa_sep[n])
        ca_g.extend(ca_sep[n])

#%%
# =============================================================================
# PLOT
# =============================================================================
plt.close("all")
x = np.linspace(xmin,xmax,1000)

fig=plt.figure()
ax=plt.subplot(111)

        
for condition in ["NPDR","PDR"]:#,"Glaucomatous"]:
    marker="+"
    if condition=="Healthy":
        indices = h_names
        fmt="g-"
    elif condition=="DR":
        indices = dr_names
        fmt="r--"
    elif condition=="NPDR":
        indices = npdr_names
        fmt="m-"    
        marker="P"
    elif condition=="PDR":
        indices = pdr_names
        fmt="r-"
        marker="o"
    elif condition=="Glaucomatous":
        indices = g_names
        fmt="b-"
    else:
        print("DONT KNOW THIS CONDITION")
    
    params = c_params[condition][:3]
    
    
    for name in indices:
        print(name)
        if len(ca_sep[name])>20:
            fig,ax = plot_linear_distribution(ca_sep[name],name="{}".format(name),n="auto",density=True,fig=fig,ax=ax,formatting={'marker':marker})
        else:
            print("Can't plot",name,"has",len(ca_sep[name]))
    #plot fit
    fit=gamma.pdf(x,*params)
    ax.plot(x,fit,fmt,label=condition)

    ax.legend(ncol=2)
    # =============================================================================
    #  FORMATTING    
    # =============================================================================
    fig.suptitle("{} {} Core Cell Areas".format(database,condition))
    ax.set_ylabel(r'Probability density')
    ax.set_xlabel(r"$ A/\langle A \rangle$")
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
for condition in ["Healthy","Glaucomatous"]:
    marker="+"
    fig=plt.figure()
    ax=fig.add_subplot(111)
    if condition=="Healthy":
        indices = h_names
        fmt="g-"
    elif condition=="DR":
        indices = dr_names
        fmt="r--"
    elif condition=="NPDR":
        indices = npdr_names
        fmt="m-"    
        marker="P"
    elif condition=="PDR":
        indices = pdr_names
        fmt="r-"
        marker="o"
    elif condition=="Glaucomatous":
        indices = g_names
        fmt="b-"
    else:
        print("DONT KNOW THIS CONDITION")
    
    params = c_params[condition][:3]
    
    
    for name in indices:
        print(name)
        if len(ca_sep[name])>20:
            fig,ax = plot_linear_distribution(ca_sep[name],name="{}".format(name),n="auto",density=True,fig=fig,ax=ax,formatting={'marker':marker})
        else:
            print("Can't plot",name,"has",len(ca_sep[name]))
    #plot fit
    fit=gamma.pdf(x,*params)
    ax.plot(x,fit,fmt,label=condition)

    ax.legend(ncol=2)
    # =============================================================================
    #  FORMATTING    
    # =============================================================================
    fig.suptitle("{} {} Core Cell Areas".format(database,condition))
    ax.set_ylabel(r'Probability density')
    ax.set_xlabel(r"$ A/\langle A \rangle$")
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
        