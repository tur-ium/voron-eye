# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 20:28:57 2019

@author: Artur Donaldson

A simple model for the effects of neovascularization of the Voronoi cell area
diagram using blood vessel branchpoints as seed points. 

The simulation is run NR times to form an ensemble of realizations of the blood
vessel network. Please be patient or adjust NR.

Plots the Voronoi cell diagrams at a series of observation points and the \
peak of the normalized area distribution over time using a moving average with a\
window width of 5

"""

from scipy.spatial import Voronoi, voronoi_plot_2d
import numpy as np
from scipy.stats import gamma
from voron_eye import compare_flat_and_curved_distributions
from shapely.geometry import Polygon, Point
import matplotlib.patches as pltpatches
from matplotlib import pyplot as plt
from generalized_analysis import generalized_ks_test
#%%
# =============================================================================
# PARAMETERS
# =============================================================================
im_w = 3500
im_h = 2336
NR = 10 #Number of realizations of the Voronoi cells to include in the ensemble
mask = Polygon([[0,0],[im_w,0],[im_w,im_h],[0,im_h]])#REGION OF THE IMAGE TO INCLUDE
Nsites = 100 #Initial number of sites 
addThisMany = 150 #Number of sites to add

epsilon = 1e-10

plot = True #Plot Voronoi cell diagrams
observation_points = [0,49,99,149] #Plot them after adding these number of sites

np.random.seed(3142) #Seed for random number generator
# =============================================================================
# END OF PARAMETERS
# =============================================================================
# =============================================================================
# FORMATTING PARAMETERS
# =============================================================================
plt.rcParams['font.size'] = 14
formatting={"marker":"+","s":20} #Formatting for scatter plot for retina
# =============================================================================
# END OF FORMATTING PARAMETERS
# =============================================================================
#plt.close("all")
#%%
if plot:
    fig1 = plt.figure() #Voronoi cells 

fig2 = plt.figure()
ax2 = plt.subplot(111)

sum_at_time = np.zeros(addThisMany)
sum_sq_at_time = np.zeros(addThisMany)
ensemble = np.array((NR,addThisMany))

for plots in range(NR):
    sites = np.random.random_sample(size=(Nsites,2))*np.array([[im_w,im_h]])
    print(plots)
    #ax = plt.subplot(111)
    #plt.scatter(sites[:,0],sites[:,1])
    #V = Voronoi(sites,incremental=True)
    #voronoi_plot_2d(V,ax=ax,points=False,vertices=False)
    unused,areas = compare_flat_and_curved_distributions(sites,mask,1,1,im_w,im_h,1,1,1,0,False,threshold=epsilon)
    # =============================================================================
    # ADD SITES
    # =============================================================================
    mean_log = list()
    mode_log = list()
    median_log = list()
    areas_log = list()
    
    observations = 0
    xmin = 0
    xmax = max(areas/areas.mean())
    for x in range(addThisMany):
        
        sites = np.append(sites,np.random.random(size=(1,2))*.5*np.array([[im_w,im_h]]),axis=0)
        areas, unused = compare_flat_and_curved_distributions(sites,mask,1,1,im_w,im_h,1,1,1,0,False,threshold=epsilon)
        
        
        mean = np.mean(areas)
        norm_areas = areas/mean
        h,bins = np.histogram(norm_areas,density=True,bins="auto")
        #METHOD 1 FITTING (BEWARE MAY NOT BE A GOOD FIT TO GAMMA!)
#        a,loc,scale = gamma.fit(norm_areas,floc=0)
#        print("Parameters for best gamma fit to (shape,loc,scale): {:.3f}, {:.3f}, {:.3f}".format(a, loc, scale))
#        fit_h = gamma.pdf(np.linspace(xmin,xmax,100),a,loc,scale)
        #mode_log.append(mode(norm_areas)[0])
#        mode_log.append((a-1)*scale)
        
        #METHOD 2 histogramming rather than fitting
        mode = bins[h.argmax():h.argmax()+1]*.5
        tot_area = sum(areas)
        mean_log.append(mean)
        mode_log.append(mode)
        median_log.append(np.median(norm_areas))
        areas_log.append(norm_areas)
        sum_at_time[x]+=mode_log[x]
        sum_sq_at_time[x]+=mode_log[x]**2
        
        if plot and x in observation_points and plots==0:
            observations+=1
            a = int(.5*len(observation_points))*100+int(.5*len(observation_points))*10+observations
            ax = fig1.add_subplot(a)
            plt.figure()
            voronoi_plot_2d(Voronoi(sites),show_points=None,show_vertices=False,ax=ax)
            ax.set_xlim([0,im_w])
            ax.set_ylim([0,im_h])
            ax.set_title("{} points".format(x+1+Nsites))
            ax.axis("off")
            #generalized_ks_test(norm_areas,gamma,floc=True,label="vasculature after {} pathological branchpoints".format(x),reps=10,N=int(1e7))
#%%
#    ax2.plot(mean_log,mode_log)
#    xmin = min(mean_log)
#    xmax = max(mean_log)*1.1
#    ax2.set_xlim((0,max(mean_log)))
#    label_xoffset = (xmax-xmin)/100
#    for x in observation_points:
#        ax2.scatter(mean_log[x],mode_log[x])
#        ax2.text(mean_log[x]-label_xoffset,mode_log[x]+.01,"{}".format(x+1+Nsites))
#    ax2.set_title("Effect on mode as sites are added (number shows number of sites)")
#    ax2.set_ylabel(r"Mode $A/\langle A \rangle$")
#    ax2.set_xlabel(r"Mean area")

#Average data
ax2.clear()
peaks = sum_at_time/NR
std = np.sqrt(sum_sq_at_time/NR-peaks**2)/np.sqrt(NR)
#Smooth data by averaging over 5 datapoints
smoothed = peaks.reshape((30,5))
means = smoothed.mean(axis=1)
#errs = smoothed.std(axis=1)
std_reshaped = std.reshape((30,5))
errs = std_reshaped.mean(axis=1)/np.sqrt(5)

#Plot
#ax2.errorbar(range(addThisMany)[::5],means,yerr=errs)
def plotCI(x,y,y_errs,ax,label,formatting={'alpha':.3,'fill':(.1,.1,.1,.4),'hatch':False}):
    Nx = len(x)
    upper_ci = y+y_errs
    lower_ci = y-y_errs
    #plt.errorbar(bin_centres,masses,yerr=std_errs.tolist())
    
    #PLOT CONFIDENCE INTERVAL
    coords=np.zeros((len(x)*2,2))
    coords[:Nx,0]=x
    coords[:Nx,1]=lower_ci
    coords[Nx:,0]=x[::-1]
    coords[Nx:,1]=upper_ci[::-1]

    p = pltpatches.Polygon(coords,closed=True,**formatting)
    ax.add_patch(p)

    ax.plot(x,upper_ci,'b--',label="$1 \sigma$")
    ax.plot(x,lower_ci,'b--')
    ax.plot(x,y,label=label)
x = range(addThisMany)[::5]
plotCI(x,means,errs,ax2,"Moving average")
fig2.suptitle("Neovascularization simulation, $N_{init}="+str(Nsites)+", N=$"+str(NR))
ax2.set_ylabel(r"$(A/\langle A\rangle)_{peak}$")
ax2.set_xlabel("Pathological branchpoints added, $t$")
ax2.set_xlim([0,addThisMany])

