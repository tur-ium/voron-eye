# -*- coding: utf-8 -*-
"""
Created on Tue Apr  9 12:51:45 2019

@author: admin
"""
from setup import areas_directory,im_directory
from estimate_area_errors import estimate_errs,plot_errorbar
from voron_eye import optimum_bin_number,plot_linear_distribution
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gamma,lognorm
from generalized_analysis import generalized_ks_test, generalized_ad_test
# =============================================================================
# FORMATTING PARAMETERS
# =============================================================================
plt.rcParams['font.size'] = 14
formatting={"marker":"+","s":20,"color":"red"} #Formatting for scatter plot for retina
#Format for Neovascularization
n_formatting = formatting.copy()
n_formatting["color"]="#411900"
# =============================================================================
# PARAMETERS
# =============================================================================
phi_N = 500
lat_N = 500
database = "HRF"
names = ["{:02d}_dr".format(n) for n in [6,12,14]]
xmax = 5
floc=True
save_data =  True
pop_params = 3.235,0,0.309
pdr_params = 1.363,0,0.734
pdr_names = ["{:02d}_dr".format(n) for n in range(1,16)]

# =============================================================================
# END OF PARAMETERS
# =============================================================================

def remove_zeros(x):
    y = x.copy()
    while y.count(0)>0:
        y.remove(0)
    return y

flat_core_cell_areas_sep = dict()
curved_core_cell_areas_sep = dict()


## =============================================================================
## JUST NEOVASCULARIZED REGIONS
## =============================================================================
database = "HRF\\neovascularized\\"
flat_core_cell_areas_sep = dict()
curved_core_cell_areas_sep = dict()

areas_directory = areas_directory+database
for name in names:
    f = open(areas_directory+"curved_core_cell_areas_"+name+"_phiN={}latN={}.txt".format(phi_N,lat_N),"r")
    
    areas = list()
    i = 0
    for line in f.readlines():
        if i > 1:
            areas = np.array([float(x) for x in line.rsplit(",")[0:-1]])
        i+=1
    areas = np.array(areas)
    mean = areas[areas.nonzero()].mean()
    curved_core_cell_areas_sep[name]=areas/mean

curved_neovascularized_areas = list()
for name in names:
    curved_neovascularized_areas.extend(curved_core_cell_areas_sep[name]/curved_core_cell_areas_sep[name].mean())
curved_neovascularized_areas = remove_zeros(curved_neovascularized_areas)

# =============================================================================
# PRINT GENERAL INFO
# =============================================================================
print("Number of non-zero Voronoi cells {}".format(len(curved_neovascularized_areas)))
print("Average cells per image {}".format(len(curved_neovascularized_areas)/len(names)))

# =============================================================================
# LOAD POPULATION DR AREAS
# =============================================================================    

curved_core_cell_areas_sep = dict()

for name in names:
    f = open(areas_directory+"curved_core_cell_areas_"+name+"_phiN={}latN={}.txt".format(phi_N,lat_N),"r")
    areas = list()
    i = 0
    for line in f.readlines():
        if i > 1:
            areas = np.array([float(x) for x in line.rsplit(",")[0:-1]])
        i+=1
    nonzero_areas = areas[areas.nonzero()]
    curved_core_cell_areas_sep[name]=nonzero_areas/nonzero_areas.mean()

pop_areas = list()
for n in names:
    pop_areas.extend(curved_core_cell_areas_sep[n])
pop_areas=remove_zeros(pop_areas)


# =============================================================================
# FITS
# =============================================================================
x = np.linspace(0,xmax,1000)
fit_pop = gamma.pdf(x,*pdr_params)

if floc:
    neo_params = gamma.fit(curved_neovascularized_areas,floc=0)
fit_neo = gamma.pdf(x,*neo_params)
# =============================================================================
# PLOT
# =============================================================================
fig = plt.figure()
ax = plt.subplot(111)
#plot_linear_distribution(pop_areas,"DR population", density=True,n="auto",fig=fig,ax=ax)

plt.plot(x,fit_pop,label="PDR whole image",c="red")
plt.plot(x,fit_neo,label="Neovascularized region only",c="#411900")

plot_linear_distribution(pop_areas,"PDR", density=True,n="auto",fig=fig,ax=ax,formatting=formatting)
plot_linear_distribution(curved_neovascularized_areas,"Neovascularized", density=True,n="auto",fig=fig,ax=ax,formatting=n_formatting)

fig.suptitle("Core cell area distribution (DR)")
ax.set_xlabel(r"$A/\langle A \rangle$")
ax.set_xlim([0,xmax])
ax.set_ylim([0,.8])
#%%
database = "Neovascularized"
if save_data:
    datasets = [curved_neovascularized_areas]
    dataset_names = ["neovascularized"]
    fitfuncs = [gamma,lognorm]
    fitfunc_names=["Gamma","Lognorm"]
    fwrite = open("{}_fit_results.txt".format(database),"w")
    for d in range(len(datasets)):
        fwrite.write("{},shape,loc,scale,k value,p value,ad k value,ad sig level,Ncells\n".format(dataset_names[d]))
        areas_norm = np.array(datasets[d])
        N = len(datasets[d])
        
        ax.set_xlabel(r"$A/\langle A\rangle$")
        for i in range(len(fitfuncs)):
            k, p, k_err, p_err = generalized_ks_test(areas_norm,fitfuncs[i],label=dataset_names[d]+" "+fitfunc_names[i],floc=floc,ax=ax)
            ad_k, siglevel, stat_err, ad_k_err = generalized_ad_test(areas_norm,fitfuncs[i],label=dataset_names[d]+fitfunc_names[i],floc=floc)
            if floc:
                params = fitfuncs[i].fit(areas_norm,floc=0)
            else:
                params = fitfuncs[i].fit(areas_norm)
            formatters = ["{:.3f}" for x in range(len(params))]
            fwrite.write("{},".format(fitfunc_names[i])+",".join(formatters).format(*params)+",{:.3f},{:.5f},{:.3f},{:5f},{:03d}\n".format(k,p,ad_k,siglevel, N))
        ax.legend()
    fwrite.close()