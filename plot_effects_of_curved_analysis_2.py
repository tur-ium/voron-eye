# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 12:34:47 2019

@author: admin
"""
from setup import areas_directory,im_directory
import numpy as np
import matplotlib.pyplot as plt
from voron_eye import plot_linear_distribution
from generalized_analysis import generalized_ad_test,generalized_ks_test
from scipy.stats import gamma,lognorm
from scipy.stats import ks_2samp
from estimate_area_errors import plot_errorbar,estimate_errs
# =============================================================================
# FORMATTING PARAMETERS
# =============================================================================
plt.rcParams['font.size'] = 14
formatting={"marker":"+","s":20} #Formatting for scatter plot for retina
fmt_h,fmt_g,fmt_dr = formatting.copy(),formatting.copy(),formatting.copy()
fmt_h["color"]="green"
fmt_dr["color"]="red"
fmt_g["color"]="blue"
# =============================================================================
# PARAMETERS
# =============================================================================
plot_gamma = True
plot_lognorm = False
floc = True
phi_N = 500
lat_N = 500

savedata = False
#Limits on axes
xmin = 0
xmax = 4
ymin = 0
ymax = 1

#HRF or DRIVE
database="HRF"

f_params = [
        [3.657,0.000,0.273,0.010,0.95], #H
        [3.442,0.000,0.291,0.016,0.41], #G
        [2.732,0.000,0.366,0.021,0.24]  #DR
]
c_params = [
        [3.504,0,0.285,0.012,0.85],  #H
        [3.321,0,0.301,0.018,0.23],  #G
        [3.235,0,0.309,0.024,0.26],  #NPDR
        [1.363,0,0.734,0.036,0.36]   #PDR

]

# =============================================================================
# END OF PARAMETERS
# =============================================================================
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

def remove_zeros(x):
    y = x.copy()
    while y.count(0)>0:
        y.remove(0)
    return y

ca_dr=remove_zeros(ca_dr)
ca_g=remove_zeros(ca_g)
ca_h=remove_zeros(ca_h)
#%%
# =============================================================================
# PLOT
# =============================================================================
#plt.close("all")
fig5 = plt.figure()
ax5 = fig5.add_subplot(111)

# =============================================================================
# HEALTHY
# =============================================================================
print("---\nHEALTHY\n---")
fmt=fmt_h.copy()

#GAMMA
data = fa_h
condition = "Healthy flat"
fmt.update({"marker":"x"})
if plot_gamma:
    a,loc,scale = gamma.fit(data,floc=0)
    generalized_ks_test(data,gamma,label="Gamma",floc=floc)
    fit_h = gamma.pdf(np.linspace(xmin,xmax,100),a,loc,scale)
    ax5.plot(np.linspace(xmin,xmax,100),fit_h,"g--",label="Gamma fit for {}".format(condition))
if plot_lognorm:
    #LOGNORM
    if floc: 
        a,loc,scale = lognorm.fit(data,floc=0) 
    else: 
        a,loc,scale = lognorm.fit(data)
    generalized_ks_test(data,lognorm,label="Lognorm",floc=floc)
    fit_h = lognorm.pdf(np.linspace(xmin,xmax,100),a,loc,scale)
    ax5.plot(np.linspace(xmin,xmax,100),fit_h,"g--",label="fit for {}".format(condition))
plot_linear_distribution(data,condition,n="auto",density=True,fig=fig5,ax=ax5,formatting=fmt)

data = ca_h
condition = "Healthy curved"
fmt.update({"marker":"o"})
if plot_gamma:
    a,loc,scale = gamma.fit(data,floc=0)
    generalized_ks_test(data,gamma,label="Gamma",floc=floc)
    fit_h = gamma.pdf(np.linspace(xmin,xmax,100),a,loc,scale)
    ax5.plot(np.linspace(xmin,xmax,100),fit_h,"g-",label="Gamma fit for {}".format(condition))
if plot_lognorm:
    #LOGNORM
    if floc: 
        a,loc,scale = lognorm.fit(data,floc=0) 
    else: 
        a,loc,scale = lognorm.fit(data)
    generalized_ks_test(data,lognorm,label="Lognormal",floc=floc)
    fit = lognorm.pdf(np.linspace(xmin,xmax,100),a,loc,scale)
    ax5.plot(np.linspace(xmin,xmax,100),fit,"g-",label="fit for {}".format(condition))
#plot_linear_distribution(data,condition,n="auto",density=True,fig=fig5,ax=ax5,formatting=fmt)
norm_areas, errs= estimate_errs(h_names,areas_directory,database="HRF")
plot_errorbar(np.array(norm_areas),np.array(errs),label=condition,ax=ax5,color=fmt["color"])

#%%
# =============================================================================
# DR
# =============================================================================
print("---\nDR\n---")
fmt=fmt_dr
data = fa_pdr
condition = "PDR flat"
fmt.update({"marker":"x"})
if plot_gamma:
    a,loc,scale = gamma.fit(data,floc=0)
    generalized_ks_test(data,gamma,label="Gamma",floc=floc)
    fit_h = gamma.pdf(np.linspace(xmin,xmax,100),a,loc,scale)
    ax5.plot(np.linspace(xmin,xmax,100),fit_h,"r--",label="Gamma fit for {}".format(condition))
if plot_lognorm:
    #LOGNORM
    if floc: 
        a,loc,scale = lognorm.fit(data,floc=0) 
    else:
        a,loc,scale = lognorm.fit(data)
    generalized_ks_test(data,lognorm,label="Lognormal",floc=floc)
    fit_dr = lognorm.pdf(np.linspace(xmin,xmax,100),a,loc,scale)
    ax5.plot(np.linspace(xmin,xmax,100),fit_dr,"r--",label="fit for {}".format(condition))
plot_linear_distribution(data,condition,n="auto",density=True,fig=fig5,ax=ax5,formatting=fmt)

data = ca_pdr
condition = "PDR curved"
fmt.update({"marker":"o"})
if plot_gamma:
    a,loc,scale = gamma.fit(data,floc=0)
    generalized_ks_test(data,gamma,label="Gamma",floc=floc)
    fit_h = gamma.pdf(np.linspace(xmin,xmax,100),a,loc,scale)
    ax5.plot(np.linspace(xmin,xmax,100),fit_h,"r-",label="Gamma fit for {}".format(condition))
if plot_lognorm:
    #LOGNORM
    if floc: 
        a,loc,scale = lognorm.fit(data,floc=0) 
    else:
        a,loc,scale = lognorm.fit(data)
    generalized_ks_test(data,lognorm,label="Lognormal",floc=floc)
    fit = lognorm.pdf(np.linspace(xmin,xmax,100),a,loc,scale)
    ax5.plot(np.linspace(xmin,xmax,100),fit,"r-",label="fit for {}".format(condition))
#plot_linear_distribution(data,condition,n="auto",density=True,fig=fig5,ax=ax5,formatting=fmt)
norm_areas, errs= estimate_errs(pdr_names,areas_directory,database="HRF")
plot_errorbar(np.array(norm_areas),np.array(errs),label=condition,ax=ax5,color=fmt["color"])

#%%
if not DRIVE:
    # =============================================================================
    # GLAUCOMATOUS
    # =============================================================================
    print("---\nGLAUCOMATOUS\n---")
    fmt=fmt_g
    data = fa_g
    condition = "Glaucomatous flat"
    fmt.update({"marker":"x"})
    if plot_gamma:
        a,loc,scale = gamma.fit(data,floc=0)
        generalized_ks_test(data,gamma,label="Gamma",floc=floc)
        fit_h = gamma.pdf(np.linspace(xmin,xmax,100),a,loc,scale)
        ax5.plot(np.linspace(xmin,xmax,100),fit_h,"b--",label="Gamma fit for {}".format(condition))
    if plot_lognorm:
        #LOGNORM
        if floc: 
            a,loc,scale = lognorm.fit(data,floc=0) 
        else: 
            a,loc,scale = lognorm.fit(data)
        generalized_ks_test(data,lognorm,label="Lognormal",floc=floc)
        fit_g = lognorm.pdf(np.linspace(xmin,xmax,100),a,loc,scale)
        ax5.plot(np.linspace(xmin,xmax,100),fit_g,"b--",label="fit for {}".format(condition))
    plot_linear_distribution(data,condition,n="auto",density=True,fig=fig5,ax=ax5,formatting=fmt)
    
    data = ca_g
    condition = "Glaucomatous curved"
    fmt.update({"marker":"o"})
    if plot_gamma:
        a,loc,scale = gamma.fit(data,floc=0)
        generalized_ks_test(data,gamma,label="Gamma",floc=floc)
        fit_h = gamma.pdf(np.linspace(xmin,xmax,100),a,loc,scale)
        ax5.plot(np.linspace(xmin,xmax,100),fit_h,"b-",label="Gamma fit for {}".format(condition))
    if plot_lognorm:
        #LOGNORM
        if floc: 
            a,loc,scale = lognorm.fit(data,floc=0) 
        else:
            a,loc,scale = lognorm.fit(data)
        generalized_ks_test(data,lognorm,label="Lognormal",floc=floc)
        fit_g = lognorm.pdf(np.linspace(xmin,xmax,100),a,loc,scale)
        ax5.plot(np.linspace(xmin,xmax,100),fit_g,"b-",label="fit for {}".format(condition))
    #plot_linear_distribution(data,condition,n="auto",density=True,fig=fig5,ax=ax5,formatting=fmt_g)
    norm_areas, errs= estimate_errs(g_names,areas_directory,database="HRF")
    plot_errorbar(np.array(norm_areas),np.array(errs),label=condition,ax=ax5,color=fmt["color"])
#%% 
# =============================================================================
# FORMATTING
# =============================================================================
fig5.suptitle("HRF Curved versus Flat Analysis")
ax5.set_ylabel(r'Probability density')
ax5.set_xlabel(r"$ A/\langle A \rangle$")
ax5.set_xlim(xmin,xmax)
ax5.set_ylim(ymin,ymax)
ax5.legend()

#%%
# =============================================================================
# SAVE DATA
# =============================================================================
#from scipy.stats import alpha,anglit,arcsine,beta,betaprime,bradford,burr,cauchy,chi,chi2,cosine,dgamma,dweibull,erlang,expon,exponweib,exponpow,f,fatiguelife,fisk,foldcauchy,foldnorm,frechet_r,frechet_l,genlogistic,genpareto,genexpon,genextreme,gausshyper,gamma,gengamma,genhalflogistic,gilbrat,gompertz,gumbel_r,gumbel_l,halfcauchy,halflogistic,halfnorm,hypsecant,invgamma,invgauss,invweibull,johnsonsb,johnsonsu,ksone,kstwobign,laplace,logistic,loggamma,loglaplace,lognorm,lomax,maxwell,mielke,nakagami,ncx2,ncf,nct,norm,pareto,pearson3,powerlaw,powerlognorm,powernorm,rdist,reciprocal,rayleigh,rice,recipinvgauss,semicircular,t,triang,truncexpon,truncnorm,tukeylambda,uniform,vonmises,wald,weibull_min,weibull_max,wrapcauchy

if savedata:
    if HRF:
        datasets = [fa_h,fa_dr,fa_g]#[ca_h,ca_dr,ca_npdr,ca_pdr,ca_g]
        dataset_names = ["healthy","DR","glaucomatous"]#["healthy","diabetic retinopathy","NPDR","PDR","glaucomatous"]
    if DRIVE:
        datasets = [ca_h,ca_dr]
        dataset_names = ["healthy","diabetic retinopathy"]
    fitfuncs = [gamma,lognorm]
    fitfunc_names=["Gamma","Lognorm"]
    '''
    fitfuncs = [alpha,anglit,arcsine,beta,betaprime,bradford,burr,cauchy,chi,chi2,cosine,dgamma,dweibull,erlang,expon,exponweib,exponpow,f,fatiguelife,fisk,foldcauchy,foldnorm,frechet_r,frechet_l,genlogistic,genpareto,genexpon,genextreme,gausshyper,gamma,gengamma,genhalflogistic,gilbrat,gompertz,gumbel_r,gumbel_l,halfcauchy,halflogistic,halfnorm,hypsecant,invgamma,invgauss,invweibull,johnsonsb,johnsonsu,ksone,kstwobign,laplace,logistic,loggamma,loglaplace,lognorm,lomax,maxwell,mielke,nakagami,ncx2,ncf,nct,norm,pareto,pearson3,powerlaw,powerlognorm,powernorm,rdist,reciprocal,rayleigh,rice,recipinvgauss,semicircular,t,triang,truncexpon,truncnorm,tukeylambda,uniform,vonmises,wald,weibull_min,weibull_max,wrapcauchy]
    fitfunc_names = ['alpha',
 'anglit',
 'arcsine',
 'beta',
 'betaprime',
 'bradford',
 'burr',
 'cauchy',
 'chi',
 'chi2',
 'cosine',
 'dgamma',
 'dweibull',
 'erlang',
 'expon',
 'exponweib',
 'exponpow',
 'f',
 'fatiguelife',
 'fisk',
 'foldcauchy',
 'foldnorm',
 'frechet_r',
 'frechet_l',
 'genlogistic',
 'genpareto',
 'genexpon',
 'genextreme',
 'gausshyper',
 'gamma',
 'gengamma',
 'genhalflogistic',
 'gilbrat',
 'gompertz',
 'gumbel_r',
 'gumbel_l',
 'halfcauchy',
 'halflogistic',
 'halfnorm',
 'hypsecant',
 'invgamma',
 'invgauss',
 'invweibull',
 'johnsonsb',
 'johnsonsu',
 'ksone',
 'kstwobign',
 'laplace',
 'logistic',
 'loggamma',
 'loglaplace',
 'lognorm',
 'lomax',
 'maxwell',
 'mielke',
 'nakagami',
 'ncx2',
 'ncf',
 'nct',
 'norm',
 'pareto',
 'pearson3',
 
 'powerlaw',
 'powerlognorm',
 'powernorm',
 'rdist',
 'reciprocal',
 'rayleigh',
 'rice',
 'recipinvgauss',
 'semicircular',
 't',
 'triang',
 'truncexpon',
 'truncnorm',
 'tukeylambda',
 'uniform',
 'vonmises',
 'wald',
 'weibull_min',
 'weibull_max',
 'wrapcauchy']'''
    
    fwrite = open("{}_fit_results.txt".format(database),"w")
    for d in range(len(datasets)):
        fwrite.write("{},shape,loc,scale,k value,p value,ad k value,ad sig level,Ncells\n".format(dataset_names[d]))
        areas_norm = np.array(datasets[d])
        N = len(datasets[d])
        #mean = areas_array.mean()
        #areas_norm = areas_array/mean
        fig, ax = plot_linear_distribution(areas_norm,dataset_names[d],n="auto",density=True)
        ax.set_xlabel(r"$A/\langle A\rangle$")
        for i in range(len(fitfuncs)):
            #print(generalized_ks_test(areas_norm,fitfuncs[i],label=dataset_names[d]+" "+fitfunc_names[i],floc=floc))
            k, p, k_err, p_err = generalized_ks_test(areas_norm,fitfuncs[i],label=dataset_names[d]+" "+fitfunc_names[i],floc=floc,ax=ax,N=int(1e7),reps=2)
            ad_k, siglevel, stat_err, ad_k_err = generalized_ad_test(areas_norm,fitfuncs[i],label=dataset_names[d]+fitfunc_names[i],floc=floc,N=int(1e7),reps=2)
            if floc:
                params = fitfuncs[i].fit(areas_norm,floc=0)
            else:
                params = fitfuncs[i].fit(areas_norm)
            formatters = ["{:.3f}" for x in range(len(params))]
            fwrite.write("{},".format(fitfunc_names[i])+",".join(formatters).format(*params)+",{:.3f},{:.5f},{:.3f},{:5f},{:03d}\n".format(k,p,ad_k,siglevel, N))
        ax.legend()
    fwrite.close()
    
#%%
#fig = plt.figure()
#ax = fig.add_subplot(111)
#x = np.linspace(xmin,xmax,1000)
#colors = ["g","b","r"]
#dataset_names = ["healthy","diabetic retinopathy","glaucomatous"]
#
#i=0
#fmt=fmt_h
#fmt.update({"marker":"+"})
#fit = gamma.pdf(x,*f_params[i][:3])
#ax.plot(x,fit,colors[i]+"--",label=dataset_names[i]+" flat")
#plot_linear_distribution(ca_h,condition,n="auto",density=True,fig=fig,ax=ax,formatting=fmt)
#
#fmt.update({"marker":"o"})
#fit = gamma.pdf(x,*c_params[i][:3])
#ax.plot(x,fit,colors[i]+"-",label=dataset_names[i]+" curved")
#plot_linear_distribution(ca_h,condition,n="auto",density=True,fig=fig,ax=ax,formatting=fmt)
#
#ax.legend()