# -*- coding: utf-8 -*-
"""
Created on Sun Mar 17 19:48:38 2019

@author: admin
"""
import numpy as np
from scipy.stats import ks_2samp,anderson_ksamp

def generalized_ks_test(data,model,reps=10,N=int(1e6),label="%UNKNOWN%",ax=None,floc=False):
    """Test how well a given scipy.stats distribution function fits a given set\
data using the Kolmogorov-Smirnoff (KS) test
    PARAMS
    ---
    data: list
        A sample of data
    model: scipy.stats continuous distribtion object
    reps: number of repetitions of samples to take when testing the fit (default: 10)
    N: Number of data points to use in the model (default: 1 000 000)
    label: Name when printing results to console
    ax: Axes
    floc: if True forces the location parameter to be 0 otherwise unconstrained
   RETURNS
   ---
   mean k-statistic, KS p-value, std of mean k-statistic, std of KS p-value
"""
    xmin = min(data)
    xmax = max(data)
    if floc:
        params = model.fit(data,floc=0)
    else:
        params = model.fit(data)
    formatters = ["{:.3f}" for x in range(len(params))]
    print("PARAMS OF FIT for {}: ".format(label)+",".join(formatters).format(*params))
    
    x = np.linspace(xmin,xmax,1000)
    fit = model.pdf(x,*params)
    if ax:
        ax.plot(x,fit,label=label+", floc {}".format("On" if floc else "Off"))
    kstats,ps = list(), list()
    for i in range(reps):
        sample = model.rvs(*params,size=N)
        kstat,p = ks_2samp(sample,data)
        #print("RESULTS OF KS TEST,it.",i,": (kstat,p) = ", kstat,p)
        kstats.append(kstat)
        ps.append(p)
    print("2 sample Kolmogorov-Smirnoff test: (k,pval,k_std,p_std): {:.3f}, {:.5f}, {:.5f}, {:.5f}".format(np.mean(kstats),np.mean(ps),np.std(kstats),np.std(ps)))
    return (np.mean(kstats),np.mean(ps),np.std(kstats),np.std(ps))

def generalized_ad_test(data,model,reps=10,N=int(1e6),label="%UNKNOWN%",ax=None,floc=False):
    """Test how well a given scipy.stats distribution function fits a given set\
data using the Anderson-Darling (AD)
    PARAMS
    ---
    data: list
        A sample of data
    model: scipy.stats continuous distribtion object
    reps: number of repetitions of samples to take when testing the fit (default: 10)
    N: Number of data points to use in the model (default: 1 000 000)
    label: Name when printing results to console
    ax: Axes
    floc: if True forces the location parameter to be 0 otherwise unconstrained
   RETURNS
   ---
   mean statistic, AD siglevel, std of mean statistic, std of AD p-value
"""
    xmin = min(data)
    xmax = max(data)
    if floc:
        params = model.fit(data,floc=0)
    else:
        params = model.fit(data)
    formatters = ["{:.3f}" for x in range(len(params))]
    print("PARAMS OF FIT for {}: ".format(label)+",".join(formatters).format(*params))
    x = np.linspace(xmin,xmax,100)
    fit = model.pdf(x,*params)
    if ax:
        ax.plot(x,fit,label=label)
    stats,siglevels = list(),list()
    for i in range(reps):
        sample = model.rvs(*params,size=N)
        stat,threshs,siglevel = anderson_ksamp((sample,data))
        stats.append(stat)
        siglevels.append(siglevel)
    print("Results of Anderson-Darling test: (stat,siglevel,stat_std,siglevel_std): {:.3f}, {:.5f}, {:.5f}, {:.5f}".format(np.mean(stats),np.mean(siglevels),np.std(stats),np.std(siglevels)))
    return np.mean(stats),np.mean(siglevels),np.std(stats),np.std(siglevels)