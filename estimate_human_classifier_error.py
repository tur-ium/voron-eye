# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 19:15:13 2019

@author: Artur Donaldson

Estimates human error in classification of branchpoints
"""
from setup import bp_directory,areas_directory,im_directory
from read_branchpoints import read_branchpoints
import numpy as np
names = ["01_dr","01_dr","01_g","02_g"]
ad=list()
ss = list()
for name in names:
    directory = bp_directory+"HRF\\"
    
    branchpoints_path = directory+name+"_final_BPs.txt"
    sites_ss = read_branchpoints(branchpoints_path)
    ss.extend(sites_ss)
    
    branchpoints_path = directory+name+"_final_BPs_ad.txt"
    sites_ad = read_branchpoints(branchpoints_path)
    ad.extend(sites_ad)
all_bps = ad.copy()
all_bps.extend(ss)
print(len(all_bps))
x=0
all_bps.sort()
while x < len(all_bps)-1:
    if all_bps[x]==all_bps[x+1]:
        all_bps.pop(x+1)
    else:
        x+=1

def set_recur(l):
    s = set()
    for x in l:
        s.add(frozenset(x))
    return s
    
ad_set=set_recur(ad)
ss_set=set_recur(ss)
bp_set = set_recur(all_bps)
print("All points labelled BPs by either classifier",len(all_bps))
print("Classifier 1: ",len(ad))
print("Classifier 2: ",len(ss))
print("Classifiers differed in the labelling of {} points".format(len(ad_set.difference(ss_set))))
#ss = np.array(ss)
#ad = np.array(ad)
#
#N = len(ad)
#
#diff = (ss!=ad).sum()