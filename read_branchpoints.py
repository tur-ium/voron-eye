# -*- coding: utf-8 -*-
"""
Created on Thu Jan 17 10:12:39 2019

@author: admin
"""
import csv
def read_branchpoints(f):
    data = open(f,"r")
    
    bp_reader = csv.reader(data,delimiter=',')
    
    sites = []
    
    for row in bp_reader:
        coords = [int(float(x)) for x in row]
        sites.append(coords)
    return sites

def read_metadata(f):
    data = open(f,"r")
    
    bp_reader = csv.reader(data,delimiter=',')
    
    