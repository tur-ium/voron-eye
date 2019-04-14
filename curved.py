# -*- coding: utf-8 -*-
"""
Created on Fri Nov  9 14:54:46 2018

@author: Artur Donaldson
"""
import numpy as np

def cartesian_to_polar(x,y,r):
    """
    Convert projection of coordinates (x,y) onto cartesian plane to spherical polar 
        (phi,lat)
    
    PARAMS
    ---
    x,y: Cartesian coords
    r: Radius of curvature
    
    RETURNS
    ---
    (phi,lat): angle from x-axis in the x-y plane, latitude
    """
    
    lat= np.arccos((x**2+y**2)**.5/r)
    if x>0 and y>0:
        phi = np.arctan(y/x)
    elif x<0 and y>0:
        phi = np.arctan(y/x)+np.pi
    elif x<0 and y<0:
        phi = -np.pi+np.arctan(y/x)
    elif x>0 and y<0:
        phi = np.arctan(y/x)
    else:
        phi = 0
    return (phi,lat)

def magnification_factor(R_im,R_eye,FOV):
    M = R_im/(R_eye*np.sin(.5*FOV))   #Magnification factor
    return M

def lat(r_im,R_eye,M):
    """Returns latitude of point on eye at a distance on the image R_im, given
        eye radius R_eye and magnification factor M using Donder's reduced eye model"""
    return np.arccos(r_im/(M*R_eye))

def jacobian(r_im,R_eye,M):
    """Jacobian for transformation from 2D Euclidean space to 2D spherical space"""
    return ((M*R_eye)**2-r_im**2)**-.5*r_im**-1

def pixel_area_on_retina(x,y,R_eye,R_image,FOV,f,J):
    """Finds area corresponding to a pixel at (x,y) on the image to the
    corresponding are on the surface of the retina. R = Radius of eye, R_image is radius of image.
    f is the mapping from the distance from the centre of the image to the latitude on the image

    PARAMETERS:
    ---
    x, y: float
        Coordinates of pixel
    R_eye: float
        Radius of eye
    R_image: float
        Radius of image
    f: function
        Mapping from the distance from the centre of the image to the latitude on the image
    J: function
        Magnitude of the determinant of the Jacobian for the mapping
    """
    r_im = (x**2+y**2)**.5
    M = magnification_factor(R_image,R_eye,FOV)

    lat = f(r_im,R_eye,M)
    J = J(r_im,R_eye,M)

    return R_eye**2*np.cos(lat)*J

def distance_on_retina(phi1,lambda1,phi2,lambda2,R_eye):
    """Returns distance between two points (phi1,lambda1),(phi2,lambda2) on the surface of a spherical eye\
    of radius R_eye"""
    alpha = np.arccos(np.sin(lambda1)*np.sin(lambda2)+np.cos(lambda1)*np.cos(lambda2)*np.cos(np.abs(phi2-phi1)))
    return R_eye*alpha
