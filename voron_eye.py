# -*- coding: utf-8 -*-
"""
Created on Thu Dec 20 12:31:22 2018

@author: admin
"""
#OWN CODE
from curved import pixel_area_on_retina,lat,jacobian,magnification_factor

#Required packages
from scipy.spatial import Voronoi, voronoi_plot_2d
from matplotlib import pyplot as plt
import numpy as np
from shapely.geometry import Polygon, Point
from descartes import PolygonPatch
from curved import cartesian_to_polar, distance_on_retina
from numpy import pi, sin, cos
from mpl_toolkits.mplot3d import Axes3D

epsilon = 1e-10 #Small value for floating point comparison
 

# =============================================================================
# TESTING: RANDOM POINTS ON SQUARE IMAGE
# =============================================================================
##Seed for demonstration
#np.random.seed(seed=14)
#number_of_sites = 100
#
#im_w = 2000
#FOV = 45*np.pi/180
#
#R_im = int(.5*im_w)
#R_eye = R_im/np.cos(.5*np.pi-.5*FOV)
#
#sites = list()
#for n in range(number_of_sites):
#    site = np.random.randint(-R_im,R_im,size=2)
#    sites.append(site)
#sites = np.array(sites)
#
##Will close all matplotlib windows
#plt.close("all")
#
##Create mask
#square_mask = Polygon([[0,0],[im_w,0],[im_w,im_w],[0,im_w]])
#th = np.linspace(0,2*np.pi,100)
#x,y = np.cos(th), np.sin(th)
#x_T  = R_im*x.reshape((1,len(x)))
#y_T = R_im*y.reshape((1,len(y)))
#point_array = np.concatenate((x_T,y_T),0).T
#circular_mask = Polygon(point_array)
#

def compare_flat_and_curved_distributions(voronoi_sites,mask,R_eye,R_im,im_w,im_h,FOV,x_0,y_0,R_optic,plot=True,threshold = epsilon):
    """ Calculates the distribution of the areas of Voronoi cells which are not\
    on the boundary of the image. Requires a mask as a Shapely Polygon object.
    
    PARAMETERS
    ---
    voronoi_sites: nested list
        List of 2d coordinates of Voronoi sites of the form [[x0,y0],[x1,y1]]
    mask: Shapely Polygon
        A polygon whose enclosed region represents the field of view of the camera
    R_eye: float
        Radius of eye
    R_im: float
        Radius of image (in px)
    FOV: float
        Field of View (in radians)
    x_0,y_0,R_optic: int, int, float
        x,y position of centre of optic disc (from bottom left corner) and the radius of the optic disc
    plot: boolean
        Plots the Voronoi diagram if True
    threshold: float
        Threshold for comparison of the boundary sites. If the ratio of the difference between
the areas of the cells to the cell area is less than this, the full region will be included in the analysis
    RETURNS
    ---
    flat_core_cell_area_dist, curved_core_cell_area_dist: np.array, np.array
        Core Voronoi cell areas for flat and curved geometries as a numpy array
        
    """
    flat_area_dict = dict()
    curved_area_dict = dict() #Maps Voronoi region index to area
    #Using a dictionary here to maintain one-to-one correspondence between V.regions \
    # and the areas. Boundary cells will not be included in here
    patches = list()

    if plot:
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)

    site_array = np.array(voronoi_sites)       
    # =============================================================================
    # REMOVE OPTIC DISC
    # =============================================================================
    
    th = np.linspace(0,2*np.pi,10)
    x,y =  R_optic*np.cos(th)+x_0,  R_optic*np.sin(th)+y_0
    x_T  = x.reshape((1,len(x)))
    y_T = y.reshape((1,len(y)))
    point_array = np.concatenate((x_T,y_T),0).T
        
    #Add boundary of the optic disc to the list of vascular branchpoints so as \
#to exclude this area from Voronoi cells
    voronoi_sites = np.concatenate((site_array,point_array))
    
    optic_disc = Polygon(point_array)
    mask = mask.difference(optic_disc)
    
    if plot:
        optic_disc_poly = PolygonPatch(optic_disc,facecolor=(1,0,0),alpha=.9)
        ax1.add_patch(optic_disc_poly)
        
    # =============================================================================
    # GENERATE VORONOI DIAGRAM ON 2D IMAGE
    # =============================================================================
    V = Voronoi(voronoi_sites,incremental=True)
    N_r = len(V.regions)
    
    def indicesV_to_vertices(voronoi,indices):
        return [voronoi.vertices[i].tolist() for i in indices]
    
    def region_to_site(voronoi):
        """Returns a dictionary mapping Voronoi regions on the image to the Voronoi site"""
        return {voronoi.point_region[i]: i for i in range(len(voronoi.points))}
    
    #Plot Voronoi diagram
    if plot:
        voronoi_plot_2d(V,ax1,show_vertices=False)
    
    # =============================================================================
    # GET CORE CELLS
    # =============================================================================
    regions_to_sites = region_to_site(V)
    for r in range(N_r):
        vertices = indicesV_to_vertices(V,V.regions[r])
    
        #Check #1 for whether the region is on the boundary of the image
        if len(vertices)>2 and -1 not in V.regions[r]:
            voronoi_polygon = Polygon(vertices)
            #Check #2 for if the region is on the boundary of the image
            poly = voronoi_polygon.intersection(mask)
            voronoi_site = regions_to_sites[r]
            x,y = V.points[voronoi_site]
#            if plot:
#                ax1.text(x+.3,y+.3,voronoi_site)
            
            #print(abs(poly.area - voronoi_polygon.area)/voronoi_polygon.area)
            if abs(poly.area - voronoi_polygon.area)/voronoi_polygon.area < threshold:
                #Find area on image
                area_on_image = poly.area
                
                flat_area_dict[r] = area_on_image
                #print(x,y)
                S = pixel_area_on_retina(x,y,R_eye,R_im,FOV,lat,jacobian)
                
                area_on_eye = area_on_image*S
                curved_area_dict[r]=area_on_eye
                if plot: 
                    patches.append(poly)                    
                    polygon = PolygonPatch(poly,facecolor=(0,np.random.random()*.8+.2,np.random.random()*.8+.2),alpha=.5)
                    ax1.add_patch(polygon)
            elif poly.area > 0:
                if plot:
                    polygon = PolygonPatch(poly,facecolor=(1,0,0),alpha=.5)
                    ax1.add_patch(polygon)
    
    
    flat_core_cell_areas = np.array(list(flat_area_dict.values()))
    curved_core_cell_areas = np.array(list(curved_area_dict.values()))
    
    total_core_cell_area = sum(curved_core_cell_areas)
    #Analytical solution for total area on retina in Field of View
    magnification = magnification_factor(R_im,R_eye,FOV)
    A_a = 2*np.pi*(magnification*R_eye)**2*(1-np.cos(.5*FOV))
    
    if plot:
        mask_polygon = PolygonPatch(mask,fc="#ffffff00",ec="#ff0000ff")
        ax1.add_patch(mask_polygon)
        ax1.set_aspect("equal")
        ax1.set_xlim(-.5*im_w,.5*im_w)
        ax1.set_ylim(-.5*im_h,.5*im_h)
        fig1.suptitle("Voronoi cells")
        ax1.set_xlabel(r"$x$ (px)")
        ax1.set_ylabel(r"$y$ (px)")
        fig1.subplots_adjust(bottom=.25)
        fig1.text(0,0.01,"Core cell area: {:.0f} px ({:.2%}).".format(total_core_cell_area,total_core_cell_area/A_a))
        fig1.text(0,0.06,"Green/blue cells are core cells")
        fig1.text(0,0.11,"Red cells are boundary cells, white cells are outside of the field of view")
    return flat_core_cell_areas, curved_core_cell_areas
def curved_core_cell_areas(voronoi_sites,mask,R_eye,FOV,x_0,y_0,R_optic,phi_N=100,lat_N=100,plot=True,threshold=None,verbose=False):
    """
    voronoi_sites: nested list
        List of 2d coordinates of Voronoi sites of the form [[x0,y0],[x1,y1]]
    mask: Shapely Polygon
        A polygon whose enclosed region represents the field of view of the camera
    R_eye: float
        Radius of eye
    R_im: float
        Radius of image (in px) (NOT IMPLEMENTED)
    FOV: float
        Field of View (in radians)
    phi_N,lat_N: int
        Number of phi and latitude gridlines to use
    x_0,y_0,R_optic: int, int, float
        x,y position of centre of optic disc (from bottom left corner) and the radius of the optic disc
    plot: boolean
        Plots the Voronoi diagram if True
    threshold: float
        Threshold for comparison of the boundary sites. If the ratio of the difference between
the areas of the cells to the cell area is less than this, the full region will be included in the analysis
    RETURNS
    ---
    curved_core_cell_area_dist: np array
        Core Voronoi cell areas for curved geometry as a numpy array"""
    # =============================================================================
    # REMOVE OPTIC DISC
    # =============================================================================
    site_array = np.array(voronoi_sites)       
   
    th = np.linspace(0,2*np.pi,10)
    x,y =  R_optic*np.cos(th)+x_0,  R_optic*np.sin(th)+y_0
    x_T  = x.reshape((1,len(x)))
    y_T = y.reshape((1,len(y)))
    optic_disc_boundary = np.concatenate((x_T,y_T),0).T
        
    #Add boundary of the optic disc to the list of vascular branchpoints so as \
#to exclude this area from Voronoi cells
    voronoi_sites = np.concatenate((site_array,optic_disc_boundary))
    
    optic_disc = Polygon(optic_disc_boundary)
    mask = mask.difference(optic_disc)
    

    # =============================================================================
    # CALCULATED CONSTANTS
    # =============================================================================
    lat_min = .5*pi-.5*FOV
    lat_max = .5*pi
    d_phi = 2*pi/phi_N
    d_lat = (lat_max-lat_min)/lat_N
    
    R_im = R_eye * cos(lat_min)
    # =============================================================================
    # PLOT PARAMS
    # =============================================================================
    print("Plotting: {}".format(plot))
    if plot:
        plt.rcParams.update({"font.size":12})
    
    s_sites = list()
    # =============================================================================
    # CONVERT CARTESIAN COORDS TO POLAR COORDS
    # =============================================================================
    for c_coord in voronoi_sites:
    
        s_coord = list(cartesian_to_polar(*c_coord,R_eye))
        s_sites.append(s_coord)
    # =============================================================================
    # CALCULATE Z COORDS OF POINTS
    # =============================================================================
    z_sites = R_eye*sin(np.array(s_sites)[:,1])
    
    # =============================================================================
    # Define grid
    # =============================================================================
    phis,lats = np.ogrid[-pi:pi:phi_N*1j,lat_min:lat_max:lat_N*1j]
    regions = np.zeros((phi_N,lat_N))
    
    # =============================================================================
    # CONVERT SPHERICAL POLAR GRID TO CARTESIAN (FOR COMPARISON WITH MASK AND PLOTTING)
    # =============================================================================
    def polar_to_cartesian(phi,lat,r):
        x=r*cos(lat)*cos(phi)
        y=r*cos(lat)*sin(phi)
        z=r*sin(lat)
        return (x,y,z)
    R_c = np.ndarray((phi_N,lat_N,3)) #Cartesian coordinates for each coordinate in the grid
    xs,ys,zs = list(),list(),list()
    for i1 in range(phi_N):
        for i2 in range(lat_N):
            x,y,z = polar_to_cartesian(phis[i1,0],lats[0,i2],R_eye)
            R_c[i1,i2,:] = x,y,z
            xs.append(x)
            ys.append(y)
            zs.append(z)
    
    def check_dist(coord1,coords2,dist_func=distance_on_retina):
        """coord1: (2,1) np.array
            Coordinate to be checked
        coord2: (2xN) np.array
            Coordindates to be checked"""
        min_d = np.inf
        closest_site = None
        N = coords2.shape[0]
        for i in range(N):
            d = abs(dist_func(coord1[0],coord1[1],coords2[i,0],coords2[i,1],R_eye))
            if d < min_d:
                min_d = d
                closest_site = i
        return closest_site
    
    for i1 in range(phi_N):
        for i2 in range(lat_N):
            min_d = np.inf
            closest_site = None
            #print(R_c[i1,i2,0:2])
            if mask.contains(Point(R_c[i1,i2,0:2])):
                for i3 in range(len(s_sites)):
                    d = abs(distance_on_retina(phis[i1,0],lats[0,i2],*s_sites[i3],R_eye))
                    if d < min_d:
                        min_d = d
                        closest_site = i3
                if closest_site is not None:
                    regions[i1,i2]=closest_site
    # =============================================================================
    # PLOT REGIONS IN 3D
    # =============================================================================
        
    if plot:
        fig = plt.figure()
        ax = plt.subplot(111,projection="3d")
        ax.scatter(xs,ys,zs,c=regions.ravel(),s=1)
        
        # =============================================================================
        # PLOT SITES
        # =============================================================================
        for i in range(voronoi_sites.shape[0]):
            ax.scatter(voronoi_sites[i,0],voronoi_sites[i,1],z_sites[i],marker="x",label="Site "+str(i))
        
        # =============================================================================
        # PLOT SETTINGS
        # =============================================================================
        ax.set_aspect("equal")
        ax.set_xlim([-R_eye,R_eye])
        ax.set_ylim([-R_eye,R_eye])
        ax.set_zlim([-R_eye,R_eye])
        ax.set_xlabel(r"$x$")
        ax.set_ylabel(r"$y$")
        ax.set_zlabel(r"$z$")
    
    # =============================================================================
    # CALCULATE AREAS
    # =============================================================================
    def area_on_eye(lam,R,d_phi,d_angle):
        return R**2*cos(lam)*d_phi*d_lat
    
    region_areas = np.zeros(len(s_sites))
    print(region_areas.shape)
    for i1 in range(phi_N):
        for i2 in range(lat_N):
            if mask.contains(Point(R_c[i1,i2,0:2])):
                i3 = int(regions[i1,i2])
#                if i3 not in region_areas:
#                    region_areas[i3] = 0
                region_areas[i3] += area_on_eye(lats[0,i2],R_eye,d_phi,d_lat)
    if verbose:
            
        # =============================================================================
        # TEST 2: CHECK TOTAL AREA OF CONIC SECTION = ANALYTICAL VALUE
        # =============================================================================
        A_analytic = 2*pi*R_eye**2*(sin(lat_max)-sin(lat_min)) #Analytic solution
        A = sum(region_areas.values())
        print("Total area")
        print("Numerical solution: ", A)
        print("Analytic solution: ",A_analytic)        
        
        # =============================================================================
        # TEST 3: CHECK AGAINST SPHERICAL VORONOI PACKAGE
        # =============================================================================
        #if plot:
            #from scipy.spatial import SphericalVoronoi
            #site_coords = np.concatenate((np.array(voronoi_sites),z_sites.reshape(len(voronoi_sites),1)),axis=1)
            #V = SphericalVoronoi(site_coords,R_eye)
            #ax.scatter(V.vertices[:, 0], V.vertices[:, 1], V.vertices[:, 2],c='g')
    
    return regions,region_areas
def core_cell_distribution(voronoi_sites,mask,R_eye,R_im,im_w,im_h,FOV,x_0,y_0,R_optic,plot=True):
    """ Calculates the distribution of the areas of Voronoi cells which are not\
    on the boundary of the image. Requires a mask as a Shapely Polygon object.
    PARAMETERS
    ---
    voronoi_sites: nested list
        List of 2d coordinates of Voronoi sites of the form [[x0,y0],[x1,y1]]
    mask: Shapely Polygon
        A polygon whose enclosed region represents the field of view of the camera
    plot: boolean
        Plots the Voronoi diagram if True
    R_eye: float
        Radius of eye
    R_im: int
        Radius of retina on image (in px)
    im_w, im_h: int, int
        Width and height of image
    FOV: float
        Field of View (in radians)
    x_0, y_0, R_optic: int, int, float
        Position and radius of the optic disc
    plot: boolean
       If True plots Voronoi cell diagram 
    RETURNS
    ---
    Core Voronoi cell areas as a numpy array
    """
    area_dict = dict() #Maps Voronoi region index to area
    #Using a dictionary here to maintain one-to-one correspondence between V.regions \
    # and the areas. Boundary cells will not be included in here
    patches = list()

    if plot:
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)

    site_array = np.array(voronoi_sites)       
    # =============================================================================
    # REMOVE OPTIC DISC
    # =============================================================================
    #Shift origin
    x_0 = x_0 -.5*im_w
    y_0 = .5*im_h - y_0
    
    th = np.linspace(0,2*np.pi,20)
    x,y = np.cos(th), np.sin(th)
    x_T  = R_optic*x.reshape((1,len(x)))+x_0
    y_T = R_optic*y.reshape((1,len(y)))+y_0
    
    point_array = np.concatenate((x_T,y_T),0).T
    
    #Add boundary of the optic disc to the list of vascular branchpoints so as \
#to exclude this area from Voronoi cells
    voronoi_sites = np.concatenate((site_array,point_array))
    
    optic_disc = Polygon(point_array)
    mask = mask.difference(optic_disc)
    
    if plot:
        optic_disc_poly = PolygonPatch(optic_disc,facecolor=(1,0,0),alpha=.9)
        ax1.add_patch(optic_disc_poly)
        
    # =============================================================================
    # GENERATE VORONOI DIAGRAM ON 2D IMAGE
    # =============================================================================
    V = Voronoi(voronoi_sites,incremental=True)
    N_r = len(V.regions)
    
    def indicesV_to_vertices(voronoi,indices):
        return [voronoi.vertices[i].tolist() for i in indices]
    
    def region_to_site(voronoi):
        """Returns a dictionary mapping Voronoi regions on the image to the Voronoi site"""
        return {voronoi.point_region[i]: i for i in range(len(voronoi.points))}
    
    #Plot Voronoi diagram
    if plot:
        voronoi_plot_2d(V,ax1)
    
    # =============================================================================
    # GET CORE CELLS
    # =============================================================================
    for r in range(N_r):
        vertices = indicesV_to_vertices(V,V.regions[r])
        
        #Check #1 for whether the region is on the boundary of the image
        if len(vertices)>2 and -1 not in V.regions[r]:
            voronoi_polygon = Polygon(vertices)
            #Check #2 for if the region is on the boundary of the image
            poly = voronoi_polygon.intersection(mask)
            voronoi_site = region_to_site(V)[r]
            x,y = V.points[voronoi_site]
            #if plot:
                #ax1.text(x+.3,y+.3,voronoi_site)
            
            #print(abs(poly.area - voronoi_polygon.area)/voronoi_polygon.area)
            if abs(poly.area - voronoi_polygon.area)/voronoi_polygon.area < threshold:
                #Find area on image
                area_on_image = poly.area
                
                #print(x,y)
                S = pixel_area_on_retina(x,y,R_eye,R_im,FOV,lat,jacobian)
                
                area_on_eye = area_on_image*S
                area_dict[r]=area_on_eye
                if plot: 
                    patches.append(poly)                    
                    polygon = PolygonPatch(poly,facecolor=(0,np.random.random()*.8+.2,np.random.random()*.8+.2),alpha=.5)
                    ax1.add_patch(polygon)
            elif poly.area > 0:
                if plot:
                    polygon = PolygonPatch(poly,facecolor=(1,0,0),alpha=.5)
                    ax1.add_patch(polygon)
    
    core_cell_areas = np.array(list(area_dict.values()))
    
    total_core_cell_area = sum(core_cell_areas)
    #Analytical solution for total area on retina in Field of View
    magnification = magnification_factor(R_im,R_eye,FOV)
    A_a = 2*np.pi*(magnification*R_eye)**2*(1-np.cos(.5*FOV))
    
    if plot:
        mask_polygon = PolygonPatch(mask,fc="#ffffff00",ec="#ff0000ff")
        ax1.add_patch(mask_polygon)
        ax1.set_aspect("equal")
        ax1.set_xlim(-.5*im_w,.5*im_w)
        ax1.set_ylim(-.5*im_h,.5*im_h)
        fig1.suptitle("Voronoi cells")
        ax1.set_xlabel(r"$x$ (px)")
        ax1.set_ylabel(r"$y$ (px)")
        fig1.subplots_adjust(bottom=.25)
        fig1.text(0,0.01,"Core cell area: {:.0f} px ({:.2%}).".format(total_core_cell_area,total_core_cell_area/A_a))
        fig1.text(0,0.06,"Green/blue cells are core cells")
        fig1.text(0,0.11,"Red cells are boundary cells, white cells are outside of the field of view")
    return core_cell_areas
#
#core_cell_areas = core_cell_distribution(sites,circular_mask)
#av_cell_area = sum(core_cell_areas)/len(core_cell_areas)
#core_cell_areas_norm = core_cell_areas/av_cell_area
#print("Average core Voronoi cell area: {:.0f} px".format(av_cell_area))

# =============================================================================
# PLOT DISTRIBUTION
# =============================================================================

def plot_log_distribution(array,name,n=100,norm=False,fig=None,ax=None):
    '''Plots loglog frequency and cumulative distributions for a quantity
    PARAMETERS
    ---
    d: np.array
        Numpy array of values to be histogrammed
    name: string
        Name of the quantity represented by the values of the dictionary
    n: int
        Number of bins
    norm: boolean
        Normalize the distribution by the average value
    fig: plt.figure object
        Figure to plot on
    ax: plt.axis
        Axis to plot on
    ''' 

    if norm:
        average = sum(array)/len(array)
        array = array/average
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
    log_max_degree = np.log10(max(array))
    bin_edges = np.logspace(0,log_max_degree,n-1)
    bin_edges = list(bin_edges)
    #bin_edges.insert(0,0.)
    
    degree_dist,bins = np.histogram(array,bins=bin_edges,density=False)
    prob_density,bins = np.histogram(array,bins=bin_edges,density=True)
    #prob_dist = degree_dist/max(degree_dist)
    #print(prob_dist)
    bin_centres = [.5*sum(bin_edges[i:i+2]) for i in range(0,len(bin_edges)-1)]
    #plt.semilogx(bin_centres,prob_dist,'b+',label="Mass density")
    
    ax.semilogx(bin_centres,prob_density,'o',label="Prob. density {}".format(name.capitalize()))
    ax.set_ylim(0,max(prob_density))
    ax.grid()
    fig.suptitle('{} distribution'.format(name.capitalize()))
    ax.set_xlabel('log({})'.format(name))
    ax.set_ylabel(r'Probability density')
    ax.legend()
    
    return fig,ax

def plot_linear_distribution(array,name,n=100,norm=False,density=False,fig=None,ax=None):
    '''Plots linear frequency and cumulative distributions for a quantity
    PARAMETERS
    ---
    d: np.array
        Numpy array of values to be histogrammed
    name: string
        Name of the quantity represented by the values of the dictionary
    n: int
        Number of bins
    norm: boolean
        Normalize the distribution by the maximum
    fig: plt.figure object
        Figure to plot on
    ax: plt.axis
        Axis to plot on
    RETURNS
    ---
    fig, ax: figure and axis of the plot
    ''' 
    if not ax:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
    if norm:
        average = sum(array)/len(array)
        array = array/average
    
    max_degree = max(array)
    bin_edges = np.linspace(0,max_degree,n-1)
    bin_edges = list(bin_edges)
    #bin_edges.insert(0,0.)
    
    #PLOT PROBABILITY DENSITY
    if density:
        prob_density,bins = np.histogram(array,bins=bin_edges,density=True) #Probability density
        bin_centres = [.5*sum(bin_edges[i:i+2]) for i in range(0,len(bin_edges)-1)]
        ax.set_ylabel(r'Probability Density')
        ax.plot(bin_centres,prob_density,'o',label="Prob. density {}".format(name.capitalize()))
        ax.set_ylim(0,max(prob_density)*1.1)
        ax.grid()
    #PLOT NORMALIZED PROBABILITY
    if norm:
        degree_dist,bins = np.histogram(array,bins=bin_edges,density=False)
        area = sum([degree_dist[i]*(bin_edges[i+1]-bin_edges[i]) for i in range(len(degree_dist))])
        #prob_dist = degree_dist/max(degree_dist) #Normalized probability distribution
        prob_dist = degree_dist/area
        bin_centres = [.5*sum(bin_edges[i:i+2]) for i in range(0,len(bin_edges)-1)]
        #plt.semilogx(bin_centres,prob_dist,'b+',label="Mass density")
        ax.set_ylabel(r'Probability/Max. probability')
        ax.plot(bin_centres,prob_dist,'o',label="{}".format(name.capitalize()))
        ax.set_ylim(0,max(prob_dist)*1.1)
        ax.grid()
    fig.suptitle('{} distribution'.format(name.capitalize()))
    ax.set_xlabel('{}'.format(name))
    
    ax.legend()
    
    return fig,ax

def calculate_and_compare_distributions(array,ref_dist,name,siglevel=.05,n=100,norm=False,density=False,plot=False,fig=None,ax=None):
    '''Calculates histogram for an array of values `array` and compares it with the reference distribution `ref_dist` using the \
Kolmogorov-Smirnoff test
    PARAMETERS
    ---
    d: np.array
        Numpy array of values to be histogrammed
    ref_dist: function
        Reference distribution function
    name: string
        Name of the quantity represented by the values of the dictionary
    n: int
        Number of bins
    norm: boolean
        Normalize the distribution by the maximum
    ''' 
    if norm:
        average = sum(array)/len(array)
        array = array/average
    max_degree = max(array)
    bin_edges = np.linspace(0,max_degree,n)
    bin_edges = list(bin_edges)
    degree_dist,bins = np.histogram(array,bins=bin_edges,density=False)
    area = sum([degree_dist[i]*(bin_edges[i+1]-bin_edges[i]) for i in range(len(degree_dist))])
    prob_dist = degree_dist/area
    
    bin_centres = np.array([.5*sum(bins[i:i+2]) for i in range(0,len(bins)-1)])
#    values_for_gamma = np.array([sum(bins[i:i+1]) for i in range(0,len(bins)-1)])
    #scipy_ks = ks_2samp(prob_dist,ref_dist(values_for_gamma))
    
#    if scipy_ks.pvalue>siglevel:
#        print("p-value: {} > {}, so the two distributions match".format(scipy_ks.pvalue,siglevel))
#    else:
#        print("p-value: {} < {}, so the samples are drawn from different distributions".format(scipy_ks.pvalue,siglevel))
#    
    if plot:
        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        
        ax.scatter(bin_centres*2,prob_dist,label="{}".format(name))
        smooth_xs = 2*np.linspace(min(bin_centres),max(bin_centres),50)
        ax.plot(smooth_xs,ref_dist(smooth_xs),"g--",label="Random")
        ax.set_ylim((0,1))
        fig.suptitle("{}".format(name))
        ax.set_ylabel(r'Normalized probability, p s.t.$\int p = 1$')
        plt.xlabel(r"Normalized area $A/\langle A\rangle$")
        ax.legend()
        
def kolmogorov_smirnoff_test(d1,d2,siglevel,n=0,m=0):
    """Calculates Kolmogorov-Smirnoff statistic of n randomly-selected samples from \
distribution d1 and m samples from d2. n=0 (m=0) indicates to use all values from d1 (d2)
RETURNS
---
    True if the two distributions differ by more than the significance level, otherwise\
returns False
"""
    if n==0:
        n = len(d1)
    if m==0:
        m=len(d2)
    if n!=m:
        raise Exception("Not Implemented")
    stat = max(abs(d1-d2))
    def c(alpha):
        return (-.5*np.log(alpha))**.5
    
    print("Sig level: ",siglevel)
    print("K-S stat: ",stat)
    print(" Threshold: ", c(siglevel)*(n+m)**.5/(n*m)**.5)
    if stat > c(siglevel)*(n+m)**.5/(n*m)**.5:
        return False
    else:
        return True