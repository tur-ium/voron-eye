
# Voron-Eye 
## A Simple Package for Voronoi Analysis of Retinal Blood Vessels

Voron-Eye is a Python package which enables calculation of Voronoi cell area distributions in spherical geometry using masking. This software was written by Artur Donaldson as part a project on analysis of transport networks. The project used voronoi analysis of the branchpoints of retinal blood vessel system of the eye to distinguish between different conditions.

The data and code are made available here for future work. Recent studies have shown that replication of scientific work is a serious issue in science. The working hypothesis is that publishing the code and data alongside papers can help to address the issue by making it easier to peer review and give a foundation for future work.

I would like to thank my project partner, Shervin Sabeghi for taking on the data collection side of this project and my supervisors Dimitri Vvedensky (Professor in Condensed Matter Theory, Imperial College London) and Simon Morgan (visiting researcher, Imperial College London) for their time and guidance making the results accurate.
 
REQUIRED PACKAGES
-----------------
numpy
matplotlib
scipy.spatial
scipy.stats
shapely
descartes

FILE STRUCTURE
----------
*curved.py*
Contains code for transforming lengths and areas on the 2D fundus image to physical lengths and areas on the curved surface of the retina.

*generalized_analysis.py*
Includes functions for checking the goodness of fit of a scipy.stats distribution to a sample of a data using the 2-sided Kolmogorov-Smirnoff test or the Anderson-Darling test.

*neovascularization_simulation.py*
A simple simulation for neovascularization on the retina. Shows that the effect of neovascularization is to shift the peak normalized area of the distribution to lower values.

*voron_eye.py*
Computes probability distributions of Voronoi cell areas given a set of branchpoints and a mask. Also comes with a plotting function plot_linear_distribution

DIR *branchpoint_data/*
Contains data for the branchpoints of the vascular networks for two datasets - HRF and DRIVE.

DIR *voronoi_cell_areas/*
Contains voronoi cell areas calculated using link_to_ibranch.py

CONTACT
-------
For more information or to report bugs contact Artur Donaldson (@tur-ium)
