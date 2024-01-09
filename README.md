# Dune_AM
MATLAB scripts and data associated with the paper titled: "Mesoscale Structure of the Atmospheric Boundary Layer Across a Natural Roughness Transition"
What you will find here: The data associated with the alkali flat and dune field, as well as the MATLAB script used to create Figures 2B, 3A, 4B, and 4C. All the data that was used to make the observations and conduct the analysis in the paper is included. 

#------- The Data -------#

DuneField directory contains velocity data and wall shear stress data taken within the dune field portion of the domain. Data in DuneData are for the streamwise component (comp(u,0)) of the velocity and the x,y,z coordinates (README files). Files are named by x_XXXX.* giving the streamwise location in meters with reference to the start of the domain (x=0).
In the WSS directory there is a compressed file called WallStress.tar.xz which contains wall shear stress data taken along the center line of the domain. The corresponding README file contains the probe x,y,z coordinates for each data point. The SWSS_xXXXX.* files contain wall shear stress data with x constant and varying in y along the spanwise direction. This is done to get a sense of the average wall stress across multiple dunes. However, a more comprehensive version of this data is in the Slice directory, which contains many more x-locations, but the same type of information. This is used to create the blue triangles in Fig 3A.

In the AlkaliFlat directory we have streamwise component of the velocity data. However, this was done in a weird manner with probes constant along z instead of x. So each file is labeled z_XXX.X.* which gives the data at a constant elevation. The data was originally grabbed at MANY locations in x, which is contained in OrigData, please make sure to unzip all the data from z_XXX.X.zip. The data was then consolidated to 4 data points, x = 1000m, 1250, 1500m, 1750m. Both sets of data are used in the analysis, as the OrigData is loaded, but only the data corresponding to the 4 x-locations are grabbed and used. Unfortunately this makes the code a little slow to start since the file sizes are large. Sorry, but I didn't want to mess around with the data or code after the plots were completed, so sit back, enjoy some tea or coffee, and just wait the few minutes for the data to load.

File Formatting:

README FILES: 
First Col: Probe Number 
Second Col: X-Coordinate
Third Col: Y-Coordinate 
Fourth Col: Z-Coordinate

comp(u,0) FILES:
First Col: Step number
Second Col: Simulation Time
Third Col: Total number of data points
Fourth Col-End Col: Data corresponding to either x varying or z varying data


#------- The Script -------#

The script is a MATLAB script compatible with MATLAB 2019b and subsequent releases. I tried to comment it as much as I could to ensure that anyone could read the code and understand what was going on. If the code doesn't run on the first pass, please read the error codes from MATLAB and try and figure out what might be the issue. If things persist, please then feel free to contact me at justincooke96@gmail.com. At the end of the code, the Figures used in the manuscript are output. Other figures included in the code are for sanity checks and early mock-ups of figures that could be interesting to you. 

Thanks for your interest in my work! If you found any of this data or the analysis insightful, please make sure to cite the paper! 
