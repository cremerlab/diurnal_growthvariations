We did most of the data processing, analysis, and plotting with Python scriipts we developed for this purpose. 

Data to start with
--------------------
The scripts are based on types of data exported from COMSOL.

1. The full spatio-temporal data on bacteria/nutrient concentration, fermentation products (if appicable), and the velocities (radial and along z). The variables are named c2,c3,c4, u, and w. We exported these files from COMSOL using the export function. To use the data we then transformed them into a simpler to read DataFrame (Script 1 below) with the data exported from COMSOL saved in the folder data_fromcomsol/full_variables_spatiotemporal_resolution

DESCRIBE EXPORT SETTINGS HERE.  

2. The variation of different integrated model variables over time. 

COMSOL allows for the convinient integration of variables over different simulation domains. Examples include the total amount of bacteria or nutrients in the system, or the flux of nutrients and bacteria entering or exiting via the influx valve and the distal end of the simulation domain. These files are read in and used in different plotting striptcs with the data exported from COMSOL saved in the folder data_fromcomsol/systems_variables


Scripts
--------
To run these scripts with explemlarity simulation data you can download the simulation data from the data repository of this work, having the different subfolders (data_fromcomsol, data_kymographss, etc in the same folder than the script files.)


Script 1: read_data.py
-----------------------
-----------------------
Transform the export from COMSOL to a easy to read dataframe. The resulting dataframes are stored in folder data_simpleformat.

Script 2: generate_averages_and_kymograph_data.py
-------------------------------------------------
--------------------------------------------------
Generate kymograph data (calculated radially averaged values of variables or take values at centerline). Output is stored in data_kymographs.

Script 3: plot_data_kymographs_.ipynb
-------------------------------------
-------------------------------------
Plot kymographs, snaptshows, and trends over time.

Script 4: plot_data_movies.py
--------------------------------
--------------------------------
Plot sequence of images for videos. All images for one video will be stored in a subfolder under vidos. The script gerate_movies_from_images.ipynb will be used to use ffmpg and then generate a video file.

Script 5: plot_data_trends_SCFA_and_cecum.ipynb
-----------------------------------------------
-----------------------------------------------
Plot different trends over time

Script 6: plot_data_kymographs-differentgrowthrate.ipynb
--------------------------------------------------------
--------------------------------------------------------
Plot change with growth rate.



