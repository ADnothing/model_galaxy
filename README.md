# model_galaxy
Modelisation on galaxy evolution using PEGASE-HR

This is a labwork done in internship at IAP directed by Damien Le Borgne (https://github.com/damleborgne)
Thoses codes use his work (https://github.com/damleborgne/PEGASE-HR and https://github.com/damleborgne/pypeg)

Formation_Galaxies : is a work (mostly in french) about the formation of galaxies

## Folder pypeg_modeling

lib : contains all the python's library used in the code

models : defins the diffrent model used for this project for scipy

modeling : contains the function 'plots' that plot the datas with the chosed models with optimized
           parameters and their variances.
           contains also the function 'plotQEXPY' which uses the qexpy library to fit the datas.
           
best_scen : code written by Damien Le Borgne with a modification to save in .csv the data collected (don't work alone, need other files and code from Damien Le Borgne)

main : if you have the data.csv file in the right path, ask you the model you want to test
       and uses the function 'plots'
       
alternatuve : is an alternative main, it uses the 'plotQEXPY' function of the file modeling

automat : is a code to use best_scen_MCMC automaticly to compute with a huge amount of steps and as much data as we want

best_scen_MCMC : is a modification of best_scen but it uses the MCMC method
       
The data.csv file here is just an example of the format used same thing for data_err.csv but with the uncertainty,
the column infall contains the infall time of a galaxy in Gyr,
the column sf contains the star formation time of a galaxy in Gyr,
the column winds contains the age of the galaxy where the galactitc winds begin and stop the star formation.

The last column is not used at the moment but using the PEGASE-HR simulation we still got it. We may use it later so we keep saving it.

The Param.csv file contains initials parameters and aimed parameters for an MCMC algorithm.

## Folder cosmo

In here, we use the COSMO catalogues (https://cosmos2020.calet.org/catalogues/). We use the file : COSMOS2020_CLASSIC_R1_v2.0.fits.

To understand the catalogues and how to use it, you should go to https://github.com/cosmic-dawn/cosmos2020-readcat and also read https://cosmos2020.calet.org/catalogues/COSMOS2020_CLASSIC_R1_v2.0.header.

libraies : contains all the libraries we use here.

select : selects the data of the subaru telescop for filters B, V and IB827 (http://svo2.cab.inta-csic.es/svo/theory/fps3/index.php?id=Subaru/Suprime.IB827&&mode=browse&gname=Subaru&gname2=Suprime). The datas are separated for diffrent redshift intervals (0<z<0.5 ; 0.5<z<1 ; ...) you don't have to run it if you want just dowload each zn.npz.

models : contains diffrent models we crated using pypeg.

ploting : plots the observations with the models

binning : file that Damien Le Borgne gaved us to plot the oberved galaxies as bins.
