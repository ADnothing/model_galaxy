# model_galaxy
Modelisation on galaxy evolution using PEGASE-HR

This is a labwork done in internship at IAP directed by Damien Le Borgne (https://github.com/damleborgne)
Thoses codes use his work (https://github.com/damleborgne/PEGASE-HR and https://github.com/damleborgne/pypeg)

lib : contains all the python's library used in the code

models : defins the diffrent model used for this project

modeling : contain the function 'plots' that plot the datas with the chosed models with optimized
           parameters and their variances
           
best_scen : code written by Damien Le Borgne with a modification to save in .csv the data collected (don't work alone, need other files and code from Damien Le Borgne)

main : if you have the data.csv file in the right path, ask you the model you want to test
       and uses the function 'plots'
       
The data.csv file here is just an example of the format used,
the column infall contains the infall time of a galaxy in Gyr,
the column sf contains the star formation time of a galaxy in Gyr,
the column winds contains the age of the galaxy where the galactitc winds begin and stop the star formation.

The last column is not used at the moment but using the PEGASE-HR simulation we still got it. We may use it later so we keep saving it.
