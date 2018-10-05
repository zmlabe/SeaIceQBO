# QBO and Arctic sea ice - perturbation study 

###### Under construction... ```[Python 3.6]```

## Contact
Zachary Labe - [Research Website](http://sites.uci.edu/zlabe/) - [@ZLabe](https://twitter.com/ZLabe)

## Description
Here we show that the phase of the Quasi-Biennial Oscillation (QBO) modulates the atmospheric response to Arctic sea ice loss. We conducted idealized experiments using WACCM4 to composite the phases of the QBO (westerly, easterly, and neutral) and assess the importance (mechanisms) of the stratosphere-troposphere pathway. The QBO in WACCM4 is prescribed by relaxing equatorial zonal winds between 86 and 4 hPa to observed radiosonde data (28-month period). We conduct a series of large ensemble simulations (7 experiments at 200 years each) to increase the signal-to-noise ratio in the stratosphere.

+ ```Scripts/```: Main [Python](https://www.python.org/) scripts/functions used in data analysis and plotting
+ ```requirements.txt```: List of environments and modules associated with the most recent version of this project. A Python [Anaconda3 Distribution](https://docs.continuum.io/anaconda/) was used for our analysis. All AGCM experiments were processed through resources on CISL's [Cheyenne](https://www2.cisl.ucar.edu/resources/computational-systems/cheyenne) supercomputer. Tools including [NCL](https://www.ncl.ucar.edu/), [CDO](https://code.mpimet.mpg.de/projects/cdo), and [NCO](http://nco.sourceforge.net/) were also used for initial data manipulation.

## Data
+ CESM Large Ensemble Project (LENS) : [[DATA]](http://www.cesm.ucar.edu/projects/community-projects/LENS/data-sets.html)
    + Kay, J. E and Coauthors, 2015: The Community Earth System Model (CESM) large ensemble project: A community resource for studying climate change in the presence of internal climate variability. Bull. Amer. Meteor. Soc., 96, 1333–1349, doi:10.1175/BAMS-D-13-00255.1 [[Publication]](http://journals.ametsoc.org/doi/full/10.1175/BAMS-D-13-00255.1)
+ Whole Atmosphere Community Climate Model (WACCM4) : [[DATA]](http://www.cesm.ucar.edu/working_groups/Whole-Atmosphere/code-release.html)
    + Marsh, D. R., M. J. Mills, D. E. Kinnison, J.-F. Lamarque, N. Calvo, and L. M. Polvani, 2013: Climate change from 1850 to 2005 simulated in CESM1 (WACCM). J. Climate, 26, 7372–7391, doi:10.1175/JCLI-D-12-00558.1 [[Publication]](http://journals.ametsoc.org/doi/abs/10.1175/BAMS-D-13-00255.1)
    + Smith, K. L., R. R. Neely, D. R. Marsh, and L. M. Polvani (2014): The Specified Chemistry Whole Atmosphere Community Climate Model (SC-WACCM), Journal of Advances in Modeling Earth Systems, 6(3), 883–901, doi:10.1002/2014MS000346 [[Publication]](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2014MS000346)


## Publications
+ **Labe, Z.M.**, Y. Peings, and G. Magnusdottir, 2018. Importance of the QBO on the large-scale atmospheric response to loss of Arctic sea ice, *in prep*


## Conferences
+ [2] **Labe, Z.M.**, G. Magnusdottir, and Y. Peings. Linking the Quasi-Biennial Oscillation and Projected Arctic Sea-Ice Loss to Stratospheric Variability in Early Winter, *20th Conference on Middle Atmosphere*, Phoenix, AZ (Jan 2019). [[Abstract]](https://ams.confex.com/ams/2019Annual/meetingapp.cgi/Paper/352664)
+ [1] Magnusdottir, G., **Z.M. Labe**, and Y. Peings. The role of the stratosphere, including the QBO, in Arctic to mid-latitude teleconnections associated with sea-ice forcing, *2018 American Geophysical Union Annual Meeting*, Washington, DC (Dec 2018). [[Abstract]](https://agu.confex.com/agu/fm18/meetingapp.cgi/Paper/399117)
