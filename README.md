phcfM R package
===============

`phcfM` is an R package for modelling anthropogenic deforestation. It
was initially developed to obtain REDD+ baseline scenarios of
deforestation for the *programme holistique de conservation des forêts
à Madagascar* (from which the package is named after). It includes
two main functions:

1. `demography()`, to model the population growth with time in a
hierarchical Bayesian framework using population census data and
Gaussian linear mixed models.

2. `deforestation()`, to model the deforestation process in a
hierarchical Bayesian framework using land-cover change data and
Binomial logistic regression models with variable time-intervals
between land-cover observations.

Code and manual
---------------

The last stable version of the `phcfM` R package (`phcfM_1.2`) is
officially available for several operating systems (Unix, Windows and
Mac OSX) on the Comprehensive R Archive Network ([CRAN](http://cran.r-project.org/web/packages/phcfM/index.html)).

The source code and manual for the testing version (`phcfM_1.2-1`) is available here:

* Package source: [phcfM_1.2-1.tar.gz](http://sourceforge.net/projects/phcfm/files/phcfM_1.2-1.tar.gz/download)
* Reference manual: [phcfM_1.2-1.pdf](http://sourceforge.net/projects/phcfm/files/phcfM_1.2-1.pdf/download)

Example
-------

A GRASS location (`phcfM_SM`) and two mapsets with geographical data
layers (`PERMANENT` and `study_area_4`) are available to illustrate the
use of the `phcfM` R package. Associated to the GRASS location, a
directory (`./scripts`) includes the data and the R/GRASS scripts used
for the demographic and deforestation models.

* GRASS location: [phcfM_SM.zip](http://sourceforge.net/projects/phcfm/files/phcfM_SM.zip/download)
* ./scripts: [scripts.zip](http://sourceforge.net/projects/phcfm/files/scripts.zip/download)

The following scheme illustrates the structure of the `./scripts` folder.

./scripts  
|-- 0_demography_SAs.R  
|-- 1_deforestation_dataset_SA4.R  
|-- 2_deforestation_model_SA4.R  
|-- 3_deforestation_forecast_PS4.R  
|-- data  
|-- |-- 0_data_demography_SAs.txt  
|-- |-- 1_data_deforestation_SA4.txt  
|-- fragindex  
|-- |-- recl.txt  
|-- |-- r.forestfrag  
|-- outputs  
|-- |-- demography  
|-- |-- forecast  
|-- |-- model

Four scripts are available: `0_demography_SAs.R` is a R script to
estimate the parameters of the demographic model using the function
`demography()` of the `phcfM` R package,
`1_deforestation_dataset_SA4.R` is a R/GRASS script to prepare the
data-set for the deforestation model from the geographic data layers
of the `phcfM_SM` GRASS location, `2_deforestation_model_SA4.R` is a R
script to estimate the parameters of the deforestation model using the
function `deforestation()` of the `phcfM` R package and
`3_deforestation_forecast_PS4.R` is a R/GRASS script to forecast
deforestation and carbon dioxide emissions.

The `./scripts/data` folder includes the population census data file
`0_data_demography_SAs.txt` used for the demographic model. Running
the `1_deforestation_dataset_SA4.R` script generates the file
`1_data_deforestation_SA4.txt` (the data-set used for the
deforestation model) in the `./scripts/data` folder.  The
`./scripts/fragindex` folder includes the GRASS function
`r.forestfrag` to compute the forest fragmentation index. The
`./scripts/outputs` folder includes the outputs for the demographic
model, the deforestation model and the deforestation forecast.


Related publications
--------------------

**Vieilledent G., Grinand C., Vaudry R.** 2013. Forecasting
  deforestation and carbon emissions in tropical developing countries
  facing demographic expansion: a case study in Madagascar. *Ecology
  and Evolution*. DOI: 10.1002/ece3.550.

Contact
-------

Ghislain Vieilledent	
Cirad, UPR BSEF		
Campus de Baillarguet	
TA C-105/D		
34398 Montpellier cedex 5	
FRANCE

Skype: ghislain.vieilledent   
Email: ghislain(dot)vieilledent(at)cirad(dot)fr   
WWW: <http://ghislain.vieilledent.free.fr>


