phcfM R package
===============

`phcfM` is an R package for modelling anthropogenic deforestation. It
was initially developed to obtain REDD baseline scenarios of
deforestation for the *programme holistique de conservation des forêts
à Madagascar* (from which the package was named after). It includes
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

* Package source: [phcfM_1.0.tar.gz](http://sourceforge.net/projects/phcfm/files/phcfM_1.0.tar.gz/download)
* Windows binary (R 2.14.2): [phcfM_1.0.zip](http://sourceforge.net/projects/phcfm/files/phcfM_1.0.zip/download)
* Reference manual: [phcfM.pdf](http://sourceforge.net/projects/phcfm/files/phcfM.pdf/download)

Example
-------

A GRASS location (`phcfM_SM`) and two mapsets with geographical data
layers (`PERMANENT` and `study_area_4`) are available to illustrate the
use of the `phcfM` R package. Associated to the GRASS location, a
directory (`./scripts`) includes the data and the R/GRASS scripts used
for the demographic and deforestation models.

* GRASS location: [phcfM_SM.zip](http://sourceforge.net/projects/phcfm/files/phcfM_SM.zip/download)
* ./scripts: [scripts.zip](http://sourceforge.net/projects/phcfm/files/scripts.zip/download)

Related publications
--------------------

**Vieilledent G., Grinand C., Vaudry R.** Forecasting anthropogenic deforestation and carbon emissions in tropical forest. In prep.

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


