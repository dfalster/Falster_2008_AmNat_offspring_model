# A general model for the scaling of offspring size and adult size

[Daniel Falster](http://danielfalster.com/)
[Angela Moles](http://www.bees.unsw.edu.au/angela-moles)
[Mark Westoby](http://bio.mq.edu.au/research/groups/ecology//westoby/mark.htm)

This repository contains all the code used in the manuscript:
  
Falster, D.S., Moles, A.T. & Westoby, M. (2008) A general model for the scaling of offspring size and adult size. *American Naturalist* **172**: 299–317. DOI: [10.1086/589889](http://doi.org/10.1086/589889)

**Abstract:** Understanding evolutionary co-ordination among different life-history traits is a key challenge for ecology and evolution. Here, we develop a general quantitative model predicting how offspring size should scale with adult size by combining a simple model for life-history evolution with a frequency-dependent survivorship model. The key innovation is that larger offspring are afforded three different advantages during ontogeny: higher survivorship per time, a shortened juvenile phase, and advantage during size-competitive growth. In this model, it turns out that size-asymmetric advantage during competition is the factor driving evolution towards larger offspring sizes. For simplified and limiting cases the model is shown to produce the same predictions as the previously existing theory it is founded on. The explicit treatment of different survival advantages has biologically important new effects, mainly through an interaction between total maternal investment in reproduction and the duration of competitive growth. This leads on to explain alternative allometries between log offspring size and log adult size, as observed in mammals (slope=0.95) and plants (slope=0.54). Further, it suggests how these differences relate quantitatively to specific biological processes during recruitment. In these ways the model generalizes across previous theory and provides explanations for some differences between major taxa.

## Running the code

All analyses were done in `R`. You can reproduce the results for this project by running the code in the `analysis.R` file. Figures will be output to a directory called `output`. 

If reproducing these results on your own machine, you much first install the package [smatr](cran.r-project.org/package=smatr), described in [Warton et al 2013](http://doi.org/10.1111/j.2041-210X.2011.00153.x):
  
```
install.packages('smatr')
```

You can access an interactive RStudio session with the required software pre-installed by opening a container hosted by [Binder](http://mybinder.org): 

[![Launch Rstudio Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dfalster/Falster_2008_AmNat_offspring_model/master?urlpath=rstudio)

To ensure long-term [computational reproducibility](https://www.britishecologicalsociety.org/wp-content/uploads/2017/12/guide-to-reproducible-code.pdf) of this work, we have created a [Docker](http://dockerhub.com) image to enable others to reproduce these results on their local machines using the same software and versions we used to conduct the original analysis. Instructions for reproducing this work using the docker image are available at bottom of the page. 

## Material included in the repository include:

- `data/`: Raw data
- `R/` directory containing functions used in analysis
- `DECRIPTION`: A machine-readable [compendium]() file containing key metadata and dependencies 
- `license.md`: License for the materials
- `Dockerfile` & `.binder/Dockerfile`: files used to generate docker containers for long-term reproducibility
- `analytical_solutions.m`: contains matlab code used to derive the analytical solution reported in the paper

## Further details

- The equations in the Figures 4-6 were added outside of R so will not be visible in the generated figures.
- This code was archived post-publication, so there is no link to the code at the journal's website.
- The code included in this release has been slightly restructured to improve readability and code quality. Original code was written between 2006-2008 by Daniel Falster, then updated in 2014. The version run at time of publication is available [here](https://github.com/dfalster/Falster_2008_AmNat_offspring_model/tree/c06c5de3a54b4589581b2f74b2c9fc5d1529fd6d) (NB - does not include data). While I have endeavoured to make the code more readable, it's still a little arcane - my apologies for that.
- More information about the data used in this analysis is available in the file [data/README.md](data/README.md)

## Citation

For archival purposes, the code used to produce figures for publication has been lodged with figshare [here](http://dx.doi.org/10.6084/m9.figshare.1094315). This is the same as code included in [this github release](https://github.com/dfalster/Falster_2008_AmNat_offspring_model/releases/tag/1.0).

To cite this code:
  
```
Falster, Daniel (2014): Code from: A general model for the scaling of offspring size and adult size. figshare. http://dx.doi.org/10.6084/m9.figshare.1094315
```

## Running via Docker

If you have Docker installed, you can recreate the compute environment as follows. 

First clone this repository:

```
git clone https://github.com/traitecoevo/Falster_2008_AmNat_offspring_model.git
```

Then fetch the container:

```
docker pull traitecoevo/falster_2008_amnat_offspring_model
```

From within your downloaded repo, launch the container using the following code (it will map your current working directory inside the docker container): 

```
docker run --user root -v $(pwd):/home/rstudio/ -p 8787:8787 -e DISABLE_AUTH=true traitecoevo/falster_2008_amnat_offspring_model
```

The code above initialises a docker container, which runs an rstudio session, which is accessed by pointing your browser to [localhost:8787](http://localhost:8787). For more instructions on running docker, see the info from [rocker](https://hub.docker.com/r/rocker/rstudio).

### NOTE: Building the docker image

For posterity, the docker image was built off [`rocker/verse:3.6.1` container](https://hub.docker.com/r/rocker/verse) via the following command, in a terminal contained within the downloaded repo:

```
docker build -t traitecoevo/falster_2008_amnat_offspring_model .
```

and was then pushed to [dockerhub](https://cloud.docker.com/u/traitecoevo/repository/docker/traitecoevo/falster_2008_amnat_offspring_model). The image used by binder builds off this container, adding extra features needed bi binder, as described in [rocker/binder](https://hub.docker.com/r/rocker/binder/dockerfile).

