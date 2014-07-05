# A general model for the scaling of offspring size and adult size

[Daniel Falster](http://danielfalster.com/)
[Angela Moles](http://www.bees.unsw.edu.au/angela-moles)
[Mark Westoby](http://bio.mq.edu.au/research/groups/ecology//westoby/mark.htm),

This repository contains all the code used in the manuscript:

Falster, D.S., Moles, A.T. & Westoby, M. (2008) A general model for the scaling of offspring size and adult size. *American Naturalist* **172**: 299–317. DOI: [10.1086/589889](http://doi.org/10.1086/589889)

**Abstract:** Understanding evolutionary co-ordination among different life-history traits is a key challenge for ecology and evolution. Here, we develop a general quantitative model predicting how offspring size should scale with adult size by combining a simple model for life-history evolution with a frequency-dependent survivorship model. The key innovation is that larger offspring are afforded three different advantages during ontogeny: higher survivorship per time, a shortened juvenile phase, and advantage during size-competitive growth. In this model, it turns out that size-asymmetric advantage during competition is the factor driving evolution towards larger offspring sizes. For simplified and limiting cases the model is shown to produce the same predictions as the previously existing theory it is founded on. The explicit treatment of different survival advantages has biologically important new effects, mainly through an interaction between total maternal investment in reproduction and the duration of competitive growth. This leads on to explain alternative allometries between log offspring size and log adult size, as observed in mammals (slope=0.95) and plants (slope=0.54). Further, it suggests how these differences relate quantitatively to specific biological processes during recruitment. In these ways the model generalizes across previous theory and provides explanations for some differences between major taxa.

## Running

To reproduce the figures from the paper, open a terminal window and run the command

```
make
```

Alternative you can open `R` and type

```
source(`analysis.R`)
```

A separate file  `analytical_solutions.m` contains matlab code used to derive the analytical solution reported in the paper.

## Citation

For archival purposes, the code used to produce figures for publication has been lodged with figshare [here](). This is the same as code included in [this github release]().

To cite this code:

```

```

## Further details

- This code was archived post-publication, so there is no link to the code at the journal's website.
- The code included in this release has been slightly restructured to improve readability and code quality. Original code was written between 2006-2008 by Daniel Falster, then updated in 2014. The version run at time of publication is available [here](https://github.com/dfalster/Falster_2008_AmNat_offspring_model/tree/c06c5de3a54b4589581b2f74b2c9fc5d1529fd6d) (NB - does not include data). While I have endeavoured to make the code more readable, it's still a little arcane - my apologies for that.
- More information about the data used in this analysis is available in the file [data/README.md](data/README.md)
