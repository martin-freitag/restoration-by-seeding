<!-- README.md is generated from README.Rmd. Please edit that file -->

# Restoration of plant diversity in permanent grassland by seeding: assessing the limiting factors along land-use gradients [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4616821.svg)](https://doi.org/10.5281/zenodo.4616821)

Martin Freitag

18 March 2021

This project provides the R code resources to reproduce the analyses of
the manuscript

Freitag M, Klaus VH, Bolliger R, Hamer U, Kleinebecker T, Prati D,
Schäfer D, Hölzel N: *Restoration of plant diversity in permanent
grassland by seeding: assessing the limiting factors along land-use
gradients.* accepted in *Journal of Applied Ecology.*

# Structure

## Data

As per the Rules of Procedure of the Biodiversity Exploratories project,
the data used in this article is stored publicly in the *Biodiversity
Exploratories Information System* (BExIS) database along with the
metadata describing the data and its structure. You can download the
data as a .zip folder from
[BExIS](https://www.bexis.uni-jena.de/sws/PublicDataLink/Index), extract
the folder and put the .txt data files in the `data/` folder of this
project.

### Experimental design and land-use data

The **metadata.txt** file provides metadata to the vegetation records
sampled in the experiment over the five years. Columns correspond to
unique observations (*plotIdent*) sampled in the 73 grasslands (*plot*)
in three regions (*region*) over five years (*year*, *SADEyear*) and
two-factorial diverse-seeding and topsoil disturbance treatments
(*treatment*, dummy-coded *seed* and *dist*). Grassland-level land-use
intensity data includes fertilization (*Fert*, in kg Nitrogen
ha<sup>-1</sup>), grazing (*graz*, in livestock unit grazing days ×
ha<sup>-1</sup>) and mowing (*Mow*, average mowing frequency per year)
intensities Observation-level percentage covers of vegetation
(*cover.herbs.pc*), litter (*cover.litter.pc*) and bare soil
(*cover.bare.soil.pc*) and harvested biomass (*biomass.g*, in g biomass
m<sup>-2</sup>) is provided as well.

### Species cover data

Species cover data from vegetation records is available in long format
in the file **species-long.txt**. Each row corresponds to a *species*
and its *cover* of an observation (*plotIdent*).

### Trait data

We used plant functional trait data saved in **LEDA-traits.txt** from
the [LEDA traitbase](https://doi.org/10.1111/j.1365-2745.2008.01430.x)
(Kleyer et al. 2008). We averaged canopy height (*height*, in m), seed
mass (*seedMass*, in mg) and specific leaf area (*SLA*, in
mm<sup>2</sup> mg<sup>-1</sup>) measurements per *species*.

### Seed mixtures

The *sown-species.txt* file lists the species in the seed mixtures of
the three regions.

### Seeding density and seed germination rates

The **seeding-density-germination-rates.txt** file contains the seeding
densities of the species in the seed mixtures. Based on seed mass, we
have sown seeds in three densities (*seeddensity*, categorical, and
*seeddensity.per.m2*, numeric). We estimated seed viability by counting
the number of emerging seedlings in a germination test
(*germinationrate*, percentage of germinated seeds).

## Scripts

Note that I used Stan (<http://mc-stan.org/>) via the packages
[`rstan`](http://mc-stan.org/rstan/) and
[`brms`](https://github.com/paul-buerkner/brms) to analyse the data. To
install and use Stan, please follow the [installation
instructions](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

**01-1\_diversity-treatment-effects.R** We analysed how seeding and
topsoil disturbance influenced plant diversity over time. We used two
complementary diversity indices, namely species richness and effective
number of species S<sub>PIE</sub>.

**01-2\_land-use-effects-establishment-multivariate.R** We fitted two
separate multivariate models for seeding-only and seeding and topsoil
disturbance (combined) treatments to evaluate how grazing, mowing and
productivity determine the number of established seeding species. We use
`rstan` to fit the models, and the associated Stan script
*01-2\_delta-multivariate.stan* can be found in the *Stan\_files/*
folder.

**01-3\_establishment-traits-binomial.R** We modelled establishment of
sown species to explore how functional traits affect the establishment
of sown species five years after seeding along a productivity gradient.
We use `rstan` to fit the models, and the associated Stan script
*01-2\_traits-binomial.stan* can be found in the *Stan\_files/* folder.

**02-1\_create-figures-treatment.R**,
**02-2\_create-figures-multivariate.R** and
**02-3\_create-figures-traits.R** contain scripts to create figures for
the three above-mentioned analyses.

**03\_supplement\_delta-richness-2019-resident-richness.R** We
additionally analysed how the number of established species correlates
with the resident plant diversity (presented in the appendix of the
manuscript).

**03\_supplement\_litter-and-predictor-plots.R** Here, we create figures
of predictor distributions used in the analyses and other figures
presented in the appendix of th manuscript.

## Manuscript

The *doc/* folder contains all materials to reproduce the main text and
supplementary material. The *Freitag-et-al\_revision\_main-text.Rmd*
file contains the manuscript main text and imports figures from the
`figs/` folder. *Freitag-et-al\_revision\_supporting-information.Rmd*
contains additional figures and tables from the `tables/` folder. The
submitted manuscript PDF proofs are stored in the `doc/` folder as well.

## Package versioning

I use the `renv` package to keep track of package versions. The current
state of the project is saved in the lockfile `renv.lock`. After
installing `renv` you can restore the state of the project from the
`renv.lock` file with `renv::restore()` or update the lockfile after
updating packages with `renv::snapshot()`.
