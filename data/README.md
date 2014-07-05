# Some notes on dataset construction

## Mammals

The data set for mammals is built from data provided in (Ernest 2003). All mass data were converted to kg dry mass by 1 minus 0.605, the mean body water content for 11 mammals measured in vitro (Wang et al. 1999). Following (Purvis and Harvey 1995) and (Ernest et al. 2003), the following variables are calculated:

- Wa = adult body size (kg dry mass)
- W0 = size at weaning (kg dry mass)
- A = annual reproductive investment by an individual over one year (kg/individual/year) was calculated as: RI = fL NL W0reproductive lifespan where fL is the litter frequency per year, NL is the number of offspring per litter, and W0 is the mass of an offspring at independence (weaning). In addition, I multiplied by average pre-weaning survival of 0.5 (two refs suggest values of 0.45 (O'Donoghue 1994) and 0.63 (Gaillard et al. 2000) are appropriate ).
- M	= per month mortality, from Purvis 1995. I used this value of M to estimate average in preference to the values given in Ernest, as these are maximum lifespans.
- Average lifespans (E; yrs) was calculated using adult average mortality data supplied in  (Purvis and Harvey 1995)  [E = 1/(1-EXP(-M))/12] (Note: Usually would just convert as Ea = 1/M/12 but here they use this other formula. Maybe because survival was calculated on a yearly basis and want to get instantaneous rate).

Although consistent with relationships, data for the Blue whale (Balaenoptera musculus) was excluded because it lay sufficiently far from remainder of data that it would have an undue influence on regression parameters.

## Plants

I primarily used data compiled by Falster and Moles (Moles et al 2004), supplemented with data from Vile et al. (2006) and Niklas & Enquist (2004). If seed mass were missing, I sourced these from a database maintained by Angela Moles.

Compared to Moles et al. (2004), this dataset:

-	did not take species averages, data are species*site combinations
-	did not distinguish between mean and maximum parameters, where both available used means, otherwise merged
-	seed number converted via seed output (mass per year)/seed mass, or vice versa if both not given in original
-	above-ground mass adjusted to account for missing below-ground mass

For data from Niklas & Enquist (2004), we

-	Only included data which might be regarded as species’ size at adulthood
-	Mixed stands (<90% dominated by one spp), seedlings, crops excluded
-	For trees with multiple values I took the value closest to 20 years (size at maturity), rather than mass values for old trees with large accumulations of dead wood

Adjustment for roots

-	Original citations were checked to see whether mass measurement included below-ground structures with following codes
-	Analysis by (Niklas 2004) found average root:shoot ratio of 0.305 across large datasets
- Dry mass data were adjusted by multiplying by the following constants:
	* 1.0   (roots already included)
	* 1.305 (adjust for missing roots)
	* 1.15 (midpoint between above numbers, assuming 50:50 distribution)

### Cleaning

I removed some data points on the following grounds

- reproductive output was greater than plant mass (6 points)
- the data were outliers in relationship, and X and Y variables were from different studies, meaning I had little certainty in the match (5 points)

## References

- Ernest, S. K. M. 2003. Life history characteristics of placental nonvolant mammals. Ecology 84:3402-3402.
- Ernest, S. K. M., B. J. Enquist, J. H. Brown, E. L. Charnov, J. F. Gillooly, V. M. Savage, E. P. White, F. A. Smith, E. A. Hadly, J. P. Haskell, S. K. Lyons, B. A. Maurer, K. J. Niklas, and B. Tiffney. 2003. Thermodynamic and metabolic effects on the scaling of production and population energy use. Ecology Letters 6:990-995.
- Gaillard, J. M., M. Festa-Bianchet, N. G. Yoccoz, A. Loison, and C. Toigo. 2000. Temporal variation in fitness components and population dynamics of large herbivores. Annual Review of Ecology and Systematics 31:367-393.
- Henery, M., and M. Westoby. 2001. Seed mass and seed nutrient content as predictors of seed output variation between species. Oikos 92:479-490.
- Lord, J. M., and M. Westoby. 2006. Accessory costs of seed production. Oecologia 150:310-317.
- Moles, A., D. S. Falster, M. Leishman, and M. Westoby. 2004. Small-seeded plants produce more seeds per square metre of canopy per year, but not per individual per lifetime. Journal of Ecology 92:384–396.
- Niklas, K. J. 2004. Plant allometry: is there a grand unifying theory? Biological Reviews 79:871-889.
- Niklas, K. J., and B. J. Enquist. 2004. Biomass Allocation and Growth Data of Seeded Plants. Data set. Available on-line [http://www.daac.ornl.gov] from Oak Ridge National Laboratory Distributed Active Archive Center, Oak Ridge, Tennessee, U.S.A.
- O'Donoghue, M. 1994. Early Survival of Juvenile Snowshoe Hares. Ecology 75:1582-1592.
- Purvis, A., and P. H. Harvey. 1995. Mammal Life-History Evolution - a Comparative Test of Charnovs Model. Journal of Zoology 237:259-283.
- Vile, D., B. Shipley, and E. Garnier. 2006. A structural equation model to integrate changes in functional strategies during old-field succession. Ecology 87:504-517.
- Wang, Z., P. Deurenberg, W. Wang, A. Pietrobelli, R. N. Baumgartner, and S. B. Heymsfield. 1999. Hydration of fat-free body mass: review and critique of a classic body-composition constant. Am J Clin Nutr 69:833-841.
