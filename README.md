# Repository for corroboree frog (*Pseudophryne pengilleyi*) analysis

This repo contains four scripts for analyses contained in the manuscript. All analyses were conducted with `nimble` 0.13 in `R` 4.3.

  1.  `corroboree-analysis-age.R` fits negative binomial distributions to the age data from two extant and one museum population to estimate the mean age structure.
  2.  `corroboree-analysis-prevalence.R` fits a logistic regression to swab samples analysed for the presence of *Bd*. There are two intercept corresponding to the two sites and four random year effects (shared across sites; note that site Warogong-Cotterill just has one year of data). We then compute marginal means for each site and plot the predicted prevalence for each year.
  3.  `corroboree-analysis-mortality.R` estimated baseline annual mortality rates using the museum age structure which is then used to specify the intercept of a geometric survival model to estimate the effect of *Bd* infection on the mortality rates of laboratory-exposed frogs.
  4.  `corroboree-analysis-elasticity.R` contains the elasticity analysis populated with the survival estimates from (3) above to investigate the effect of a delayed age of maturation on the elasticity of adult survival.
  
Data for these files are contained within the `data` folder. MCMC output for the models in (1--3) above are contained in the `mcmc` folder. Figure output is in the `figs` folder.