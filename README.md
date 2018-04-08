# plague2017

Data and R code for Epidemics submission 2018

This repository is provided as a companion to the article on Epidemics  journal. The following code can be found:

## Two R code implementations estimating reproduction number

- Without exponential assumption (assumed binomial processes)[^3] using command `R0_Nishiura`

```R
# Pneumonic cases from 22-09--03-10-2017
y <- c(1, 3, 6, 3, 5, 9, 7, 20, 25, 39, 53, 61)
# Use the same discrete generation time gamma(5.4, 0.9)
genTime <- generation.time("gamma", c(5.4, 0.9))
# Run 
R0_Nishiura(y, genTime, init=1.5, lo=1.001)
# Profiling...
#       R0    2.5 %   97.5 % 
# 6.984100 5.206784 9.791975 
```

- Sub-exponential assumption[^1],[^2], using command `R0_Viboud_Chowell`

> Note that using daily data with this method lead to huge values of R0 during the first few days

```R
# Pneumonic cases from 22-09--03-10-2017
y <- c(1, 3, 6, 3, 5, 9, 7, 20, 25, 39, 53, 61)
# Use the same discrete generation time gamma(5.4, 0.9)
genTime <- generation.time("gamma", c(5.4, 0.9))
# Run 
out <- R0_Viboud_Chowell(y, genTime, init=c(r=1, p=.5), lo=1e-8, up=c(1e+2, 1))
#> out
#$r
#        r 
#0.4228841 #

#$p
#        p 
#0.9085438 #

#$Rt
#           2            3            4            5            6            7 
#8.947769e+09 8.131145e+04 3.565584e+02 2.820609e+01 9.957149e+00 7.068513e+00 
#           8            9           10           11           12 
#6.342856e+00 5.970402e+00 5.664230e+00 5.391341e+00 5.145617e+00 

```

## Python notebook stochastic simulation of the plague transmission model

[Stochastic simulation](https://github.com/systemsmedicine/plague2017/blob/master/plague.ipynb)

[^1]: Viboud, C., Simonsen, L., & Chowell, G. (2016). A generalized-growth model to characterize the early ascending phase of infectious disease outbreaks. Epidemics, 15, 27–37. http://doi.org/10.1016/j.epidem.2016.01.002

[^2]: Chowell, G., Viboud, C., Simonsen, L., & Moghadas, S. M. (2016). Characterizing the reproduction number of epidemics with early subexponential growth dynamics. Journal of the Royal Society Interface, 13(123), 20160659. http://doi.org/10.1098/rsif.2016.0659

[^3]: Nishiura, H. (2010). Correcting the Actual Reproduction Number: A Simple Method to Estimate R0 from Early Epidemic Growth Data. International Journal of Environmental Research and Public Health, 7(1), 291–302. http://doi.org/10.3390/ijerph7010291
