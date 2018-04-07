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

```R
# Pneumonic cases from 22-09--03-10-2017
y <- c(1, 3, 6, 3, 5, 9, 7, 20, 25, 39, 53, 61)
# Use the same discrete generation time gamma(5.4, 0.9)
genTime <- generation.time("gamma", c(5.4, 0.9))
# Run 
R0_Viboud_Chowell(y, genTime, init=c(r=1, p=.5), lo=1e-8, up=c(1e+2, 1))
```

## Python notebook stochastic simulation of the plague transmission model

[Stochastic simulation](https://github.com/systemsmedicine/plague2017/blob/master/plague.ipynb)

[^1]: Viboud, C., Simonsen, L., & Chowell, G. (2016). A generalized-growth model to characterize the early ascending phase of infectious disease outbreaks. Epidemics, 15, 27–37. http://doi.org/10.1016/j.epidem.2016.01.002
[^2]: Chowell, G., Viboud, C., Simonsen, L., & Moghadas, S. M. (2016). Characterizing the reproduction number of epidemics with early subexponential growth dynamics. Journal of the Royal Society Interface, 13(123), 20160659. http://doi.org/10.1098/rsif.2016.0659
[^3]: Nishiura, H. (2010). Correcting the Actual Reproduction Number: A Simple Method to Estimate R0 from Early Epidemic Growth Data. International Journal of Environmental Research and Public Health, 7(1), 291–302. http://doi.org/10.3390/ijerph7010291
