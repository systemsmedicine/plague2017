# R0 sub-exponentioal
# -------------------------------------------------------------------------
# 
# Viboud, C., Simonsen, L., & Chowell, G. (2016). A generalized-growth model to characterize the early ascending phase of infectious disease outbreaks. Epidemics, 15, 27–37. http://doi.org/10.1016/j.epidem.2016.01.002
# Chowell, G., Viboud, C., Simonsen, L., & Moghadas, S. M. (2016). Characterizing the reproduction number of epidemics with early subexponential growth dynamics. Journal of the Royal Society Interface, 13(123), 20160659. http://doi.org/10.1098/rsif.2016.0659
# 
# Note that daily data lead to huge values of R0 during the first days
# 
# -------------------------------------------------------------------------
R0_Viboud_Chowell <- function(.dat=y, .gen.time=genTime, init=c(r=1, p=.5), lo=1e-8, up=c(1e+2, 1),...) {
	# 1. Estimate r and p
	# -----------------------------------------------------------------------
	# C(t) = (r/m*t + A)^m, m = 1/(1-p), A = C0^{1/m}, C0: 1st observation
	modNLS <- nls(.dat ~ (r/(1/(1-p))*seq_along(.dat) + .dat[1]^(1/(1/(1-p))) )^(1/(1-p)), start=init, algorithm='port', lower=lo, upper=up,...)
	# 2. Simulate epidemic curves
	# -----------------------------------------------------------------------
	yhat <- predict(modNLS)
	# 3. Apply equation R_{t_i} = I_i / sum_{j=0}^i I_{i-j} rho_j
	# -----------------------------------------------------------------------
	names(yhat) <- 1:length(yhat)
	f_mat <- function(x) yhat[which(names(yhat)==x)]
	f_sum <- function(x) {
		ylist <- sapply(x-genTime$time, f_mat)
		return(sum(unlist(Map(`*`, ylist,  genTime$GT))))
	}
	f_Rt <- function(x) yhat[x]/f_sum(x)
	return(list(r=coef(modNLS)['r'], p=coef(modNLS)['p'], Rt=sapply(seq_along(yhat), f_Rt)[-1])) # exclude day 1 (Inf R0)
}
# Example
# -------------------------------------------------------------------------
# Pneumonic cases from 22-09--03-10-2017
y <- c(1, 3, 6, 3, 5, 9, 7, 20, 25, 39, 53, 61)
# Use the same discrete generation time gamma(5.4, 0.9)
genTime <- generation.time("gamma", c(5.4, 0.9))
# Run 
out <- R0_Viboud_Chowell(y, genTime, init=c(r=1, p=.5), lo=1e-8, up=c(1e+2, 1))
# > out
# $r
#         r 
# 0.4228841 

# $p
#         p 
# 0.9085438 

# $Rt
#            2            3            4            5            6            7 
# 8.947769e+09 8.131145e+04 3.565584e+02 2.820609e+01 9.957149e+00 7.068513e+00 
#            8            9           10           11           12 
# 6.342856e+00 5.970402e+00 5.664230e+00 5.391341e+00 5.145617e+00 
