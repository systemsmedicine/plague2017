# R0 without exponential assumption
# -------------------------------------------------------------------------
# Nishiura, H. (2010). Correcting the Actual Reproduction Number: A Simple Method to Estimate R0 from Early Epidemic Growth Data. International Journal of Environmental Research and Public Health, 7(1), 291–302. http://doi.org/10.3390/ijerph7010291
# -------------------------------------------------------------------------
R0_epidemic <- function(.dat=y, .gen.time=genTime, init=2, lo=1.01, ...) {
  require('stats4')
  names(.dat) <- 1:length(.dat)
	# -----------------------------------------------------------------------
	# Calculate the total number of effective contacts made by potential primary cases in time t that have an equal probability of resulting in a secondary transmission, 
	# 
	# Data are the number of secondary cases in time t.
	# -----------------------------------------------------------------------
	f_mat <- function(x) .dat[which(names(.dat)==x)]
	f_sum <- function(x) {
    ylist <- sapply(x-genTime$time, f_mat)
    return(sum(unlist(Map(`*`, ylist,  genTime$GT))))
	}
	# Negative log likelihood of Binomial distribution
	# -----------------------------------------------------------------------
	nBnll <- function(R0) {
		f_ll <- function(x) f_sum(x)*log(1/R0) + (.dat[x]-f_sum(x))*log(1-1/R0)
		return(-sum(sapply(seq_along(.dat), f_ll)))
	} 
	# -----------------------------------------------------------------------
	fit <- mle(nBnll, start=list(R0=init), method='L-BFGS-B', lower=lo, ...)
	return(c(coef(fit), confint(fit)))
}

# Example
# -------------------------------------------------------------------------
# Pneumonic cases from 22-09--03-10-2017
y <- c(1, 3, 6, 3, 5, 9, 7, 20, 25, 39, 53, 61)
# Use the same discrete generation time gamma(5.4, 0.9)
genTime <- generation.time("gamma", c(5.4, 0.9))
# Run 
R0_epidemic(y, genTime, init=1.5, lo=1.001)
# Profiling...
#       R0    2.5 %   97.5 % 
# 6.984100 5.206784 9.791975 
