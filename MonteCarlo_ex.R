# Michael Gunnarson
# Summer 2022

rm(list = ls())


# forty RTD's placed in 0.0 C ice bath, room temp = 30 C

R = c(110.2, 112.4, 113.6, 108.5, 116.0, 112.8, 111.1, 111.2, 
	111.2, 110.7, 111.3, 110.5, 110.8, 112.1, 113.4, 111.8, 
	111.2, 110.6, 110.7, 115.6, 112.2, 111.3, 111.3, 112.7, 
	111.0, 109.8, 112.9, 110.5, 112.0, 110.6, 113.7, 111.3, 
	110.9, 111.2, 109.6, 110.2, 108.7, 110.0, 111.6, 111.3)


dev.new()
plot(R, xlab = "RTD Number", ylab = "Resistance (Ohms)", main = "RTD")


### Functions



Newt <- function(ma) {
	# create function "func" st
	# if y = f(x)
	# func = f(x)-y, where y is fixed and x is the desired solution
	# found iteratively

	# initialize

	count = 0
	res = 1
	
	while(res > 1e-6 & count < 10000) {
		count = count +1
	
	
	
		if(count == 1) {
			ma_last = ma - 0.1 # manufacture one point
			f1 = func(ma_last)
		}

		f2 = func(ma)
		df = (f2-f1)/(ma-ma_last)
		ma_new = ma - f2/df
	
		res = abs(f2)

		# shift variables
		f1 = f2
		ma_last = ma
		ma = ma_new
		# at end of it all, desired mach number is ma_last
	}
	return(list(ma_last,res,count))
}


# test data:
# g = seq(from = 0, to = 4, length.out = 1000)
# f = x^2


trapz  <- function(x,y) {
	# trapezoidal integrator based on already known x,y data
	# x must be of equal distance, i don't think R could handle
	# trying to integrate anything else.  It couldn't handle the check
	# last time I tried it.
	dx = x[2]-x[1]
	end = length(x)

	A = dx/2*y[1] + sum(dx*y[2:(end-1)]) + dx/2*y[end]
	return(A)
	
}

# trapz(g,f)


mean_p_boot <- function(x,b,P) {
	# parametric bootstrap based off of code from professor shaw
	# at present, only takes on b value
	# beta values assumed normal distrbution
	# parametric is where the curve is assumed normal pdf distribution
	# nonparametric is where the resampling creates its own pdf distribution

	M = 1e4
	S = numeric(M)
	N = length(x)
	xbar = mean(x)
	sx = sd(x)
	for (m in 1:M) {
		xSample = rnorm(n = N, mean = xbar, sd = sx)
		beta = rnorm(n=1, mean=0, sd = b)
		xSample = xSample + beta
		S[m] <- mean(xSample)
	}
	Pmin = (1-P)/2
	Pmax = (1+P)/2
	
	return(list(mean(S),quantile(S, probs = c(Pmin,Pmax))))
	# quantile is 95% confidence level
	# mean is mean(S) with a 95% confidence level that it lies between
	# quantile values		
}


sd_p_boot <- function(x,b,P) {
	# parametric bootstrap based off of code from professor shaw
	# at present, only takes on b value
	# beta values assumed normal distrbution

	M = 1e4
	S = numeric(M)
	N = length(x)
	xbar = mean(x)
	sx = sd(x)
	for (m in 1:M) {
		xSample = rnorm(n = N, mean = xbar, sd = sx)
		beta = rnorm(n=1, mean=0, sd = b)
		xSample = xSample + beta
		S[m] <- sd(xSample)
	}
	Pmin = (1-P)/2
	Pmax = (1+P)/2
	
	return(list(mean(S),quantile(S, probs = c(Pmin,Pmax))))
	# quantile is 95% confidence level
	# the standard deviation is sd(S) with a 95% confidence level that it lies between
	# quantile values		
}



mean_np_boot <- function(x,b,P) {
	# nonparametric bootstrap based off of code from professor shaw
	# at present, only takes on b value
	# beta values assumed normal distrbution

	M = 1e4
	S = numeric(M)
	N = length(x)
	xbar = mean(x)
	sx = sd(x)
	for (m in 1:M) {
		xSample = sample(x, size = N, replace = TRUE)
		beta = rnorm(n=1, mean=0, sd = b)
		xSample = xSample + beta
		S[m] <- mean(xSample)
	}
	Pmin = (1-P)/2
	Pmax = (1+P)/2
	
	return(list(mean(S),quantile(S, probs = c(Pmin,Pmax))))
	# quantile is 95% confidence level
	# mean is mean(S) with a 95% confidence level that it lies between
	# quantile values		
}


sd_np_boot <- function(x,b,P) {
	# parametric bootstrap based off of code from professor shaw
	# at present, only takes on b value
	# beta values assumed normal distrbution

	M = 1e4
	S = numeric(M)
	N = length(x)
	xbar = mean(x)
	sx = sd(x)
	for (m in 1:M) {
		xSample = sample(x, size = N, replace = TRUE)
		beta = rnorm(n=1, mean=0, sd = b)
		xSample = xSample + beta
		S[m] <- sd(xSample)
	}
	Pmin = (1-P)/2
	Pmax = (1+P)/2
	
	return(list(mean(S),quantile(S, probs = c(Pmin,Pmax))))
	# quantile is 95% confidence level
	# the standard deviation is sd(S) with a 95% confidence level that it lies between
	# quantile values		
}





### 1.  Identify/quantify errors, random, systematic, in-evaluatable

random_errors = c('Resolution')
systematic_errors = c('Accuracy', '0.9% + 2')
inevaluatable_errors = c('Linearity')

random_errors
systematic_errors
inevaluatable_errors

### 2.  Identify outliers using Chauvenet's criterion, remove at end

m = mean(R)
N = length(R)
sdR = sd(R)

# translate to z:

z = (R-m)/sdR
P = 1-1/(2*N)

func <- function(zeta) {
	# set output of zeta function to be P/2
	# write function such that it is equal to zero (g = f(x) - y)
	int_fun <- function(zeta) {
		# integral part of the funciton, to be called by
		# trapezoidal rule function
		int = exp(-zeta^2/2)
		return(int)
	}
	x = seq(from = 0, to = zeta, length.out = 10000)
	y = int_fun(x)
	zer = 1/(2*pi)^0.5 * trapz(x,y) - P/2

}

z_max = Newt(0.2)[[1]] # newton's method to solve for Z_max

remove = R[abs(z) >= z_max]
chauv = R[abs(z) < z_max]

mchv = mean(chauv)
Nchv = length(chauv)
sdchv  = sd(chauv)

m = mchv
N = Nchv
sdR = sdchv

### 3.  Evaluate 95% uncertainty interval for the mean of population using
      # non parametric M.C. approach.  Parametric MC can be used for systematic errors
	

# with code from professor Shaw:

# given that:
# Accuracy = beta = 0.9% + 2 of reading


find_beta <- function(R) {
	# 0.9% + 2

	end = length(R)
	betaMin = numeric(end)
	betaMax = betaMin
	

	for (i in 1:length(R)) {
		# since beta is a percentage of each reading, each reading has
		# a corresponding beta.  Need to loop to find them all	
		
		betaMax[i] = R[i] + (0.9/100*(R[i]+0.2)+0.2)
		betaMin[i] = R[i] - (0.9/100*(R[i]+0.2)+0.2)
		# https://www.es.co.th/Schemetic/PDF/FLUKE-114-117.PDF
		# https://www.fluke.com/en-us/learn/blog/digital-multimeters/accuracy-precision
	
	}
	
	return(list(betaMax,betaMin))

}




beta = find_beta(R)
P = 0.95

b = c(beta[[1]],beta[[2]])
b = sd(b)

par_mean = mean_p_boot(R,b,P)
npar_mean = mean_np_boot(R,b,P)

par_mean
npar_mean


### 4.  Evaluate 95% uncertainty interval for the standard deviation of the population
	# use nonparametric MC approach.  Use parametric MC for systematic errors


par_sd = sd_p_boot(R,b,P)
npar_sd = sd_np_boot(R,b,P)

par_sd
npar_sd