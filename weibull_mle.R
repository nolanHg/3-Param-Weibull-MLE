###################################################################
#								  #
# file:		weibull_mle.R					  #
# authors:	Nolan Gagnon, Puspanjali Subudhi 		  #
# version:	V1.0						  #
# date:		02/04/2018					  #
#								  #
# brief:	Calculates the maximum likelihood estimates       #
#		for the parameters of a three-parameter           #
#		Weibull distribution using the method described   #
#		in Panchang and Gupta (1989) [henceforth P&G].	  #
#								  #
###################################################################


###################################################################
#			   IMPORTS		                  #
###################################################################
library(FAdist) # Provides access to the rweibull3() function.


###################################################################
#			   R OPTIONS				  #
###################################################################
options(digits = 5) # Display numbers with 5 digits.


###################################################################
#			   FUNCTIONS				  #
###################################################################

##
## brief: Compute the next estimate for the shape parameter using Newton-Raphson.
## param y: Modified Weibull sample x_j - u_i.
## param prev_w: Previous estimate for the shape parameter.
## returns: Next estimate for the shape parameter.
##
get_next_w <- function(y, prev_w) 
{
	## Sums used in Newton-Raphson method.	
	S1 <- sum(log(y)) 
	S2 <- sum(y ^ prev_w) 
	S3 <- sum((y ^ prev_w) * log(y)) 
	S4 <- sum(((log(y))^2) * (y^prev_w)) 
	
	## Next Newton-Raphson estimate; found using equation 2.1 in P&G.
	## For reference:  Next Newton-Raphson estimate is x_(n+1) = x_n - f(x_n) / f'(x_n).
	next_w <- prev_w + ((1 / prev_w) + (S1 / length(y)) - (S3 / S2)) / ((1 / (prev_w^2)) + S4 / S2 - (S3 / S2)^2)
	return(next_w)	
}

##
## brief: Compute the maximized log-likelihood function using u_i, v_i, w_i,
##	   found using the method described in P&G.
## param v: Scale parameter.
## param w: Shape parameter.
## param y: Modified Weibull sample x_j - u_i
## returns: L(u_i, v_i, w_i), i.e., the maximized log-likelihood function for u_i, v_i, w_i.
##
max_log_likelihood <- function(v, w, y) 
{
	S <- sum(log(y)) 
	L <- length(y) * (log(w / v) - 1) + (w - 1) * S
	return(L)
}	

##
## brief: Compute the maximized log-likelihood function using u_i = x_min,
##	  v_i = (1 / N) * SUMMATION(x_j - x_min), w_i = 1.
## param v: Scale parameter with value equal to v_i, mentioned in the brief.
## param n: Sample size.
##
max_log_likelihood_alt <- function(v, n)
{
	L <- n * (-log(v) - 1)
	return(L)
}

##
## brief: Compute MLEs for the parameters of a three-parameter Weibull distribution 
##	  using the method in P&G.
## param data: Sample from Weibull distribution. 
## param epsilon: Small number to prevent argument of logarithm
##	          in equation 1.3 (from P&G) from being 0.
## param u_divs: Number regions to divide u-space into.
## param NR_iters: Number of iterations to use in Newton-Raphson.
##
get_weibull_mles <- function(data, epsilon, u_divs, NR_iters) 
{
	## Sample size.
	N <- length(data)

	## Get the minimum value from the sample (x_min in the paper).
	data_min <- min(data)

	## Generate divided u-space
	u_space <- seq(0, data_min - epsilon, len = u_divs)
	
	## Compute estimate for w_1 (when u_i = 0).
	S <- sum(log(data))
	w_estimate <- N / (N * log(max(data)) - S)
	
	## Allocate memory space for v-parameter and w-parameter vectors.
	v_vect <- rep(NA, length(u_space))
	w_vect <- rep(NA, length(u_space))
	
	## Number of iterations to perform in Newton-Raphson.	
	iters <- NR_iters
	
	## Build the v-parameter and w-parameter vectors.
	for (i in 1 : length(u_space)) {
		
		## Initial Newton-Raphson estimate for w_i (when i != 1).
		if (i != 1) {	
			w_estimate_prev <- w_estimate
			w_estimate <- w_estimate - (w_estimate_prev - w_estimate)
		}
	
		## Compute modified Weibull sample.
		y <- data - u_space[i]
		
		## Calculate final Newton-Raphson estimate for shape parameter w_i.
		for (j in 1 : iters) { w_estimate <- get_next_w(y, w_estimate) }
		
		## Populate v-parameter and w-parameter vectors.	
		if (i == length(u_space)) {
			## When u_space[i] == data_min, w_i and v_i
			## are given by equation 2.3 from P&G.
			v_vect[i] <- (1 / N) * sum(data - data_min)
			w_vect[i] <- 1
		} else {
			## w_i is given by the Newton-Raphson estimate above.
			w_vect[i] <- w_estimate
			
			## v_i is calculated using w_i and u_space[i] in equation
			## 1.4d from P&G.	
			v_vect[i] <- (1 / N) * sum(y ^ w_estimate)
		}
	}
	
	## Allocate memory space for maximized log-likelihood vector (L).
	L <- rep(NA, length(u_space)) 

	## Build the L-vector.
	for (i in 1 : length(u_space)) {

		y <- data - u_space[i]
		v <- v_vect[i] 
		w <- w_vect[i] 
		
		## Use equation 2.3 from P&G when u_space[i] = data_min.	
		if (i == length(u_space)) { L[i] <- max_log_likelihood_alt(v, N) }
		
		## Otherwise, use equation 2.2 from P&G.
		else { L[i] <- max_log_likelihood(v, w, y) }
	}

	## Find the index of the largest value in the maximized log-likelihood vector.
	index_of_max <- which(L == max(L))
	
	## Get the MLEs for u, v, and w.	
	u_mle <- u_space[index_of_max]
	v_mle <- v_vect[index_of_max]
	w_mle <- w_vect[index_of_max]
	mles <- c(u_mle, v_mle, w_mle)
	
	ret_vect <- vector(mode = "list", 2)	
	
	ret_vect[[1]] <- mles
	ret_vect[[2]] <- c(u_space, L)

	## Return the MLEs and (u_space, L) in a list of vectors.
	return(ret_vect)
}

#plot_mllf <- function(data

w_param <- 1.5 
v_param <- 1 
u_param <- 10	
sample_size <- 100
sample_name <- paste("Random Sample with N = ", sample_size, " from Weibull(u = ", u_param, ", v = ", v_param, ", w = ", w_param, ")", sep = "")
w3samp <- rweibull3(sample_size, shape = w_param, scale = v_param, thres = u_param)

W <- get_weibull_mles(w3samp, 10e-6, 500, 5)
MLEs <- W[[1]]
u_and_L <- W[[2]]
print(MLEs)

#jpeg(file = paste(gsub(" ", "_", sample_name), ".jpeg", sep =""))
#plot(u_space, L, main = sample_name, xlim = c(min(u_space), max(u_space)), ylim = c(min(L), max(L)), type = "l", sub =  paste("u = ", u_space[index_of_max], ", v = ", v_vect[index_of_max], ", w = ", w_vect[index_of_max], sep = ""))
#dev.off()
