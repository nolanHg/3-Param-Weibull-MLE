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
#		in Panchang and Gupta (1989).			  #
#								  #
###################################################################


###################################################################
#			   IMPORTS		                  #
###################################################################
library(FAdist) # Provides access to the rweibull3() function


###################################################################
#			   R OPTIONS				  #
###################################################################
options(digits = 5) # Display numbers with 5 digits


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
	
	S1 <- sum(log(y)) 
	S2 <- sum(y ^ prev_w) 
	S3 <- sum((y ^ prev_w) * log(y)) 
	S4 <- sum(((log(y))^2) * (y^prev_w)) 
	
	## Next Newton-Raphson estimate; found using equation 2.1 in Panchang and Gupta (1989)
	## For reference:  Next Newton-Raphson estimate is x_(n+1) = x_n - f(x_n) / f'(x_n).
	next_w <- prev_w + ((1 / prev_w) + (S1 / length(y)) - (S3 / S2)) / ((1 / (prev_w^2)) + S4 / S2 - (S3 / S2)^2)
	return(next_w)	
}

##
## brief: Compute the maximized log-likelihood function using u_i, v_i, w_i,
##	   found using the method described in Panchang and Gupta (1989).
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

## Small positive number used to prevent
## the argument of the log function in Equation
## 1.3 from Panchang and Gupta (1989) from 
## being zero.
epsilon <- 10e-6

## Generate a random sample of size N from Weibull(u = 50, v = 3, w = 2)
## x_j sample
w_param <- 1.5 
v_param <- 1 
u_param <- 10	
sample_size <- 100
sample_name <- paste("Random Sample with N = ", sample_size, " from Weibull(u = ", u_param, ", v = ", v_param, ", w = ", w_param, ")", sep = "")
w3samp <- rweibull3(sample_size, shape = w_param, scale = v_param, thres = u_param)

## Sample from Adatia and Chan (1985)
#sample_name <- "Adatia and Chan (1985)"
#w3samp <- c(10.0805, 10.0990, 10.2757, 10.6545, 10.6883, 11.0666, 11.2083, 11.2558, 11.8761, 12.2103)

## Sample from Rockette et al. (1974)
#sample_name <- "Rockette et al. (1974)"
#w3samp <- c(3.1, 4.6, 5.6, 6.8)

## Sample from Petruaskas and Aagaard (1971)
#sample_name <- "Petruaskas and Aagaard (1972)"
#w3samp <- c(35, 34.7, 30.9, 29.0, 28.1, 24.9, 24, 23.9, 23.3, 22.6, 22.4, 19.8, 19.8, 19.4, 19, 17.6, 16.5, 15.9, 13.3, 12.3, 12, 12)

## Sample from Smith and Naylor (1987)
#sample_name <- "Smith and Naylor (1987)"
#w3samp <- c(0.55, 0.74, 0.77, 0.81, 0.84, 0.93, 1.04, 1.11, 1.13, 1.24, 1.25, 1.27, 1.28, 1.29, 1.30, 1.36, 1.39, 1.42, 1.48, 1.48, 1.49, 1.49, 1.50, 1.50, 1.51, 1.52, 1.53, 1.54, 1.55, 1.55, 1.59, 1.59, 1.60, 1.61, 1.61, 1.61, 1.61, 1.62, 1.62, 1.63, 1.64, 1.66, 1.66, 1.66, 1.67, 1.68, 1.68, 1.69, 1.70, 1.70, 1.73, 1.76, 1.76, 1.77, 1.78, 1.81, 1.82, 1.84, 1.84, 1.89, 2.00, 2.01, 2.24)

N <- length(w3samp)

## Get the minimum value from the sample (x_min in the paper)
w3samp_min <- min(w3samp)

## Generate divided u-space
u_space <- seq(0, w3samp_min - epsilon, len = 100)

iters <- 10
sum <- 0
for (i in 1 : N) {
	sum <- sum + log(w3samp[i])
}
w_estimate <- N / (N * log(max(w3samp)) - sum)

v_vect <- rep(NA, length(u_space))
w_vect <- rep(NA, length(u_space))

for (i in 1 : length(u_space)) {
	
	if (i != 1) {	
		w_estimate_prev <- w_estimate
		w_estimate <- w_estimate - (w_estimate_prev - w_estimate)
	}

	y <- w3samp - u_space[i]

	for (j in 1 : iters) {	
		w_estimate <- get_next_w(y, w_estimate)
	}
	
	
	if (i == length(u_space)) {
 		v_vect[i] <- (1 / N) * sum(w3samp - min(w3samp))
		w_vect[i] <- 1
	} else {
		w_vect[i] <- w_estimate
		v_vect[i] <- (1 / N) * sum(y ^ w_estimate)
	}
}

L <- rep(NA, length(u_space)) 
for (i in 1 : length(u_space)) {

	y <- w3samp - u_space[i]
	v <- v_vect[i] 
	w <- w_vect[i] 
	
	if (i == length(u_space)) {
		L[i] <- max_log_likelihood_alt(v, N)
	} else {
		L[i] <- max_log_likelihood(v, w, y)
	}
}

index_of_max <- which(L == max(L))
paste("u =", u_space[index_of_max], "v =", v_vect[index_of_max], "w =", w_vect[index_of_max], sep = " ")

jpeg(file = paste(gsub(" ", "_", sample_name), ".jpeg", sep =""))
plot(u_space, L, main = sample_name, xlim = c(min(u_space), max(u_space)), ylim = c(min(L), max(L)), type = "l", sub =  paste("u = ", u_space[index_of_max], ", v = ", v_vect[index_of_max], ", w = ", w_vect[index_of_max], sep = ""))
dev.off()
