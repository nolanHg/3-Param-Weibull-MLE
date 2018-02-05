##################################################################
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
## returns: L(u_i = x_min, v_i = (1 / N) * SUMMATION(x_j - x_min), w_i = 1).
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
## returns: A list containing MLEs, u-space vector, L-vector, v-parameter vector,
##	    and w-parameter vector.
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
	
	ret_vect <- vector(mode = "list", 5)	
	
	ret_vect[[1]] <- mles
	ret_vect[[2]] <- u_space
	ret_vect[[3]] <- v_vect
	ret_vect[[4]] <- w_vect
	ret_vect[[5]] <- L

	## Return the MLEs, u_space, v_vect, w_vect, and L in a list of vectors.
	return(ret_vect)
}

plot_mllf_u <- function(u, L, title, subtitle, x_lims, y_lims)
{
	pdf(file = paste("./weibull_plots/", gsub(" ", "_", title), "_u", ".pdf", sep = ""))
	par(bg = 'lightgray')
	plot(u, L, main = title, sub = subtitle, xlim = x_lims, ylim = y_lims, type = "l", lwd = 2,
	     col = "red", col.main = "black", col.sub = "red", ylab = "Maximized log-likelihood",
	     fg = "blue")
	grid(col = "black")
	dev.off()	
}

plot_mllf_v <- function(v, L, title, subtitle, x_lims, y_lims)
{
	pdf(file = paste("./weibull_plots/", gsub(" ", "_", title), "_v", ".pdf", sep = ""))
	par(bg = 'lightgray')
	plot(v, L, main = title, sub = subtitle, xlim = x_lims, ylim = y_lims, type = "l", lwd = 2,
	     col = "red", col.main = "black", col.sub = "red", ylab = "Maximized log-likelihood",
	     fg = "blue")
	grid(col = "black")
	dev.off()	
}

plot_mllf_w <- function(w, L, title, subtitle, x_lims, y_lims)
{
	pdf(file = paste("./weibull_plots/", gsub(" ", "_", title), "_w", ".pdf", sep = ""))
	par(bg = 'lightgray')
	plot(w, L, main = title, sub = subtitle, xlim = x_lims, ylim = y_lims, type = "l", lwd = 2,
	     col = "red", col.main = "black", col.sub = "red", ylab = "Maximized log-likelihood",
	     fg = "blue")
	grid(col = "black")
	dev.off()	
}

###################################################################
#			   VISUALIZATION	                  #
###################################################################
args <- commandArgs(trailingOnly = TRUE)

if (args[1] == "SGEN") {
	u_cmdln <- as.numeric(args[2])
	v_cmdln <- as.numeric(args[3])
	w_cmdln <- as.numeric(args[4])
	ssize_cmdln <- as.numeric(args[5])
	udivs_cmdln <- as.numeric(args[6])

	w3samp <- rweibull3(ssize_cmdln, shape = w_cmdln, scale = v_cmdln, thres = u_cmdln)
	
	calcd_weib_prods <- get_weibull_mles(w3samp, 10e-6, udivs_cmdln, 5)	

	MLEs <- calcd_weib_prods[[1]]
	u <- calcd_weib_prods[[2]]
	v <- calcd_weib_prods[[3]]
	w <- calcd_weib_prods[[4]]
	L <- calcd_weib_prods[[5]]

	plot_mllf_u(u, L, paste("Random sample from Weib(u = ", u_cmdln, ", v = ", v_cmdln, ", w = ", w_cmdln, ")",
		                 sep = ""), paste("u = ", MLEs[1], ", v = ", MLEs[2], ", w = ", MLEs[3], sep = ""),
		     c(min(u), max(u)), c(min(L), max(L)))
} else {

	data_files <- list.files(paste(getwd(), "/weibull_data", sep = ""))
		
	for (file in data_files) {

		w3samp <- as.numeric(scan(file = paste(getwd(), "/weibull_data/", file, sep = ""), 
					  what = character(), sep = ","))

		## Calculate MLEs, u-space, v-space, w-space, and L.
		calcd_weib_prods <- get_weibull_mles(w3samp, 10e-6, 100, 5)

		## Extract calculated products from calcd_weib_prods.
		MLEs <- calcd_weib_prods[[1]]
		u <- calcd_weib_prods[[2]]
		v <- calcd_weib_prods[[3]]
		w <- calcd_weib_prods[[4]]
		L <- calcd_weib_prods[[5]]
		
		if (file == "Rockette.csv") {
			x_lims <- c(0, 3)
			y_lims <- c(-6.84, -6.79)
		} else if (file == "PA.csv") {
			x_lims <- c(0, 12)
			y_lims <- c(-73.5, -71.5)
		} else if (file == "AC.csv") {
			x_lims <- c(0, 10)
			y_lims <- c(-12, -8)
		} else {
			x_lims <- c(0, 0.6)
			y_lims <- c(-24, -14)
		}

		plot_mllf_u(u, L, gsub(".csv", " Data", file), 
			     paste("u = ", MLEs[1], ", v = ", MLEs[2], ", w = ", MLEs[3], sep = ""), 
			     x_lims, y_lims)

		plot_mllf_v(v, L, gsub(".csv", " Data", file), 
			     paste("u = ", MLEs[1], ", v = ", MLEs[2], ", w = ", MLEs[3], sep = ""), 
			     c(min(v), max(v)), c(min(L), max(L)))

		plot_mllf_w(w, L, gsub(".csv", " Data", file), 
			     paste("u = ", MLEs[1], ", v = ", MLEs[2], ", w = ", MLEs[3], sep = ""), 
			     c(min(w), max(w)), c(min(L), max(L)))
	}
}
