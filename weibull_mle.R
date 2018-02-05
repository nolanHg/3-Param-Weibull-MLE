library(FAdist)

options(digits = 10)

get_next_w <- function(y, prev_w) 
{
	S1 <- sum(log(y)) 

	S2 <- sum(y ^ prev_w) 
	S3 <- sum((y ^ prev_w) * log(y)) 
	S4 <- sum(((log(y))^2) * (y^prev_w)) 

	next_w <- prev_w + ((1 / prev_w) + (S1 / length(y)) - (S3 / S2)) / ((1 / (prev_w^2)) + S4 / S2 - (S3 / S2)^2)
	return(next_w)	
}

log_likelihood <- function(v, w, y) 
{
	S <- sum(log(y)) 

	L <- length(y) * (log(w / v) - 1) + (w - 1) * S

	return(L)

}	

## Small positive number used to prevent
## the argument of the log function in Equation
## 1.3 from Panchang and Gupta (1989) from 
## being zero.
epsilon <- 10^-3

## Generate a random sample of size N from Weibull(u = 50, v = 3, w = 2)
## x_j sample
#w3samp <- rweibull3(N, shape = 2.3, scale = 3, thres = 2.8)
#w3samp <- c(35.0, 34.7, 30.9, 29.0, 28.1, 24.9, 24.0, 23.9, 23.3, 22.6, 22.4, 19.8, 19.8, 19.4, 19.0, 17.6, 16.5, 15.9, 13.3, 12.3, 12.0, 12.0)

#w3samp <- c(10.0805, 10.0990, 10.2757, 10.6545, 10.6883, 11.0666, 11.2083, 11.2558, 11.8761, 12.2103)

w3samp <- c(3.1, 4.6, 5.6, 6.8)

N <- length(w3samp)

## Get the minimum value from the sample (x_min in the paper)
w3samp_min <- min(w3samp)

## Generate divided u-space
u_space <- seq(0, w3samp_min - epsilon, len = 1000)

iters <- 4 
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

	w_vect[i] <- w_estimate
	v_vect[i] <- (1 / N) * sum(y ^ w_estimate)
}

L <- rep(NA, length(u_space)) 
for (i in 1 : length(u_space)) {

	y <- w3samp - u_space[i]
	v <- v_vect[i] 
	w <- w_vect[i] 

	L[i] <- log_likelihood(v, w, y)
}

index_of_max <- which(L == max(L))
paste("u =", u_space[index_of_max], "v =", v_vect[index_of_max], "w =", w_vect[index_of_max], sep = " ")

jpeg(file = "ll.jpeg")
plot(u_space, L, xlim = c(0, 3.5), ylim = c(-6.84, -6.79), type = "l")
dev.off()
