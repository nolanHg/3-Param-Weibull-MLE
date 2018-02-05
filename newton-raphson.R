options(digits = 21)

f_orig <- function(x) {
	return(exp(2 * x) - 2)
}

f_deriv <- function(x) {
	return(2 * exp(2 * x))
}

iters = 100 
estimate = 5 
for (i in c(0:iters)) {
	estimate = estimate - (f_orig(estimate) / f_deriv(estimate))
	print(estimate)
}
