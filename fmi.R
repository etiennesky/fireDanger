fmi <-
function (t, rh) {
	t <- t
	r <- rh
	if (length(t) != length(r)) {
		stop("Input time series of differing lengths")
	}
	if (any(is.na(t) | is.na(r))) {
		unique(c(which(is.na(t)), which(is.na(r)))) -> na
		t <- t[-na] 
		r <- r[-na] 
	}
	fmi <- rep(NA, length(t))
	for (i in 1:length(t)) {
		fmi[i] <- 10 - 0.25 * (t[i] - r[i])
	}
	return(fmi)
}
