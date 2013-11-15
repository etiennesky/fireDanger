angstronIndex <-
function (t, rh) {
	t <- t
	r <- rh
	if (length(t) != length(r)) {
		stop("Input time series of differing lengths")
	}
	if (any(is.na(t) | is.na(r))) {
		na <- unique(c(which(is.na(t)), which(is.na(r))))
		t <- t[-na]
		r <- r[-na]
	}
	AI <- rep(NA, length(t))
	for (i in 1:length(t)) {
		AI[i] <- (r[i] / 20) + ((27 - t[i])/ 10)
	}
	return(AI)
}
