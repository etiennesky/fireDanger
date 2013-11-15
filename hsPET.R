hsPET <-
function(date, lat=42, t.min, t.max, ret.rad = FALSE) {
      d <- as.POSIXlt(date)
      lat <- lat * (pi / 180)
      tn <- t.min
      tx <- t.max
      rr <- ret.rad
      varnames <- c("d", "tn", "tx")
      ind <- c()
      for (i in 1:length(varnames)) {
            vec <- get(varnames[i])
            if (any(is.na(vec))) {
                  a <- which(is.na(vec))
                  ind <- c(ind, a)
            }
      }
      if (length(ind) > 0) {
            warning("Some missing values were removed from the time series before computation")
            ind <- unique(ind)
            for (i in 1:length(varnames)) {
                  assign(varnames[i], get(varnames[i])[-ind])
            }
      }
      J <- as.numeric(format(d, format = "%j"))
      Gsc <- 0.082
      ds <- 0.409 * sin(2 * pi * J / 365 - 1.39)
      ws <- acos(-tan(lat) * tan(ds))
      dr <- 1 + 0.033 * cos(2 * pi * J / 365)
      Ra <- (24 * 60 * .082 * dr * (ws * sin(lat) * sin(ds) + cos(lat) * cos(ds) * sin(ws))) / pi
      hs <- .0023 * (Ra / 2.45) * ((tx + tn) / 2 + 17.8) * sqrt(tx - tn)
      if (isTRUE(rr)) { 
            return(list("date" = d, "Rad" = Ra, "HSevap" = hs))
      }
      else {
            return(list("date" = d, "HSevap" = hs))
      }
}
