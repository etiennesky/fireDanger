wfwi_obs <- function(dataset = "WFDEI",
                     aggr = list(FUN = NULL),
                     destdir = getwd(), 
                     dictionary = TRUE, 
                     years = NULL,
                     season = NULL, 
                     lonLim = NULL,
                     latLim = NULL,
                     return.all = FALSE, 
                     init.pars = c(85, 6, 15),
                     parallel = FALSE,
                     max.ncores = 16,
                     ncores = NULL){
      x <-  c(-90, -80, -70, -60, -50, -40, -30, -20, -10, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90)
      
      if (is.null(latLim)) latLim <-c(-90, 90)
      if (is.null(lonLim)) lonLim <- c(-180, 180)
      latindmin <- which(x < min(latLim))
      if(length(latindmin)!=0) x <- x[-latindmin]
      latindmax <- which(x > max(latLim))
      if(length(latindmax)!=0)  x <- x[-latindmax]
      lats <- seq(min(x)+5, max(x)-5, 10) 
      
      temp <- loadECOMS(dataset,  latLim = latLim, lonLim = lonLim, var = "tas", dictionary = dictionary, 
                        season = season, years = years, time = "DD",
                        aggr.d = "mean")
      a <- list()
      for(i in 1:(length(x)-1)){
            latLimchunk <- c(x[i], x[i+1])
            lat <- lats[i]
       print(i)     
      Tm.obs <- loadECOMS(dataset,  latLim = latLimchunk, lonLim = lonLim, var = "tas", dictionary = dictionary, 
                          season = season, years = years, leadMonth = leadMonth, time = "DD",
                          aggr.d = "mean")
      H.obs <- loadECOMS(dataset,  latLim = latLimchunk, lonLim = lonLim, var = "hurs", dictionary = dictionary, 
                         season = season, years = years, leadMonth = leadMonth, time = "DD",
                         aggr.d = "min")
      r.obs <- loadECOMS(dataset,  latLim = latLimchunk,  lonLim = lonLim, var = "tp", dictionary = dictionary, 
                         season = season, years = years, time = "DD",
                         aggr.d = "sum")
      W.obs <- loadECOMS(dataset,  latLim = latLimchunk, lonLim = lonLim, var = "wss", dictionary = dictionary, 
                         season = season, years = years, time = "DD",
                         aggr.d = "mean")
      W.obs$Data <- W.obs$Data*3.6
      attr(W.obs$Data, "dimensions") <- attr(r.obs$Data, "dimensions") 
      multigrid.obs <- makeMultiGrid(Tm.obs, H.obs, r.obs, W.obs)
      obsmask <- W.obs
      obsmask$Data <- obsmask$Data[1,,]
      obsmask$Data[which(!is.na(obsmask$Data))] <- 100
      obsmask$Data[which(is.na(obsmask$Data))] <- 0
      attr(obsmask$Data, "dimensions") <- c("lat", "lon")
      xx <- dim(Tm.obs$Data)[2]-1
      Tm.obs <- H.obs <- r.obs <- W.obs <- NULL
      o <- fwi(multigrid.obs, mask = obsmask, lat = lat, return.all = return.all, 
               parallel = parallel, init.pars = init.pars, 
               max.ncores = max.cores, ncores = ncores)$Data[,,1:xx,,drop=FALSE]
      print("hello")
      months <- as.integer(substr(multigrid.obs$Dates[[1]]$start, start = 6, stop = 7))
      multigrid.obs <- NULL
      month.ind <- which(months == months[1])
      a[[i]] <- o[,-month.ind,,]
}
#       c.x <- unname(do.call("abind", list(coords.x)))
#       c.y <- unname(do.call("abind", list(coords.y)))
#       xy <- cbind(c.x, c.y)

      out <- unname(do.call("abind", list(a, along = 2)))
      temp$Data <- out
      attr(temp$Data, "dimensions") <- c("time", "lat", "lon")
      temp$Variable$varName <- "fwi"
      attr(temp$Variable, "use_dictionary") <- FALSE
      attr(temp$Variable, "description") <- "Fire Weather Index"
      attr(temp$Variable, "units") <-  "none"
      attr(temp$Variable, "longname") <- "Fire Weather Index"
#       write.csv(out, file = paste(destdir, "/wfi_", dataset, "_", season[2], "_", season[length(season)], ".csv", sep=""))
#       write.csv(xy, file = paste(destdir, "/coordinates_", dataset,"_", season[2], "_", season[length(season)], ".csv", sep=""))
save(temp, file = paste(destdir, "/wfi_", dataset, "_", years, "_", season[2], "_", season[length(season)], ".rda", sep=""))
}