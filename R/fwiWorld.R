#' @title World Fire Weather Index 
#' 
#' @description Implementation of the Canadian Fire Weather Index System 
#' 
#' @param dataset A character string indicating the database to be accessed (partial matching is enabled). 
#' Currently accepted values are "WFDEI, "CFSv2_seasonal" and "S4_seasonal_15members".
#' @param dictionary A logical flag indicating if a dictionary is used for variable homogenization. Default (strongly recommended) 
#' is set to TRUE, meaning that the function will internally perform the necessary homogenization steps to return the standard 
#' variables defined in the vocabulary (e.g. variable transformation, deaccumulation...). See details on data homogenization.
#' @param members A vector of integers indicating the members to be loaded. Default to NULL, which loads the default members 
#' (see details for the particularities concerning the CFSv2 dataset). For instance, members=1:5 will retrieve the first five 
#' members of the hindcast. Discontinuous member selection (e.g. members=c(1,5,7)) is allowed. If the requested dataset is not 
#' a forecast or the requested variable is static (e.g. orography) it will be ignored.
#' @param season An integer vector specifying the desired season (in months, January = 1 ..., December = 12). Options include 
#' one to several (contiguous) months. For full year selections (not possible for all datasets, e.g. seasonal forecasts), 
#' the argument value must be set to season = 1:12. If the requested variable is static (e.g. orography) it will be ignored. 
#' See details on the definition of temporal slices.
#' @param years Optional vector of years to select. Default to all available years. If the requested variable is static (e.g. orography) 
#' it will be ignored. See details on the definition of temporal slices.
#' @param leadMonth Integer value indicating the lead forecast time, relative to the first month of season. Note that leadMonth=1 for 
#' season=1 (January) corresponds to the December initialization. Default to 1 (i.e., 1 lead month forecast). If the dataset is not a 
#' forecast or the requested variable is static (e.g. orography) it will be ignored. A message will be printed on screen in the first 
#' case if its value is different from NULL. See details on initialization times.
#' @param return.all Logical. Should all components of the FWI system be returned?. 
#' Default to FALSE, indicating that only FWI is returned.
#' @param init.pars A numeric vector of length 3 with the initialization values for the
#'  FFMC, DMC and DC components, in this order. Default values as proposed by van Wagner (1987).
#' @param parallel Logical. Should parallel execution be used?
#' @param max.ncores Integer. Upper bound for user-defined number of cores.
#' @param ncores Integer number of cores used in parallel computation. Self-selected number of
#'  cores is used when \code{ncpus = NULL} (the default), or when \code{maxcores} exceeds the default \code{ncores} value.
#' 
#' @return A vector of the same length as the input vectors (minus possible missing observations),
#' containing the requested components of the FWI system (either all or just FWI). See details.
#' 
#' @section Daylength adjustment factors: 
#' By default, the function applies the original FWI daylength adjustment factors for DC and DMC (van Wagner 1987),
#'  although it is possible to adjust them by as a function of latitude through the argument \code{lat}.
#' The reference values used for each latitudinal range are those indicated in p.71 and Tables A3.1 and A3.2 (p69) in
#' Lawson and Armitage (2008).
#' 
#' @references
#' Lawson, B.D. & Armitage, O.B., 2008. Weather guide for the Canadian Forest Fire Danger Rating System. Northern Forestry Centre, Edmonton (Canada).
#' 
#' van Wagner, C.E., 1987. Development and structure of the Canadian Forest Fire Weather Index (Forestry Tech. Rep. No. 35). Canadian Forestry Service, Ottawa, Canada.
#' 
#' van Wagner, C.E., Pickett, T.L., 1985. Equations and FORTRAN program for the Canadian forest fire weather index system (Forestry Tech. Rep. No. 33). Canadian Forestry Service, Ottawa, Canada.
#' 
#' @author J. Bedia \& M.Iturbide, partially based on the original FORTRAN code by van Wagner and Pickett (1985)
#' @export
#' @importFrom abind abind
#' @importFrom downscaleR makeMultiGrid
#' @importFrom downscaleR subsetGrid
#' @import loadeR.ECOMS



fwiWorld <- function(dataset = "WFDEI",
                     dictionary = TRUE, 
                     years = NULL,
                     season = NULL, 
                     lonLim = NULL,
                     latLim = NULL,
                     leadMonth = 0,
                     mask = NULL,
                     return.all = FALSE, 
                     init.pars = c(85, 6, 15),
                     parallel = FALSE,
                     max.ncores = 16,
                     ncores = NULL){
      x <-  c(-90, -80, -70, -60, -50, -40, -30, -20, -10, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90)
      
      if (is.null(latLim)) latLim <-c(-90, 90)
      if (is.null(lonLim)) lonLim <- c(-180, 180)
      latind <- findInterval(latLim, x)[1] : findInterval(latLim, x)[2]
      if(x[latind[length(latind)]] < latLim[2]) latind[3] <- latind[2]+1
      x <- x[latind]
      lats <- seq(min(x)+5, max(x)-5, 10) 
      if(x[length(x)] > latLim[2]) x[length(x)] <- latLim[2]
      if(x[1] < latLim[1]) x[1] <- latLim[1]
      
      message("[", Sys.time(), "] Pre-processing...")
      suppressMessages(
      temp <- loadECOMS(dataset,  latLim = latLim, lonLim = lonLim, var = "tas", dictionary = dictionary, 
                        leadMonth = leadMonth,
                        season = season[-1], years = years, time = "DD",
                        aggr.d = "mean")
      )
      latdim <- which(downscaleR:::getDim(temp) == "lat")
      dimNames <- attr(temp$Data, "dimensions")
      a <- list()
      for(i in 1:(length(x)-1)){
            latLimchunk <- c(x[i], x[i+1])
            lat <- lats[i]
            message("------------processing latitude chunk ", i, " out of ", length(x)-1, "-------------")     
            Tm <- loadECOMS(dataset,  latLim = latLimchunk, lonLim = lonLim, var = "tas", dictionary = dictionary, 
                          season = season, years = years, leadMonth = leadMonth, time = "DD",
                          aggr.d = "mean")
            H <- loadECOMS(dataset,  latLim = latLimchunk, lonLim = lonLim, var = "hurs", dictionary = dictionary, 
                         season = season, years = years, leadMonth = leadMonth, time = "DD",
                         aggr.d = "min")
            r <- loadECOMS(dataset,  latLim = latLimchunk,  lonLim = lonLim, var = "tp", dictionary = dictionary, 
                         season = season, years = years, time = "DD",
                         aggr.d = "sum")
            W <- loadECOMS(dataset,  latLim = latLimchunk, lonLim = lonLim, var = "wss", dictionary = dictionary, 
                         season = season, years = years, time = "DD",
                         aggr.d = "mean")
            W$Data <- W$Data*3.6
            attr(W$Data, "dimensions") <- attr(r$Data, "dimensions") 
            multigrid <- makeMultiGrid(Tm, H, r, W)
            if(is.null(mask) & dataset == "WFDEI"){
                  msk <- W
                  msk$Data <- msk$Data[1,,]
                  msk$Data[which(!is.na(msk$Data))] <- 100
                  msk$Data[which(is.na(msk$Data))] <- 0
                  attr(msk$Data, "dimensions") <- c("lat", "lon")
            }else if(!is.null(mask)){
                  msk <- subsetGrid(mask, latLim = latLimchunk, lonLim = lonLim, outside = T)
            }else{
                  message("The use of a land mask is recommended")
            }
            if(i != (length(x)-1)){
                  xx <- dim(Tm$Data)[latdim]-1
            }else{
                  xx <- dim(Tm$Data)[latdim]
            }
            Tm <- H <- r <- W <- NULL
            o <- lapply(1:length(years), function(x){
                  multigrid.y <- subsetGrid(multigrid, years = years[x])
                  fwi(multigrid = multigrid.y, mask = msk, lat = lat, return.all = return.all, 
                        parallel = parallel, init.pars = init.pars, 
                        max.ncores = max.cores, ncores = ncores)$Data[,,1:xx,,drop=FALSE]
            })
            o.full <-  unname(do.call("abind", list(o, along = 2)))  
            months <- as.integer(substr(multigrid$Dates[[1]]$start, start = 6, stop = 7))
            multigrid.y <- NULL
            month.ind <- which(months == months[1])
            a[[i]] <- o.full[,-month.ind,,]
      }
      out <- unname(do.call("abind", list(a, along = latdim)))
      temp$Data <- out
      attr(temp$Data, "dimensions") <- dimNames
      temp$Variable$varName <- "fwi"
      attr(temp$Variable, "use_dictionary") <- FALSE
      attr(temp$Variable, "description") <- "Fire Weather Index"
      attr(temp$Variable, "units") <-  "none"
      attr(temp$Variable, "longname") <- "Fire Weather Index"
#       if(length(years) > 1){
#             yname <- paste(years[1], "_", years[length(years)], sep = "")
#       }else{
#             yname <- years
#       }
return(temp)
}