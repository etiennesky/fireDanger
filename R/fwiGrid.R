#' @title Fire Weather Index (FWI) from multigrid
#' 
#' @description A wrapper of function \code{\link{fwi1D}} to handle multigrids
#' 
#' @param multigrid A multigrid of the variables needed to compute the FWI. See the Warning section
#' @param mask Optional. Grid of the land mask to be applied to the data.
#' @param return.all Logical. Should all components of the FWI system be returned?. 
#' Default to FALSE, indicating that only FWI is returned.
#' @param init.pars A numeric vector of length 3 with the initialization values for the
#'  FFMC, DMC and DC components, in this order. Default values as proposed by van Wagner (1987).
#' @param parallel Logical. Should parallel execution be used?
#' @param max.ncores Integer. Upper bound for user-defined number of cores.
#' @param ncores Integer number of cores used in parallel computation. Self-selected number of
#'  cores is used when \code{ncpus = NULL} (the default), or when \code{maxcores} exceeds the default \code{ncores} value.
#' 
#' @return A grid, containing the requested components of the FWI system (either all or just FWI). See details.
#' 
#' @details 
#' 
#' In order to efficiently handle large datasets, the function internally splits the spatial domain in convenient latitudinal chunks.  
#' 
#' @section Daylength adjustment factors: 
#' In this case, daylength adjustment factors are automatically set by the funciton according to the
#'  reference values indicated in p.71 and Tables A3.1 and A3.2 (p69) in Lawson and Armitage (2008).
#' 
#' @section Warning:
#' 
#' The variables composing the input multigrid need to have standard names, as defined by the dictionary
#'  (their names are stored in the \code{multigrid$Variable$varName} component).
#' These are: \code{"tas"} for temperature, \code{"tp"} for precipitation, \code{"wss"} for windspeed. In the case of relative humidity,
#' either \code{"hurs"} or \code{"hursmin"} are accepted, the latter in case of FWI calculations according to the \dQuote{proxy} version
#' described in Bedia \emph{et al} 2014.
#' 
#' Be aware of the required input units.
#' 
#' @references
#' Lawson, B.D. & Armitage, O.B., 2008. Weather guide for the Canadian Forest Fire Danger Rating System. Northern Forestry Centre, Edmonton (Canada).
#' 
#' van Wagner, C.E., 1987. Development and structure of the Canadian Forest Fire Weather Index (Forestry Tech. Rep. No. 35). Canadian Forestry Service, Ottawa, Canada.
#' 
#' van Wagner, C.E., Pickett, T.L., 1985. Equations and FORTRAN program for the Canadian forest fire weather index system (Forestry Tech. Rep. No. 33). Canadian Forestry Service, Ottawa, Canada.
#' 
#' Bedia, J., Herrera, S., Camia, A., Moreno, J.M., Gutierrez, J.M., 2014. Forest Fire Danger Projections in the Mediterranean using ENSEMBLES Regional Climate Change Scenarios. Clim Chang 122, 185--199. doi:10.1007/s10584-013-1005-z
#'
#' @seealso \code{\link{fwi1D}} 
#' @author J. Bedia \& M.Iturbide
#' @export
#' @importFrom abind abind asub
#' @importFrom downscaleR subsetGrid
#' @importFrom downscaleR getYearsAsINDEX

fwiGrid <- function(multigrid,
                    mask = NULL,
                    return.all = FALSE, 
                    init.pars = c(85, 6, 15),
                    parallel = FALSE,
                    max.ncores = 16,
                    ncores = NULL) {
      x <-  c(-90, -80, -70, -60, -50, -40, -30, -20, -10, 1, 10, 20, 30, 40, 50, 60, 70, 80, 90)
      latLim <- range(multigrid$xyCoords$y)
      lonLim <- range(multigrid$xyCoords$x)
      latind <- findInterval(latLim, x)[1]:findInterval(latLim, x)[2]
      if (x[latind[length(latind)]] < latLim[2]) latind[3] <- latind[2] + 1
      x <- x[latind]
      lats <- seq(min(x) + 5, max(x) - 5, 10) 
      if (x[length(x)] > latLim[2]) x[length(x)] <- latLim[2]
      if (x[1] < latLim[1]) x[1] <- latLim[1]
      dataset <- attr(multigrid, "dataset")
      years <- unique(getYearsAsINDEX(multigrid))
      dimNames <- downscaleR:::getDim(multigrid)
      latdim <- grep("lat", dimNames)
      a <- list()
      aux.list <- a
      message("[", Sys.time(), "] Calculating FWI..")
      for (i in 1:(length(x) - 1)) {
            latLimchunk <- c(x[i], x[i + 1])
            lat <- lats[i]
            multigrid_chunk <- subsetGrid(multigrid, latLim = latLimchunk)
            aux.list[[i]] <- multigrid_chunk$xyCoords$y
            ## This part avoids repeating latitudes in the extremes as a result of subsetGrid
            if (i > 1) {
                  lat.rep <- intersect(aux.list[[i]], aux.list[[i - 1]])
                  if (length(lat.rep) > 0) {
                        lat.rep.ind <- match(lat.rep, multigrid_chunk$xyCoords$y)
                        multigrid_chunk$Data <- asub(multigrid_chunk$Data, idx = -lat.rep.ind, dims = latdim)
                        attr(multigrid_chunk$Data, "dimensions") <- dimNames
                        multigrid_chunk$xyCoords$y <- multigrid_chunk$xyCoords$y[-lat.rep.ind]
                  }
            }
            if (is.null(mask) & dataset == "WFDEI") {
                  msk <- subsetGrid(multigrid_chunk, var = "tas")
                  msk$Data <- msk$Data[1,,]
                  msk$Data[which(!is.na(msk$Data))] <- 100
                  msk$Data[which(is.na(msk$Data))] <- 0
                  attr(msk$Data, "dimensions") <- c("lat", "lon")
            } else if (!is.null(mask)) {
                  msk <- subsetGrid(mask, latLim = latLimchunk, lonLim = lonLim, outside = TRUE)
            } else {
                  message("The use of a land mask is recommended")
                  msk <- NULL
            }
            o <- lapply(1:length(years), function(x) {
                  multigrid.y <- subsetGrid(multigrid_chunk, years = years[x])
                  suppressMessages(
                        fwi(multigrid = multigrid.y,
                            mask = msk,
                            lat = lat,
                            return.all = return.all, 
                            parallel = parallel, init.pars = init.pars, 
                            max.ncores = max.ncores, ncores = ncores)$Data
                  )
            })
            ind.time <- grep("^time", dimNames)
            a[[i]] <-  unname(do.call("abind", list(o, along = ind.time)))  
      }
      a <- lapply(a, "drop")
      message("[", Sys.time(), "] Done.")
      out <- unname(do.call("abind", list(a, along = latdim - 1)))
      temp <- subsetGrid(multigrid, var = "tas")
      temp$Data <- out
      attr(temp$Data, "dimensions") <- dimNames[-1]
      temp$Variable$varName <- "fwi"
      attr(temp$Variable, "use_dictionary") <- FALSE
      attr(temp$Variable, "description") <- "Fire Weather Index"
      attr(temp$Variable, "units") <-  "adimensional"
      attr(temp$Variable, "longname") <- "Fire Weather Index"
      return(temp)
}