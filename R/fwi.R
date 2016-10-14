#' @title Fire Weather Index applied to multigrids
#' 
#' @description Implementation of the Canadian Fire Weather Index System for multigrids
#' 
#' @param multigrid containing Tm (temperature records in deg. Celsius); H (relative humidity records in \%);
#' r (last 24-h accumulated precipitation in mm); W (wind velocity records in Km/h). See the Warning section.
#' @param mask Optional. Binary (0 an 1) Grid.
#' @param latLim Same as \code{lonLim} argument, but for latitude.
#' @param lonLim Vector of length = 2, with minimum and maximum longitude coordinates, in decimal degrees,
#' @param lat Optional. Latitude of the records (in decimal degrees). Default to 46,
#' applying the default parameters of the original FWI System, developed in Canada. See details.
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
#' @section Warning:
#' 
#' The variables composing the input multigrid need to have standard names, as defined by the dictionary
#'  (these names are stored in the \code{multigrid$Variable$varName} component).
#' These are: \code{"tas"} for temperature, \code{"tp"} for precipitation, \code{"wss} for windspeed. In the case of relative humidity,
#' either \code{"hurs"} or \code{"hursmin"} are accepted, the latter in case of FWI calculations according to the \dquote{proxy} version
#' described in Bedia \emph{et al} 2014.
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
#' @author J. Bedia & M. Iturbide, partially based on the original FORTRAN code by van Wagner and Pickett (1985)
#' @export
#' @importFrom abind abind
#' @importFrom transformeR redim
#' @importFrom transformeR getShape


fwi <- function(multigrid, mask = NULL, lonLim = NULL, latLim = NULL, lat = 46, return.all = FALSE, init.pars = c(85, 6, 15),
                parallel = FALSE,
                max.ncores = 16,
                ncores = NULL){
      months <- as.integer(substr(multigrid$Dates[[1]]$start, start = 6, stop = 7))
      parallel.pars <- parallelCheck(parallel, max.ncores, ncores)
      if (parallel.pars$hasparallel) {
            apply_fun <- function(...) {
                  parallel::parLapply(cl = parallel.pars$cl, ...)
            }
            on.exit(parallel::stopCluster(parallel.pars$cl))
      } else {
            apply_fun <- lapply      
      }
      #       if (!is.null(lonLim)){
      #             multigrid <- subsetGrid(multigrid, lonLim = lonLim)
      #       }
      #       if (!is.null(latLim)){
      #             multigrid <- subsetGrid(multigrid, lonLim = latLim)
      #       }
      cords <- getCoordinates(multigrid)
      if (!is.null(mask)) {
            mask1 <- transformeR:::redim(mask,
                                        member = FALSE,
                                        runtime = FALSE,
                                        drop = FALSE)
      }
      varnames <- multigrid$Variable$varName
      Tm1 <- subsetGrid(multigrid, var = grep("tas", varnames, value = TRUE))
      Tm1 <- transformeR:::redim(Tm1, drop = FALSE)
      H1  <- subsetGrid(multigrid, var = grep("hurs", varnames, value = TRUE))
      H1 <- transformeR:::redim(H1, drop = FALSE)
      r1  <- subsetGrid(multigrid, var = "tp")
      r1 <- transformeR:::redim(r1, drop = FALSE)
      W1  <- subsetGrid(multigrid, var = "wss")
      W1 <- transformeR:::redim(W1, drop = FALSE)
      fwigrid <- W1
      n.mem <- transformeR:::getShape(W1, "member")
      message("[", Sys.time(), "] Calculating FWI...")
      a <- apply_fun(1:n.mem, function(x) {
            Tm2 <- array3Dto2Dmat(subsetGrid(Tm1, members = x)$Data)
            H2 <- array3Dto2Dmat(subsetGrid(H1, members = x)$Data)
            r2 <- array3Dto2Dmat(subsetGrid(r1, members = x)$Data)
            W2 <- array3Dto2Dmat(subsetGrid(W1, members = x)$Data)
            if (!is.null(mask)) {
                  mskmsk <- array3Dto2Dmat(mask1$Data)[1,]
                  ind <- which(mskmsk > 0)
            } else {
                  ind <- which(!is.na(Tm2))
            }
            b <- array(dim = dim(Tm2))
            if (length(ind) != 0) {
                  for (i in 1:length(ind)) {
                        z <- tryCatch({fwi1D(months,
                                             Tm = Tm2[,ind[i]],
                                             H = H2[,ind[i]],
                                             r = r2[,ind[i]],
                                             W = W2[,ind[i]],
                                             lat = lat,
                                             return.all = return.all,
                                             init.pars = init.pars)},
                                      error = function(err){rep(NA, nrow(b))})
                        if (length(z) < length(b[,1])) z <- rep(NA, length(b[,1]))
                        b[,ind[i]] <- z
                  }
            }
            c <- mat2Dto3Darray(mat2D = b, x = cords$x, y = cords$y)
            return(c)
      })
      fwidat <- do.call("abind", list(a, along = 0))
      message("[", Sys.time(), "] Done.")
      dimNames <- attr(fwigrid$Data, "dimensions")
      fwigrid$Data <- unname(fwidat)
      attr(fwigrid$Data, "dimensions") <- dimNames
      fwigrid$Variable <- list()
      fwigrid$Variable$varName <- "fwi" 
      return(fwigrid)
}

