#' @title World Fire Weather Index for CFS hindcast
#' 
#' @description Implementation of the Canadian Fire Weather Index System 
#' 
#' @param dataset A character string indicating the database to be accessed (partial matching is enabled). 
#' Currently accepted values are "CFSv2_seasonal".
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
#' @importFrom abind abind
#' @importFrom downscaleR makeMultiGrid
#' @import loadeR.ECOMS

wfwi_cfs <- function(dataset, mask = NULL, dictionary = TRUE, 
                 members = NULL, season = NULL,
                 years = NULL, leadMonth = 0, return.all = FALSE, init.pars = c(85, 6, 15),
                 parallel = FALSE,
                 max.ncores = 16,
                 ncores = NULL){
#       coords <- loadECOMS(dataset, var = "tas", dictionary = TRUE, 
#                           members = 1, season = 1, years = 1991, leadMonth = 1, time = "DD",
#                           aggr.d = "mean", aggr.m = "mean")$xyCoords
#       
#       ind <- rep(c(1,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA), length(coords$x)/24)
#       x <- na.omit(coords$x*ind)
#       p <- which(x==0)
#       if(length(p)!=0){
#             x <- x[-(which(x==0))]
#       }
      x <-  c(-61.5, -50, -40, -30, -20, -10, 1, 10, 20, 30, 40, 50, 60, 70, 76)
      lats <- c(-55, -45, -35, -25, -15, 0, 5, 15, 25, 35, 45, 55, 65, 72) # L&A p71
      
      a <- list()
      for(i in 1:(length(x)-1)){
            latLim <- c(x[i], x[i+1])
            lat <- lats[i]
#             if(is.na(lonLim[2])) lonLim[2] <- 180
            msk <- subsetGrid(mask, latLim = latLim, outside = T)
            message("latitude ", i, ",  ", length(x) - i, " remaining")
            Tm <- loadECOMS(dataset,  latLim = latLim, var = "tas", dictionary = dictionary, 
                            members = members, season = season, years = years, leadMonth = leadMonth, time = "DD",
                            aggr.d = "mean")
            H <- loadECOMS(dataset,  latLim = latLim, var = "hurs", dictionary = dictionary, 
                           members = members, season = season, years = years, leadMonth = leadMonth, time = "DD",
                           aggr.d = "min")
            r <- loadECOMS(dataset,  latLim = latLim, var = "tp", dictionary = dictionary, 
                           members = members, season = season, years = years, leadMonth = leadMonth, time = "DD",
                           aggr.d = "sum")
            W <- loadECOMS(dataset,  latLim = latLim, var = "wss", dictionary = dictionary, 
                           members = members, season = season, years = years, leadMonth = leadMonth, time = "DD",
                           aggr.d = "mean")
            W$Data <- apply(W$Data, MARGIN = 1:4, function(x) x*3.6)
            attr(W$Data, "dimensions") <- attr(r$Data, "dimensions") 
            multigrid <- makeMultiGrid(Tm, H, r, W)
            xx <- dim(Tm$Data)[3]-1
            Tm <- H <- r <- W <- NULL
            b <- fwi(multigrid, mask = msk, lat = lat, return.all = return.all, 
                     parallel = parallel, init.pars = init.pars, 
                     max.ncores = max.cores, ncores = ncores)$Data[,,1:xx,,drop=FALSE]
            a[[i]] <- b
            save(list=c("b"), file = paste("temp/world_fwi_cfs_", i,"_", years, ".rda", sep=""), compress = "xz")
            b <- NULL
            gc()
      }
      out <- unname(do.call("abind", list(a, along = 3)))
      return(out)
}




