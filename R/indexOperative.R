#' @title Index calculation for CFS operative seasonal forecast and csv file generation
#' 
#' @description Implementation of the Canadian Fire Weather Index System 
#' 
#' @param dataset A character string indicating the database to be accessed (partial matching is enabled). 
#' Currently accepted values are "CFSv2_seasonal".
#' @param index Index to be computed, e.g. "fwi".
#' @param destdir Output directory.
#' @param climdir Directory where observational climatologies are stored.
#' @param url ncml directory.
#' @param mask Directory (including file name) of the land-mask grid (*.rda).
#' @param fwi.mask Logical. If TRUE (default), a mask is applied to exclude 
#' land areas which will most likely not burn, like deserts etc.
#' @param dictionary Path to the dictionary for variable homogenization. 
#' is set to TRUE, meaning that the function will internally perform the necessary homogenization steps to return the standard 
#' variables defined in the vocabulary (e.g. variable transformation, deaccumulation...). See details on data homogenization.
#' @param members A vector of integers indicating the members to be loaded. Default to NULL, which loads the default members 
#' (see details for the particularities concerning the CFSv2 dataset). For instance, members=1:5 will retrieve the first five 
#' members of the hindcast. Discontinuous member selection (e.g. members=c(1,5,7)) is allowed. If the requested dataset is not 
#' a forecast or the requested variable is static (e.g. orography) it will be ignored.
#' @param aggr.mem Member aggregation function. A list indicating the name of the
#'  aggregation function in first place, and other optional arguments to be passed to the aggregation function.
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
#' @param lonLim Vector of length = 2, with minimum and maximum longitude coordinates, in decimal degrees, of the bounding box selected.
#'  For single-point queries, a numeric value with the longitude coordinate. If \code{NULL} (default), the whole longitudinal range
#'   is selected (Note that this may lead to a large output object size).
#' @param latLim Same as \code{lonLim}, but for the selection of the latitudinal range.
#' @param operative Path of the R object that contains the fwi of the operative dataset. Default is NULL, so the load and fwi is computed.
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
#' @importFrom stats quantile
#' @importFrom utils write.table
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
#' @importFrom downscaleR getGrid
#' @importFrom downscaleR aggregateGrid
#' @importFrom downscaleR subsetGrid
#' @importFrom downscaleR climatology
#' @importFrom downscaleR interpGrid
#' @import loadeR.ECOMS

indexOperative <- function(dataset = "CFSv2_seasonal_operative",
                         index = "fwi", 
                         climdir = getwd(),
                         destdir = getwd(), 
                         url,
                         mask = NULL,
                         fwi.mask = TRUE,
                         dictionary = TRUE, 
                         members = 1:24, 
                         aggr.mem = list(FUN = NULL),
                         season = NULL, 
                         lonLim = NULL, 
                         latLim = NULL,
                         years = NULL, 
                         leadMonth = 0, 
                         operative = NULL,
                         init.pars = c(85, 6, 15),
                         parallel = FALSE,
                         max.ncores = 16,
                         ncores = NULL){
      if(is.null(operative)){
      operative <- wfwiOP(dataset =dataset,
             index = index, 
             url=url,
             mask = mask,
             dictionary = dictionary, 
             members = members, 
             aggr.mem = aggr.mem,
             season = season, 
             lonLim = lonLim, 
             latLim = latLim,
             years = years, 
             leadMonth = leadMonth, 
             return.all = FALSE, 
             init.pars = init.pars,
             parallel = parallel,
             max.ncores = max.ncores,
             ncores = ncores)     
#       save(operative, file = paste(destdir, "/", dataset,"_", years, "_", season[2], "_", season[length(season)], "_operative.rda"))
      save(operative, file = paste(destdir, "/", "operative_",season[2], "_", season[length(season)], ".rda", sep = ""))
      }else{
            operative <- get(load(operative))
      }
      lati <- operative$xyCoords$y
      loni <- operative$xyCoords$x
      operative <- subsetGrid(operative, lonLim = range(loni), latLim = range(lati), outside = T)
      coords <- expand.grid(operative$xyCoords$x, operative$xyCoords$y)
      lat <- coords[,2]
      lon <- coords[,1]

      cellID <- paste(as.integer(lon*100), "-", as.integer(lat*100), sep = "")
      indexes <-  c("indexmean", "index90", "index30")
      for(i in 1:3){
            index.type <- indexes[i]
            if(index.type == "indexmean"){
                  aggr.clim <- list(FUN = mean, na.rm = T)
            }else if(index.type == "index90"){
                  aggr.clim <- list(FUN = quantile, probs = .9, na.rm = T)
            }else if(index.type == "index30"){
                  aggr.clim <- list(FUN = function(x) sum(x > 30)/length(x))
            }
            dirskill <- paste(climdir, "/", index, gsub("index",replacement = "", x = index.type),"_SKILL_", season[2], "_", season[length(season)], ".rda", sep = "")
            dirClim <- paste(climdir, "/", index, gsub("index",replacement = "", x = index.type),"_WFDEI_", season[2], "_", season[length(season)], ".rda", sep = "")
            if(fwi.mask == TRUE){
                  dirMask <- paste(climdir, "/fwi_mask.rda", sep = "")
                  fwmsk <- get(load(dirMask))
                  fwmsk_sub <- subsetGrid(fwmsk, lonLim = range(lon), latLim = range(lat), outside = T)
                  fwi_Mask <- as.vector(t(fwmsk_sub$Data))
            }
            sk <- get(load(dirskill))
            sk_sub <- subsetGrid(sk, lonLim = range(lon), latLim = range(lat), outside = T)
            ###
            obs <- get(load(dirClim, verbose = T))
            obs_sub <- subsetGrid(obs, lonLim = range(lon), latLim = range(lat), outside = T, season = season[-1])#season = season[-1]
            obs_clim <- downscaleR:::redim(climatology(obs_sub, clim.fun = list(FUN=mean, na.rm = T)), drop = T)
            oper_clim <- downscaleR:::redim(climatology(operative, clim.fun = aggr.clim), drop = T)
            fwiClimatology <- as.vector(t(obs_clim$Data))
            fwiClimatology[which(is.na(fwiClimatology))] <- NA
            fwiPrediction <- as.vector(t(oper_clim$Data))
            ind <- which(is.na(fwiPrediction))
            skill.cat1 <- as.vector(t(sk_sub$Data[1,,]))
            skill.cat2 <- as.vector(t(sk_sub$Data[2,,]))
            skill.cat3 <- as.vector(t(sk_sub$Data[3,,]))
            print("hello")
            if(fwi.mask == TRUE){
                  output <- cbind(cellID, lon, lat, fwiPrediction, fwiClimatology, skill.cat1, skill.cat2, skill.cat3, fwi_Mask)[-ind,]
            }else{
                  output <- cbind(cellID, lon, lat, fwiPrediction, fwiClimatology, skill.cat1, skill.cat2, skill.cat3)[-ind,]
            }
            write.table(output, file = paste(destdir, "/", index, gsub("index",replacement = "", x = index.type), "_", dataset,"_", years, "_", season[2], "_", season[length(season)], "_lm", leadMonth, ".csv", sep=""), sep = " ", na = "", row.names = F, quote = F)
       }
}


#' @title Index calculation for CFS operative seasonal forecast
#' 
#' @description Implementation of the Canadian Fire Weather Index System 
#' 
#' @param dataset A character string indicating the database to be accessed (partial matching is enabled). 
#' Currently accepted values are "CFSv2_seasonal".
#' @param index Index to be computed, e.g. "fwi".
#' @param url ncml directory.
#' @param mask Directory (including file name) of the land-mask grid (*.rda).
#' @param dictionary A logical flag indicating if a dictionary is used for variable homogenization. Default (strongly recommended) 
#' is set to TRUE, meaning that the function will internally perform the necessary homogenization steps to return the standard 
#' variables defined in the vocabulary (e.g. variable transformation, deaccumulation...). See details on data homogenization.
#' @param members A vector of integers indicating the members to be loaded. Default to NULL, which loads the default members 
#' (see details for the particularities concerning the CFSv2 dataset). For instance, members=1:5 will retrieve the first five 
#' members of the hindcast. Discontinuous member selection (e.g. members=c(1,5,7)) is allowed. If the requested dataset is not 
#' a forecast or the requested variable is static (e.g. orography) it will be ignored.
#' @param aggr.mem Member aggregation function. A list indicating the name of the
#'  aggregation function in first place, and other optional arguments to be passed to the aggregation function.
#' @param season An integer vector specifying the desired season (in months, January = 1 ..., December = 12). Options include 
#' one to several (contiguous) months. For full year selections (not possible for all datasets, e.g. seasonal forecasts), 
#' the argument value must be set to season = 1:12. If the requested variable is static (e.g. orography) it will be ignored. 
#' See details on the definition of temporal slices.
#' @param years Optional vector of years to select. Default to all available years. If the requested variable is static (e.g. orography) 
#' it will be ignored. See details on the definition of temporal slices.
#' @param lonLim Vector of length = 2, with minimum and maximum longitude coordinates, in decimal degrees, of the bounding box selected.
#'  For single-point queries, a numeric value with the longitude coordinate. If \code{NULL} (default), the whole longitudinal range
#'   is selected (Note that this may lead to a large output object size).
#' @param latLim Same as \code{lonLim}, but for the selection of the latitudinal range.
#' @param leadMonth Integer value indicating the lead forecast time, relative to the first month of season. Note that leadMonth=1 for 
#' season=1 (January) corresponds to the December initialization. Default to 1 (i.e., 1 lead month forecast). If the dataset is not a 
#' forecast or the requested variable is static (e.g. orography) it will be ignored. A message will be printed on screen in the first 
#' case if its value is different from NULL. See details on initialization times.
#' @param return.all Logical. Should all components of the FWI system be returned?. 
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
#' @export
#' @references
#' Lawson, B.D. & Armitage, O.B., 2008. Weather guide for the Canadian Forest Fire Danger Rating System. Northern Forestry Centre, Edmonton (Canada).
#' 
#' van Wagner, C.E., 1987. Development and structure of the Canadian Forest Fire Weather Index (Forestry Tech. Rep. No. 35). Canadian Forestry Service, Ottawa, Canada.
#' 
#' van Wagner, C.E., Pickett, T.L., 1985. Equations and FORTRAN program for the Canadian forest fire weather index system (Forestry Tech. Rep. No. 33). Canadian Forestry Service, Ottawa, Canada.
#' 

wfwiOP <- function(dataset = "CFSv2_seasonal_operative",
                        index = "fwi", 
                        url,
                        mask = NULL,
                        dictionary = TRUE, 
                        members = 1:24, 
                        aggr.mem = list(FUN = NULL),
                        season = NULL, 
                        years = NULL, 
                        lonLim = NULL, 
                        latLim = NULL,
                        leadMonth = 0, 
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
      
      a <- list()
      coords.x <- list()
      coords.y <- list()
      for(i in 1:(length(x)-1)){
            latLimchunk <- c(x[i], x[i+1])
            lat <- lats[i]
            msk <- subsetGrid(mask, latLim = latLimchunk, lonLim = lonLim, outside = T)
            #SEASONAL
            message("latitude ", i, ",  ", length(x) - i, " remaining")
                  Tm <- NULL
                  REP<- 0
            while(is.null(Tm) & (REP<10)){
                  REP <- REP+1
                  suppressWarnings(
                  Tm <- tryCatch({loadECOMS(dataset,  latLim = latLimchunk, var = "tas", dictionary = dictionary, 
                            members = members, season = season, years = years, leadMonth = leadMonth, time = "DD",
                            aggr.d = "mean", url = url, lonLim = lonLim)}, error = function(err){NULL})
                  )
            }
            H <- NULL
            REP<- 0
            while(is.null(H) & (REP<10)){
                  REP <- REP+1
                  suppressWarnings(
            H <- tryCatch({loadECOMS(dataset,  latLim = latLimchunk, var = "hurs", dictionary = dictionary, 
                           members = members, season = season, years = years, leadMonth = leadMonth, time = "DD",
                           aggr.d = "min", url = url, lonLim = lonLim)}, error = function(err){NULL})
                  )
            }
            r <- NULL
            REP<- 0
            while(is.null(r) & (REP<10)){
                  REP <- REP+1
                  suppressWarnings(
            r <- tryCatch({loadECOMS(dataset,  latLim = latLimchunk, var = "tp", dictionary = dictionary, 
                           members = members, season = season, years = years, leadMonth = leadMonth, time = "DD",
                           aggr.d = "sum", url = url, lonLim = lonLim)}, error = function(err){NULL})
                  )
            }
            W <- NULL
            REP<- 0
            while(is.null(W) & (REP<10)){
                  REP <- REP+1
                  suppressWarnings(
            W <- tryCatch({loadECOMS(dataset,  latLim = latLimchunk, var = "wss", dictionary = dictionary, 
                           members = members, season = season, years = years, leadMonth = leadMonth, time = "DD",
                           aggr.d = "mean", url = url, lonLim = lonLim)}, error = function(err){NULL})
                  )
            }
            dimNames <- attr(W$Data, "dimensions")
            initatr <- W$InitializationDates
            mematr <- W$Members
            datatr <- attr(W, "dataset")
            dates <- W$Dates
#             Tm <- downscaleR:::redim(Tm)
#             H <- downscaleR:::redim(H)
#             r <- downscaleR:::redim(r)
#             W <- downscaleR:::redim(W)
#             
            W$Data <- W$Data*3.6
            attr(W$Data, "dimensions") <- attr(r$Data, "dimensions") 
            multigrid <- makeMultiGrid(Tm, H, r, W)
            
            
            xx <- downscaleR:::getShape(Tm)[["lat"]]-1
            Tm <- H <- r <- W <- NULL
            if(index == "fwi"){
            b <- fwi(multigrid, mask = msk, lat = lat, return.all = return.all, 
                     parallel = parallel, init.pars = init.pars, 
                     max.ncores = max.ncores, ncores = ncores)$Data[,,1:xx,,drop=FALSE]
            }
            gc()
           
            months <- as.integer(substr(multigrid$Dates[[1]]$start, start = 6, stop = 7))
            coords.x <- multigrid$xyCoords$x
            coords.y[[i]] <- multigrid$xyCoords$y[-(length(multigrid$xyCoords$y))]
            multigrid <- NULL
            month.ind <- which(months == months[1])
            a[[i]] <- b[,-month.ind,,]
            dates$start <- dates$start[-month.ind]
            dates$end <- dates$end[-month.ind]
      }
      c.x <- coords.x
      c.y <- unname(do.call("abind", list(coords.y)))
      xy <- list("x" = c.x, "y" = c.y)
      out <- unname(do.call("abind", list(a, along = 3)))
outgrid <- list()
outgrid$Variable <- list("varName" = index, "level" = NULL)
outgrid$Data <- out
attr(outgrid$Data, "dimensions") <- dimNames
outgrid$xyCoords <- xy
outgrid$Dates <- dates
outgrid$InitializationDates <- initatr
outgrid$Members <- mematr
attr(outgrid, "dataset") <- datatr
outgrid <- aggregateGrid(outgrid, aggr.mem = aggr.mem)

#     write.csv(out, file = paste(destdir, "/wfi_", dataset, "_", season[2], "_", season[length(season)], ".csv", sep=""))
#     write.csv(xy, file = paste(destdir, "/coordinates_", dataset,"_", season[2], "_", season[length(season)], ".csv", sep=""))
return(outgrid)
}



