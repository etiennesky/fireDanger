#' @title Skill calculation
#' 
#' @description Skill calculation with function veriApply
#' 
#' @param ... same parameters as in function veriApply from package easyVerification
#' @return A grid of the skill scores.
#'  
#' @author M.Iturbide
#' @export
#' @importFrom easyVerification veriApply


fwiSkill <- function(obs, 
                     fcst, 
                     verifun = c(""), 
                     prob, 
                     threshold = NULL, 
                     parallel = FALSE,
                     na.rm = TRUE){
            obs <- interpGrid(obs, getGrid(fcst))
            ensdim <- which(downscaleR:::getDim(fcst) == "member")
            tdim <- which(downscaleR:::getDim(fcst) == "time")
            
            score <- veriApply(verifun = verifun,
                                   fcst = fcst$Data,
                                   obs = obs$Data,
                                   prob = prob,
                                   ensdim = ensdim,
                                   tdim = tdim,
                                   parallel = parallel,
                                   na.rm = na.rm)
            score.grid <- easyVeri2grid(easyVeri.mat = score[[1]],
                                           obs.grid = obs,
                                           verifun = verifun)  
            if(verifun == "EnsRocss"){
            data <- array(dim = c(length(prob)+1, length(obs$xyCoords$y), length(obs$xyCoords$x)))
            for(i in 1:(length(prob)+1)){
                  data[i,,] <- score[[i]]
            }
            score.grid$Data <- data
            }
            attr(score.grid$Data, "dimensions") <- c("cat", "lat", "lon")
            return(score.grid)
}

#End

#' @title Plot grid returned by functionfFwiSkill 
#' 
#' @description Plot grid returned by functionfFwiSkill through function plotClimatology from package downscaleR
#' 
#' @param grid The grid object returned by function fwiSkill 
#'  
#' @author M.Iturbide
#' @export

plotFwiSkill <- function(grid){
      catdim <- which(downscaleR:::getDim(grid) == "cat")
      lc <- dim(grid$Data)[catdim]
      attr(grid$Data, "dimensions") <- c("member", "lat", "lon")
      grid.r <- downscaleR:::redim(grid)
      plotClimatology(climatology(grid.r), backdrop.theme = "coastline", par.strip.text = list(labels =paste("cat", 1:lc)))
}


