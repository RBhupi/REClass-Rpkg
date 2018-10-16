
#' Compute scan-by-scan ATWT of radar volume.
#'
#' Converts dBZ to rain rates using standard Z-R relationship.
#' This is to transform the normally distributed dBZ to gamma-like distribution.
#' @param \code{vol_data} 3D array containing radar data. Last dimension should be levels.
#' @param \code{conv_scale} scale break (in pixels) between convective and stratiform scales.
#' @return Sum of wavelets upto \code{conv_scale} for each scan.
#' @export
#' @seealso \code{\link{getWTSum}}
getWTClass <- function(dbz_data, conv_scale){
    dbz_data_t <- dbz2rr(dbz_data) #transform the dbz data
    wt_sum <- getWTSum(dbz_data_t, conv_scale)
    wt_class <- labelClasses(wt_sum, dbz_data)
    invisible(wt_class)
}


#' computes rain rate using Z-R relationship.
#'
#'Uses standard values for ZRA=200 and ZRB=1.6.
#' @param dbz array, vector or matrix of reflectivity in dBZ
#' @return rr rain rate in \code{mm/hr}
dbz2rr <- function(dbz, ZRA=200, ZRB=1.6){
    rr<-((10.0^(dbz/10.0))/ZRA)^(1.0/ZRB)
    return(rr)
}



#' returns sum of WT upto given scale.
#' @export
getWTSum <- function(vol_data, conv_scale) {
    dims <- dim(vol_data)
    
    #if data is 2d
    if(length(dims)==2){
        wt <- atwt2d(vol_data, max_scale = conv_scale)
        wt_sum<- apply(wt, MARGIN = c(2, 3), FUN = sum)
    } else { 			#else for volume data
    	num_levels <- dims[3]
    	wt_sum <- array(data=NA, dim = dims)
    
    	for(lev in seq(num_levels)){
        	if(max(vol_data[, , lev], na.rm=TRUE)<1) next() #this needs reviewing
        	wt <- atwt2d(vol_data[, , lev], max_scale = conv_scale)
        
        	#sum all the WT scales.
        	wt_sum[, , lev] <- apply(wt, MARGIN = c(2, 3), FUN = sum)
    	}
    }
    
    #negative WT do not corresponds to convection in radar
    wt_sum <- replace(wt_sum, wt_sum<0, 0) #check calling function
    invisible(wt_sum)
}



#' Lables 1. stratiform, 2. intense/heavy convective  and 3. moderate+transitional convective regions 
#' using given thersholds.
labelClasses <- function(wt_sum, vol_data) {
    conv_wt_threshold <- 5 #WT value more than this is strong convection
    tran_wt_threshold <- 2 #WT value for moderate convection
    min_dbz_threshold <- 10
    con_dbz_threshold <- 25 # pixel below this value are not convective
    
    wt_class <- ifelse(wt_sum>=conv_wt_threshold & vol_data>=con_dbz_threshold, 2, NA)
    wt_class <- ifelse(wt_sum<conv_wt_threshold & wt_sum>=tran_wt_threshold & vol_data>=con_dbz_threshold, 
                       3, wt_class)
    wt_class <- replace(wt_class, is.na(wt_class) & vol_data>=min_dbz_threshold, 1)
    invisible(wt_class)
}


#' Remove tiny fluctuations that may not be of interest.
#'
#' This may not be needed for clean dataset.
#' Use when WT has too much noise with trial and error approach.
#' @param wt_scan 2d wt image at a scale or sum of several scales.
#' @param times_sd Default value=1. pixels with WT value < mean + times_sd * SD are removed.
#' @return wt_scan with small values removed.
cleanWT <-function(wt_sum, dbz_vol){
    wt_sum<- replace(wt_sum, wt_sum < 10 , 0.0)
    return(wt_sum)
}




#' compute scale break for convection and stratiform regions.
#' 
#' WT will be computed upto this scale and features will be designated as convection.
#' @param res_km resolution of the image.
#' @param conv_scale_km expected size of spatial variations due to convection.
#' @return dyadic integer scale break in pixels.  
#' @export
getScaleBreak<- function(res_km, conv_scale_km){
    scale_break <-log((conv_scale_km/res_km))/log(2)+1
    return(round(scale_break))
}



#' 2D projection of 3D convective-stratiform classes.
#'
#' Checks verticle profile of the classification and finds continuous
#' regions of similar classification and assigns one dominent class.
#' If both classes has comparable presence, mixed class is assigned.
#' @param wt_class_3d Volume classification obtained from \code{\link{get_class}}
#' @param vert_clust Clusters saved in RDs files from functions \code{\link{clusterProfiles}}
#' @param vert_range same vert_range in as \code{\link{saveVertProf}} and \code{\link{sampleVertProf}}
#' @return 2d array of pixels labeled with three classes. 1. stratiform, 2. Convection, 3. Mixed (Need correction)
#' @export
#' @seealso \code{\link{get_class}}
class3dTo2d <- function(wt_class_3d, vert_clust, vert_range=1:30){
    wt_class_3d <- replace(wt_class_3d, is.na(wt_class_3d), 0)
    
    dims <- dim(wt_class_3d)
    
    class2d <- array(data=NA, dim = dims[1:2])
    
    #for each column
    for(x in seq(dims[1])){
        for(y in seq(dims[2])){
            vert_column <- wt_class_3d[x, y, vert_range]
            class2d[x, y] <- matchKModes(vert_column, vert_clust)
        }
    }
    
    invisible(class2d)
}


