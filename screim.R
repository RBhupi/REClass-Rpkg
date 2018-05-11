#' ---
#' title: "Identification of Stratiform and Convective Radar Echoes using Wavelets"
#' author: "Bhupendra Raut"
#' date: "February 27, 2016"
#' @description This script is based on the script written for Darwin project
#' to identify convection for tracking.
#' https://github.com/RBhupi/Darwin-Rscripts.git
#' check this for earlier version Darwin-Rscripts/cpol_wt_daily.R
#' ---

#'
#' This R script reads netCDF files contaning volume scans from Darwin radar data.
#' A trous wavelet transform (ATWT) is used to separate convective regions from the stratiform regions.
#'
#'
#'Issues:
#'
#' ToDo
#' 1. makes three classes: Convective, Stratiform and Transitional.
#' 2. check verticle structure
#+ echo=FALSE
#==========================================================================================
# Start the clock!
start_time <- proc.time()

#+ echo=TRUE, eval=FALSE, warning=FALSE, error=FALSE, message=FALSE
#' Following R packages are required.
library(ncdf4)    #Read/Write netcdf-4 files
library(RColorBrewer)
library(scales) #for alpha()

#for creating daily output files
library(stringr)
library(plyr)



#' saves clusters of verticle profiles to disk.
#' 
#' This function clusters given verticle profiles using K-Modes algorithm 
#' and saves it in RDS file.
#' @param vprof_sample
#' @param nclust
#' @export 
#' @seealso \url{https://cran.r-project.org/web/packages/klaR/index.html}{klaR} package
#' @seealso \code{\link{saveVertProf}}
clusterProfiles <- function(vprof_sample,  nclust =5) {
    dim_vprof <- dim(vprof_sample)
    #use weighed modes
    prof_clust <- kmodes(vprof_sample, iter.max = 100, modes = nclust, weighted = TRUE) 
    clust<-array(as.numeric(unlist(prof_clust$modes)), dim = c(nclust, dim_vprof[2]))
    
    fout <- paste("./kmodes-clust", nclust, ".RDS", sep="")
    #save file
    saveRDS(prof_clust, file = fout)
}





#'
#'
#'Saves verticle profiles to disk
#'@param fout output file path
#'@param dbz_vol 3d radar reflectvity volume that will be used for classification.
#'@param scale_break scale break obtained from \code{\link{getScaleBreak}}
#'@param vert_range number of verticle levels to be considered for classification
#'@param nsample only few random veticle samples will taken from each image. (Default 10)
#'@export
#'@note Provide only rainy volume scans to this function. Keep \code{nsample} small.

saveVertProf <- function(fout="./profiles.txt", dbz_vol, scale_break, vert_range=1:30, nsample=10){
    wt_class_3d <- getClass(dbz_vol, scale_break)
    vprof_samples <- sampleVertProf(wt_class_3d, nsample)
    write.table(x=vprof_samples, file = fout, append = TRUE, col.names = FALSE, row.names = FALSE)
}


#'sample verticle arrangement of Co and non-Co levels in the txt file for clustering.
#'
#'This function is required to get the verticle structure modes using k-modes algorithm.
#'The function is not used in actual classification run.
#'
#'@param wt_class_3d
#'@param n only few profiles will be saved from this volume for efficiency. default=10. 
#'@seealso \code{link{nearest_kmode}} function.
sampleVertProf <- function(wt_class_3d, n=10, vert_range=1:30){
    dims <- dim(wt_class_3d)
    vstruct <- matrix(data = NA, ncol=length(vert_range), nrow = 0, byrow = TRUE)
    nsample <- 1000
    #select 1000  random locations and sample vert profiles
    x_select <- sample(seq(dims[1]), size = nsample, replace = TRUE)
    y_select <- sample(seq(dims[2]), size = nsample, replace = TRUE)
    counter <-0
    
    #for each column
    for(i in seq(nsample)){
        vert_column <- wt_class_3d[x_select[i], y_select[i], vert_range]
        
        #if column is neither empty nor full of NAs
        if(!all(is.na(vert_column)) && !all(vert_column==0)){
            vstruct<-rbind(vstruct, vert_column)
            counter <- counter+1
            if(counter==n) break()
        }
    }
    return(vstruct)
}



#' Compute scan-by-scan ATWT of radar volume.
#'
#' Converts dBZ to rain rates using standard Z-R relationship.
#' This is to transform the normally distributed dBZ to gamma-like distribution.
#' @param \code{vol_data} 3D array containing radar data. Last dimension should be levels.
#' @param \code{conv_scale} scale break (in pixels) between convective and stratiform scales.
#' @return Sum of wavelets upto \code{conv_scale} for each scan.
#' @export
#' @seealso \code{\linc{getWTSum}}
getClass <- function(vol_data, conv_scale){
    vol_data_T <- dbz2rr(vol_data)
    wt_sum <- getWTSum(vol_data_T, conv_scale)
    #wt_sum <- refineWT_withDBZ(wt_sum, dbz_vol)
    wt_class <- ifelse(wt_sum>10, 2, NA)
    wt_class <- replace(wt_class, is.na(wt_class) & vol_data>-10, 1)
    invisible(wt_class)
}


#' returns sum of WT upto given scale.
getWTSum <- function(vol_data_T, conv_scale) {
    dims <- dim(vol_data_T)
    
    #if data is 2d
    if(length(dims)==2){
        wt <- atwt2d(vol_data_T, max_scale = conv_scale)
        wt_sum<- apply(wt, MARGIN = c(2, 3), FUN = sum)
        invisible(wtco)
    }
    
    #else for volume data
    num_levels <- dims[3]
    wt_sum <- array(data=NA, dim = dims)
    
    for(lev in seq(num_levels)){
        if(max(vol_data_T[, , lev], na.rm=TRUE)<1) next() #this needs reviewing
        wt <- atwt2d(vol_data_T[, , lev], max_scale = conv_scale)
        wt_pos <- replace(wt, wt<0, 0) # remove negative WT
        
        #sum all the WT scales.
        wt_sum[, , lev] <- apply(wt, MARGIN = c(2, 3), FUN = sum)
    }
    
    
    invisible(wt_sum)
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





#' computes rain rate using standard Z-R relationship.
#'
#' @param dbz array, vector or matrix of reflectivity in dBZ
#' @return rr rain rate in \code{mm/hr}
dbz2rr <- function(dbz){
    ZRA=200
    ZRB=1.6
    rr<-((10.0^(dbz/10.0))/ZRA)^(1.0/ZRB)
    return(rr)
}





#'Computes A trous wavelet transform (ATWT)
#'
#'Computes ATWT of the 2d array up to \code{max_scale}.
#'If \code{max_scale} is outside the boundaries, number of scales will be reduced.
#'Data is mirrored at the boundaries.'Negative WT are removed. Not tested for non-square data.
#'@param data2d 2d image as array or matrix.
#'@param max_scale computes wavelets up to \code{max_scale}. Leave blank for maximum possible scales.
#'@return array containing ATWT of input image with added 3rd dimention for scales.
#'@todo Need to break this into smaller functions.
atwt2d <- function (data2d, max_scale=-1){
    #removed missing and negative values
    data2d <- replace(data2d, is.na(data2d)|data2d<1, 0.0)
    
    dims <- dim(data2d)
    ny <- dims[1]
    nx <- dims[2]
    max_possible_scale <- getMaxScale(dims)
    
    
    if(max_scale<0 | max_possible_scale<max_scale)
        max_scale <- max_possible_scale
    
    
    wt <- array(data=0.0, dim = c(max_scale, dims))
    
    sf=c(0.0625, 0.25, 0.375) # function
    
    temp1<-array(data=0.0, dim = dims)
    temp2<-array(data=0.0, dim = dims)
    
    #start Wavelet loop
    for(scale in 1:max_scale){
        x1 <- 2^(scale-1)
        x2 <- 2 * x1
        
        #Row-wise (longitude) smoothing
        for (i in 1:nx){
            
            #find the indices for prev and next points on the line
            prev2 <- abs(i-x2)
            prev1 <- abs(i-x1)
            next1 <- (i+x1)
            next2 <- (i+x2)
            
            #If these indices are outside the image, "mirror" them
            #Sometime this causes issues at higher scales.
            if(next1 > nx) next1 <- 2*nx - next1
            
            if(next2 > nx) next2 <- 2*nx - next2
            
            if(prev1<1 | prev2 <1){
                prev1 <- next1
                prev2 <- next2
            }
            
            
            for (j in 1:(ny)) {
                #print(paste("i=", i,  "j=", j, "scale=", scale, "prev2=", prev2, "prev1=",prev1, "next1=", next1, "next2=", next2))
                left2  <-  data2d[j, prev2]
                left1  <-  data2d[j, prev1]
                right1  <-  data2d[j, next1]
                right2  <-  data2d[j, next2]
                temp1[j, i]  <-  sf[1] * (left2+right2) +
                    sf[2] * (left1 + right1) + sf[3] * data2d[j, i]
            }
        }
        
        
        #column-wise (latitude) smoothing
        for(i in 1:ny){
            
            prev2 <- abs(i-x2)
            prev1 <- abs(i-x1)
            next1 <- (i+x1)
            next2 <- (i+x2)
            
            #If these indices are outside the image use next values
            if(next1 > ny) next1 <- 2*ny - next1
            
            if(next2 > ny) next2 <- 2*ny - next2
            
            if(prev1<1 | prev2 <1){
                prev1 <- next1
                prev2 <- next2
            }
            
            
            
            
            for(j in 1:nx){
                top2  <-  temp1[prev2, j]
                top1  <-  temp1[prev1, j]
                bottom1  <-  temp1[next1, j]
                bottom2  <-  temp1[next2, j]
                temp2[i, j]  <-  sf[1] * (top2+bottom2) +
                    sf[2] * (top1 + bottom1) + sf[3] * temp1[i, j]
            }
        }
        
        wt[scale, , ] <- data2d - temp2
        data2d <- temp2
    }
    invisible(wt)
}


#' Calculate the mximum possible scale of ATWT for given dimensions.
#'
#' @param data_dim output of the \code{dim(data2d)} for given matrix or array.
#' @return integer value of the maximum scale.
getMaxScale<-function(data_dims){
    min_dim <- min(data_dims)
    max_scale <- log(min_dim)/log(2)
    return(floor(max_scale))
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
#' @param wt_class_3d Volume classification obtained from \code{\linc{get_class}}
#' @param vert_clust Clusters saved in RDs files from functions \code{\linc{clusterProfiles}}
#' @param vert_range same vert_range in as \code{\linc{saveVertProf}} and \code{\linc{sampleVertProf}}
#' @return 2d array of pixels labeled with three classes. 1. stratiform, 2. Convection, 3. Mixed (Need correction)
#' @export
#' @seealso \code{\linc{get_class}}
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


#'
matchKModes<-function(vert_prof, kclust){
    if(all(vert_prof==0)){
        return(0)
    }
    
    nclust <- max(kclust$cluster)
    kdist<-c(rep(NA, nclust))
    for(clust in 1:nclust){
        kdist[clust] <- sum(kclust$modes[clust,] != vert_prof)
    }
    
    nearest_match <- which(kdist==min(kdist))
    
    # in case there are two modes exactly same distance apart, we will assign the first one.
    return(nearest_match[1])
}



#====================================================================================#
#                     FUNCTIONS FOR CREATING OUTPUT FILE
#====================================================================================#

#' gives list of files belong to the same day as the first file in list_allFiles.
#'
#' Method: split the first file's name by "_" and extract 'dates' string,
#' then search for the same 'dates' in the list and return them back.
#' @param list_allFiles output of Glob or list.files()
#' @param datestr_revpos the position of date in file name from the reverse direction.
#' @export
get_1dayFiles <- function (list_allFiles, datestr_revpos) {
    firstFile <- list_allFiles[1]
    fname_split <- unlist(strsplit(firstFile, "_"))
    date_str <- fname_split[length(fname_split)-datestr_revpos]
    selectFiles <- str_detect(list_allFiles, date_str)
    flist_select <- list_allFiles[selectFiles]
    print(paste(length(flist_select), "file(s) found on", date_str))
    return(flist_select)
}


#reads time from a single netcdf file. Using with laply.
#I had to write this complecated function  to convert variable time epoch to standard unix epoch.
get_nctime <- function(filename){
    ncfile <- nc_open(filename)

    #the time units is actual time of the scan
    tunit <- ncatt_get(ncfile, varid="time", attname = "units")

    #split it to get the date-time string
    tunits_str <- unlist(str_split(tunit$value, " "))
    # the last string is the date-time
    origin_t <- tunits_str[length(tunits_str)]

    #split date-time string into date and time using "T" as separator
    origin_t <- unlist(str_split(origin_t, "T"))

    #convert time from clock format to seconds
    time_sec <- as.difftime(origin_t[2], format = "%H:%M:%S", units = "s")

    #convert to posix time
    time_posix <- as.POSIXct(as.numeric(time_sec), origin=origin_t[1], tz="UTC")

    #return it as seconds from unix epoch
    return(as.numeric(time_posix))
}




#' @export
get_outFileName <- function(sample_fName){
    fname_split <- unlist(strsplit(sample_fName, "_"))
    date_str <- fname_split[length(fname_split)-3]
    ofName<-paste("CPOL_ATWT_ECHO_CLASSIFICATION_", date_str, ".nc", sep="")
    return(ofName)
}


#'creates output netCDF file usin variables and dims from Steiner classification file.
#'
#'This is done to keep consistency between the files for easier comparison.
#'I dont know if this is still valid.
#'@export
create_outNC<-function(outfPath, time_seconds, sample_file){
    nc_sample <- nc_open(sample_file)
    x_dim <- nc_sample$dim$x
    y_dim <- nc_sample$dim$y

    t_dim <- ncdim_def(name = "time", units = "seconds since 1970-01-01 00:00:00 UTC", calendar = "gregorian",
                       vals= time_seconds, longname = "Time of the scan", unlim = TRUE)


    class_var<-ncvar_def(name="ATWT_ECHO_CLASSIFICATION", units = "", dim=list(y_dim, x_dim, t_dim), missval = -1,
                         longname = "echo classification based on Bhupendra Raut`s method", prec = "integer", shuffle = TRUE,
                         compression = 7, chunksizes = c(x_dim$len, y_dim$len, 1))
    
    #Fastest changing dimention has to be written last in R so it is alway (lon, lat, time)
    lat_dim <- ncvar_def(name = "lat0", units = nc_sample$var$latitude$units, dim = list(y_dim, x_dim),
                         prec = "float", longname = "latitude")

    lon_dim <- ncvar_def(name = "lon0", units = nc_sample$var$longitude$units, dim = list(y_dim, x_dim),
                         prec = "float", longname = "longitude")

    proj_var  <- ncvar_def(name = "grid_mapping_0", units = "", dim=list())

    ofile <- nc_create(filename = outfPath, vars =list(class_var, lat_dim, lon_dim, proj_var))

    ncatt_put(nc = ofile, varid = "ATWT_ECHO_CLASSIFICATION", attname = "coordinates", attval ="lon0 lat0")
    ncatt_put(ofile, varid = "ATWT_ECHO_CLASSIFICATION", attname = "grid_mapping", attval = proj_var$name)

    attnames <- c("grid_mapping_name",
                  "latitude_of_projection_origin",
                  "longitude_of_projection_origin",
                  "CoordinateAxes",
                  "CoordinateAxesTypes",
                  "false_easting",
                  "false_northing",
                  "inverse_flattening",
                  "semi_major_axis")
    attvals <- list("azimuthal_equidistant",
                 -12.2488888888889,
                 131.044444444444,
                 "x y",
                 "GeoX GeoY",
                 0.0,
                 0.0,
                 298.25,
                 6370997.0)
    attvals

    for(i in 1:length(attnames)){
        ncatt_put(nc = ofile, varid = "grid_mapping_0", attname = attnames[i], attval =attvals[[i]])
    }


    #read lat lon from sample file
    lat <- ncvar_get(nc_sample, varid =nc_sample$var$latitude$name, start = c(1, 1, 1, 1), count=c(-1, -1, 1, 1))
    lon <- ncvar_get(nc_sample, varid =nc_sample$var$longitude$name, start = c(1, 1, 1, 1), count=c(-1, -1, 1, 1))



    ncvar_put(ofile, varid = "lat0", vals = lat, start = c(1, 1), count = dim(lat))
    ncvar_put(ofile, varid = "lon0", vals = lon, start = c(1, 1), count = dim(lon))

    nc_close(nc_sample)


    # add flags for echo_type
    flags_val <- c(1, 2, 3, 4, 5)
    flags_meaning <- c("Cg-Cb Cb+Anvil low-St deep-St high-St")
    ncatt_put(nc = ofile, varid = "ATWT_ECHO_CLASSIFICATION", attname = "flag_values", attval =flags_val)
    ncatt_put(nc = ofile, varid = "ATWT_ECHO_CLASSIFICATION", attname = "flag_meanings", attval =flags_meaning)
    ncatt_put(nc = ofile, varid = "ATWT_ECHO_CLASSIFICATION", attname = "valid_range", attval =c(1,5))
    return(ofile)
}



#'compute convective-stratiform scale break for given radar resolution.
res_km <- 2.5
conv_scale_km <- 10
scale_break <- getScaleBreak(res_km, conv_scale_km)
vert_range <-1:30 #levels to be consider for classification.

#get all input file names

setwd("~/projects/screim/data/testdata/")
vclust_object <- readRDS(file="~/projects/screim/R/kmodes-clust5.RDS") #this loads cluster data in an object calle "prof_clust"
indir<-"./"

#read all file names recursively, give correct patterns
flist_all<-list.files(indir, pattern="*_GRIDS_2500m.nc", recursive = T, full.names = T)
print(paste(length(flist_all), "file(s) in the folder."))

outdir <- paste("./atwt_types_test/", sep = "")

dir.create(outdir)


flist_all <- flist_all[1]

while(length(flist_all)>0) {
    flist <- get_1dayFiles(flist_all, 3)
    flist_all <- flist_all[flist_all!=flist] #remaining files are stored for next iteration
    time_seconds<- laply(flist, get_nctime)
    daily_ofname <- paste(outdir, get_outFileName(flist[1]), sep="")
    outNC <- create_outNC(daily_ofname, time_seconds, flist[1])
    var_dim <- outNC$var$ATWT_ECHO_CLASSIFICATION$varsize
    empty_array <- array(data=0, dim=var_dim[1:2])
    
    file_counter <-0
    for(afile in flist){
        file_counter <- file_counter+1
        print(paste("processing file", file_counter, basename(afile)))
        
        nc_file <- nc_open(afile)
        dbz_vol <- ncvar_get(nc_file, varid="corrected_reflectivity", start=c(1, 1, 1, 1), count=c(-1, -1, max(vert_range), -1))
        nc_close(nc_file)
        
        # put zero for insignificant  reflectivity. set this tlimit low after testing is done.
        if(all(is.na(dbz_vol)) | max(dbz_vol, na.rm = TRUE)<10){
            ncvar_put(outNC, varid = "ATWT_ECHO_CLASSIFICATION", vals = empty_array,
                      start = c(1, 1, file_counter), count=c(dim(empty_array), 1))
            next()
        }
        
        wt_class_3d <- getClass(dbz_vol, scale_break)
        if(all(wt_class_3d<1, na.rm = TRUE)) next()
        
        #save_vert_prof(wt_class_3d)
        
        class2d <- class3dTo2d(wt_class_3d, vclust_object)
        ncvar_put(outNC, varid = "ATWT_ECHO_CLASSIFICATION", vals = class2d,
                  start = c(1, 1, file_counter), count=c(dim(class2d), 1))
    }
    nc_close(outNC)
    
    
}

# Check time.
end_time <- proc.time()
print(end_time - start_time)





