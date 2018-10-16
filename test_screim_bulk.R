#' title: "A script for classification of radar echoes using SCREIM package."
#' author: "Bhupendra Raut"
#' date: "October 5, 2018"

#' ??This R script reads netCDF files contaning volume scans from Darwin radar data.
#' ??A trous wavelet transform (ATWT) is used to separate convective regions from the stratiform regions.
#'
#'
#'Issues:
#'
#' ToDo
#' 1. 
#+ echo=FALSE
#==========================================================================================
# Start the clock!
start_time <- proc.time()
library(screim)
library(stringr)
library(plyr)


#' reads time from a single netcdf file. Using with laply.
#' 
#' Titan/CF radial generated files do not have constatnt time epoch. The time units changes for every file.
#' Therefore, I had to write this complicated function  to convert variable time epoch to standard unix epoch.
#' @note This function has sevarl assumptions about the way time units are written. 
#' It may not work with files other than CF-radial.
#' @param filename file name 
#' @return time in POSIX standard with "1970-01-01" as origin.
#' @export
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
    nc_close(ncfile)
    #return it as seconds from unix epoch
    return(as.numeric(time_posix))
}


bulk_proc_func <- function(dir_path, out_dir) {
    flist <- list.files(dir_path, pattern = ".nc", full.names = TRUE)
    all_file_times <- unlist(lapply(flist, FUN = get_nctime))
    
    out_fname <- stringr::str_replace(basename(flist[1]), ".nc", "_class.nc")
    outnc <- create_outNC_class(outfPath = file.path(out_dir, out_fname), all_file_times, sample_file = flist[1])
    var_dim <- outnc$var$ATWT_ECHO_CLASSIFICATION$varsize
    empty_array <- array(data=0, dim=var_dim[1:2])
    
    
    
    file_counter <-0
    for(afile in flist){
        file_counter <- file_counter+1
        print(paste("processing file", file_counter, basename(afile)))
        
        nc_file <- nc_open(afile)
        dbz_vol <- ncvar_get(nc_file, varid="corrected_reflectivity", 
                             start=c(1, 1, 6, 1), count=c(-1, -1, 1, -1))
        nc_close(nc_file)
        
        
        # put zero for insignificant  reflectivity. set this tlimit low after testing is done.
        if(all(is.na(dbz_vol)) | max(dbz_vol, na.rm = TRUE)<10){
            ncvar_put(outnc, varid = "ATWT_ECHO_CLASSIFICATION", vals = empty_array,
                      start = c(1, 1, file_counter), count=c(dim(empty_array), 1))
            next()
        }
        
        wt_class <- getWTClass(dbz_vol, conv_scale = scale_break)
        
        if(all(wt_class<1, na.rm = TRUE)) next()
        
        #save_vert_prof(wt_class_3d)
        
        ncvar_put(outnc, varid = "ATWT_ECHO_CLASSIFICATION", vals = wt_class,
                  start = c(1, 1, file_counter), count=c(dim(wt_class), 1))
    }
    nc_close(outnc)
    
}






#'compute convective-stratiform scale break for given radar resolution.
res_km <- 2.5
conv_scale_km <- 20
scale_break <- getScaleBreak(res_km, conv_scale_km)
vert_range <-1:30 #levels to be consider for classification.

setwd("~/projects/screim/data/Darwin/201701/")
dir_list <- list.dirs(getwd(), full.names = TRUE)
dir_list<-dir_list[-1]


out_dir <-file.path(getwd(), "outClass")
dir.create(out_dir)

l_ply(dir_list, .fun = bulk_proc_func, out_dir)
