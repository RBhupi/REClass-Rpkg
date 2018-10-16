

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
    selectFiles <- stringr::str_detect(list_allFiles, date_str)
    flist_select <- list_allFiles[selectFiles]
    print(paste(length(flist_select), "file(s) found on", date_str))
    return(flist_select)
}



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
    tunits_str <- unlist(stringr::str_split(tunit$value, " "))
    # the last string is the date-time
    origin_t <- tunits_str[length(tunits_str)]
    
    #split date-time string into date and time using "T" as separator
    origin_t <- unlist(stringr::str_split(origin_t, "T"))
    
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
create_outNC_class<-function(outfPath, time_seconds, sample_file){
    nc_sample <- nc_open(sample_file)
    x_dim <- nc_sample$dim$x
    y_dim <- nc_sample$dim$y
    
    t_dim <- ncdim_def(name = "time", units = "seconds since 1970-01-01 00:00:00 UTC", calendar = "gregorian",
                       vals= time_seconds, longname = "Time of the scan", unlim = TRUE)
    
    
    class_var<-ncvar_def(name="ATWT_ECHO_CLASSIFICATION", units = "", dim=list(y_dim, x_dim, t_dim), missval = -1,
                         longname = "echo classification based on Bhupendra Raut`s method", prec = "integer", shuffle = TRUE,
                         compression = 7, chunksizes = c(x_dim$len, y_dim$len, 1))
    
    #Fastest changing dimention has to be written last in R so it is alway (lon, lat, time)
    if(!is.null(nc_sample$var$latitude$units)){
        lat_units <- nc_sample$var$latitude$units
        lon_units <- nc_sample$var$longitude$units
        lat_name <- nc_sample$var$latitude$name
        lon_name <- nc_sample$var$longitude$name
    } else if(!is.null(nc_sample$var$lat0$units))
    {
        lat_units <- nc_sample$var$lat0$units
        lon_units <- nc_sample$var$lon0$units
        lat_name <- nc_sample$var$lat0$name
        lon_name <- nc_sample$var$lon0$name
        
    }
    
    else
        stop("sample file has different naming conventions than expected.")
    
    lat_dim <- ncvar_def(name = "lat0", units = lat_units, dim = list(y_dim, x_dim),
                         prec = "float", longname = "latitude")
    
    lon_dim <- ncvar_def(name = "lon0", units = lon_units, dim = list(y_dim, x_dim),
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
    lat <- ncvar_get(nc_sample, varid =lat_name, start = c(1, 1, 1, 1), count=c(-1, -1, 1, 1))
    lon <- ncvar_get(nc_sample, varid =lon_name, start = c(1, 1, 1, 1), count=c(-1, -1, 1, 1))
    
    
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
