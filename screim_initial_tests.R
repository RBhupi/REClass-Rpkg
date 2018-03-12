#' ---
#' title: "Improved Stratiform and Convective Radar Echo Identification Method"
#' author: "Bhupendra Raut"
#' date: "February 27, 2018"
#' ---

#'
#' This R script reads netCDF files 
#' 
#' make three classes: Convective, Stratiform and Transitional
#'
#'Issues: 
#'
#' ToDo
#' 1.
#+ echo=FALSE
#==========================================================================================
# Start the clock!
start_time <- proc.time()

#+ echo=TRUE, eval=FALSE, warning=FALSE, error=FALSE, message=FALSE
#' Following R packages are required.
library(ncdf4)    #Read/Write netcdf-4 files
library(plot3D)
library(RColorBrewer)
library(scales) #for alpha()

#+ echo=FALSE
#----------------------------------------------------------------------functions

#'Function computes WT of the 2d array up to given max_scale.
#'Negative WT are removed. Not tested for non-square data.
#'@param data_matrix
#'@param max_scale
wavelet <- function (data, max_scale){
    #make an output wt array
    dim_data <- dim(data)
    nlat <- dim_data[1]
    nlon <- dim_data[2]
    min_window_size <- 2^(max_scale-1)
    
    if(min(c(nlat, nlon))<min_window_size){
        stop(paste("Error: data array is not large enough to comput ",
                   max_scale, "scales!"))
    }
    
    wt <- array(data=0.0, dim = c(max_scale, nlon, nlat))
    
    sf=c(0.0625, 0.25, 0.375) #weighing function
    
    temp1<-array(dim = dim(data))
    temp2<-array(dim = dim(data))
    
    #start Wavelet loop
    for(scale in 1:max_scale){
        x1 <- 2^(scale-1)
        x2 <- 2 * x1
        
        #Row-wise (longitude) smoothing
        for (i in 1:nlon){
            
            #find the indices for prev and next points on the line
            prev2 <- abs(i-x2)
            prev1 <- abs(i-x1)
            next1 <- (i+x1)
            next2 <- (i+x2)
            
            #If these indices are outside the image, "mirror" them
            #May cause issues at first pixel.
            if(next1 > nlon) next1 <- 2*nlon - next1
            
            if(next2 > nlon) next2 <- 2*nlon - next2
            
            if(prev1<1 | prev2 <1){
                prev1 <- next1
                prev2 <- next2
            }
            
            
            for (j in 1:(nlat)) {
                #print(paste("i=", i,  "j=", j, "scale=", scale, "prev2=", prev2, "prev1=",prev1, "next1=", next1, "next2=", next2))
                left2  <-  data[j, prev2]
                left1  <-  data[j, prev1]
                right1  <-  data[j, next1]
                right2  <-  data[j, next2]
                temp1[j, i]  <-  sf[1] * (left2+right2) +
                    sf[2] * (left1 + right1) + sf[3] * data[j, i]
            }
        }
        
        
        #column-wise (latitude) smoothing
        for(i in 1:nlat){
            
            prev2 <- abs(i-x2)
            prev1 <- abs(i-x1)
            next1 <- (i+x1)
            next2 <- (i+x2)
            
            #If these indices are outside the image use next values
            if(next1 > nlat) next1 <- 2*nlat - next1
            
            if(next2 > nlat) next2 <- 2*nlat - next2
            
            if(prev1<1 | prev2 <1){
                prev1 <- next1
                prev2 <- next2
            }
            
            
            
            
            for(j in 1:nlon){
                top2  <-  temp1[prev2, j]
                top1  <-  temp1[prev1, j]
                bottom1  <-  temp1[next1, j]
                bottom2  <-  temp1[next2, j]
                temp2[i, j]  <-  sf[1] * (top2+bottom2) +
                    sf[2] * (top1 + bottom1) + sf[3] * temp1[i, j]
            }
        }
        
        wt[scale, , ] <- data - temp2
        data <- temp2
    }
    invisible(wt)
}




get_wt_vol<-function(vol_data, conv_scale){
    #removed missing and negative values
    vol_data <- dbz2rr_std(vol_data)
    vol_data <- replace(vol_data, is.na(vol_data)|vol_data<1, 0.0)
    
    
    dims <- dim(vol_data)
    max_scale <- get_max_scale(dims[1:2])
    num_levels <- dims[3]  
    wtco <- array(data=NA, dim = dims)
    wtst <- array(data=NA, dim = dims)
    wttr <- array(data=NA, dim = dims)
    
    #wt <- array(data=NA, dim=c(max_scale, dims))
    
    for(lev in seq(num_levels)){
        wt <- wavelet(vol_data[, , lev], max_scale)
        wtco[, , lev] <- apply(wt[1:conv_scale, , ], 
                               MARGIN = c(2, 3), FUN = sum)
        wttr[, , lev] <- wt[(conv_scale+1), , ]
        
        wtst[, , lev] <- apply(wt[(conv_scale+2):max_scale, , ], 
                               MARGIN = c(2, 3), FUN = sum)
    }
    wtco <- replace(wtco, wtco<=0, NA)
    wttr <- replace(wttr, wttr<=0, NA)
    wtst <- replace(wtst, wtst<=0|vol_data<=0, NA)
    
    #wtco <- replace(wtco, wtco < mean(wtco)+sd(wtco), 0.0)
    #wtst <- replace(wtst, wtst < mean(wtst)+1.5*sd(wtst), 0.0)
    
    wttr <- replace(wttr, !is.na(wtco), NA)
    wtst <- replace(wtst, !is.na(wtco), NA)
    
    invisible(list(wtco=wtco, wtst=wtst, wttr=wttr))
}

get_max_scale<-function(data_dims){
    min_dim <- min(data_dims) #find the smallest dimention of the image
    max_scale <- log(min_dim)/log(2) 
    return(floor(max_scale))   #convert to smaller integer
    
}

dbz2rr_std <- function(dbz){
    ZRA=200
    ZRB=1.6
    rr<-((10.0^(dbz/10.0))/ZRA)^(1.0/ZRB)
    return(rr)
}


#'compute convective-stratiform scale break for given radar resolution.
res_km <- 2.5
conv_scale_km <- 10
scale_break <-log((conv_scale_km/res_km))/log(2)+1
scale_break <- round(scale_break)

setwd("~/projects/screim/data/testdata/")
flist <- Sys.glob("*.nc")

ifile<-flist[1]


nc_file <- nc_open(ifile)
dbz_vol <- ncvar_get(nc_file, varid="corrected_reflectivity")



wt_vol <- get_wt_vol(dbz_vol, scale_break)

cs_class <- array(data=NA, dim = dim(dbz_vol))
cs_class <- ifelse(!is.na(dbz_vol), 1, NA)
cs_class <- ifelse(!is.na(wt_vol$wtco), 2, cs_class)

# read Steiner and Merhala class
fname_steiner <- "./darwin_class/CPOL_STEINER_ECHO_CLASSIFICATION_20161210_level2.nc"
fname_thurai <- "./darwin_class/CPOL_THURAI_ECHO_CLASSIFICATION_20161210_level2.nc"
fname_hydro <- "./darwin_class/CPOL_RADAR_ECHO_CLASSIFICATION_20161210_level2.nc"

nc_steiner <- nc_open(fname_steiner)
nc_thurai <- nc_open(fname_thurai)
nc_hydro <- nc_open(fname_hydro)
steiner <- ncvar_get(nc_steiner, varid = "steiner_echo_classification", 
                     start = c(1, 1, 46), count = c(-1, -1, 1))
thurai <- ncvar_get(nc_thurai, varid = "thurai_echo_classification", 
                    start = c(1, 1, 46), count = c(-1, -1, 1))
hydro <- ncvar_get(nc_hydro, varid = "radar_echo_classification", 
                   start = c(1, 1, 46), count = c(-1, -1, 1))

x<-1:117
y<-1:117
z<-1:41


cols2 <- brewer.pal(name="Set3", n=10)


y_cross=90
x_cross=37
pdf(paste("../../plots/wt_screim_compare_0730_rr_sc", conv_scale_km, "ycross", y_cross, ".pdf", sep=""), width = 10, height = 6, bg = "white")
par(mfrow=c(2, 4))
#image2D(wt_vol$wtco, x=x, y=y, main="convection");grid()
#image2D(wt_vol$wtst, x=x, y=y, main="stratiform");grid()

image2D(dbz_vol[, , 1], x=x, y=y, main="dBZ");grid()
#abline(h=y_cross, lwd=1, col=alpha(colour = "black", 0.5))
abline(h=y_cross, lwd=1, col=alpha(colour = "black", 0.5))

image2D(cs_class[, , 1], breaks=c(0, 1, 2), col=cols2[1:2], x=x, y=y, main="classification");grid()
#abline(h=y_cross, lwd=1, col=alpha(colour = "black", 0.5))
abline(h=y_cross, lwd=1, col=alpha(colour = "black", 0.5))

image2D(steiner,  breaks=c(0, 1, 2), col=cols2[1:2], x=x, y=y, main="Steiner");grid()
abline(h=y_cross, lwd=1, col=alpha(colour = "black", 0.3))

image2D(thurai,  breaks=c(0, 1, 2, 3), col=cols2[1:3], x=x, y=y, main="Thurai");grid()
abline(h=y_cross, lwd=1, col=alpha(colour = "black", 0.3))

image2D(hydro,  breaks=seq(0:10), col=cols2[1:10], x=x, y=y, main="HydroMeteor");grid()
abline(h=y_cross, lwd=1, col=alpha(colour = "black", 0.3))


image2D(dbz_vol[, y_cross, ], x=x, y=z, main="cross section dbz");grid()
image2D(wt_vol$wtco[, y_cross,  ], x=x, y=z, main="cross section conv");grid()
#image2D(wt_vol$wtst[, y_cross, ], x=x, y=z, main="cross section strat");grid()

image2D(cs_class[, y_cross, ], breaks=c(0, 1, 2), col=cols2[1:2], x=x, y=z, main="cross section class");grid()

dev.off()

stop("stopped by user command")

pdf("../../plots/wt_scale_test.pdf", width = 4, height = 6, bg = "white")
par(mfrow=c(3, 2))
for(plt in 1:3){
    image2D(wt_vol_pos[plt, , ], x=x, y=y, main=paste("scale", plt))
    grid(lty=2, lwd=0.5)
}

image2D(wt_sum, x=x, y=y, main="sum of scale 1to2")
grid(lty=2, lwd=0.5)
abline(a = y_cross, b=0, lwd=0.5, col=alpha(colour = "black", 0.5))

image2D(dbz_vol[, , 1], x=x, y=y, main="Radar scan [dbz]")
grid(lty=2, lwd=0.5)
abline(a = y_cross, b=0,lwd=0.5, col=alpha(colour = "black", 0.5))

image2D(dbz_vol[, y_cross, ], x=x, y=z, main="cross section")
grid(lty=2, lwd=0.5)

dev.off()


