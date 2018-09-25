#' title: "Test script for SCREIM"
#' author: "Bhupendra Raut"
#' date: "September 24, 2018"
#' @description This script is for testing SCREIM package with Darwin radar data.
#' ---

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
library(plot3D)
library(RColorBrewer)

#'compute convective-stratiform scale break for given radar resolution.
res_km <- 2.5
conv_scale_km <- 20
scale_break <- getScaleBreak(res_km, conv_scale_km)
vert_range <-1:30 #levels to be consider for classification.

setwd("~/projects/screim/data/testdata/")
#vclust_object <- readRDS(file="~/projects/screim/R/kmodes-clust5.RDS") #this loads cluster data in an object calle "prof_clust"


dbz_file <- nc_open("./CPOL_20161210_1320_GRIDS_2500m.nc")
dbz_vol <- ncvar_get(dbz_file, varid="corrected_reflectivity", 
                     start=c(1, 1, 1, 1), count=c(-1, -1, max(vert_range), -1))

dbz_vol <- replace(dbz_vol, dbz_vol<10, NA)
nc_close(dbz_file)

wt_class_3d <- getWTClass(dbz_vol, scale_break)
dbz_vol_t <- dbz2rr(dbz_vol)
wt_sum_dbz <- getWTSum(dbz_vol_t, scale_break)
#wt_sum_dbz <- replace(wt_sum_dbz, wt_sum_dbz<10, NA)

hyd_file <- nc_open("./darwin_class/CPOL_RADAR_ECHO_CLASSIFICATION_20161210_level2.nc")
stn_file <- nc_open("./darwin_class/CPOL_STEINER_ECHO_CLASSIFICATION_20161210_level2.nc")

t_ind<-81 #46 #81
steiner <- ncvar_get(stn_file, varid = "steiner_echo_classification", 
                     start = c(1, 1, t_ind), count = c(-1, -1, 1))
steiner <- replace(steiner, steiner<1, NA)
hydro <- ncvar_get(hyd_file, varid = "radar_echo_classification", 
                   start = c(1, 1, t_ind), count = c(-1, -1, 1))





cols <- brewer.pal(name="Set3", n=11)
cols1 <- brewer.pal(name="Accent", 6)
cols_dbz <-brewer.pal(name="YlGnBu", 9)

cols1 <-rev(cols1[c(T, T, T, F, T, T)])
col_2class <- cols[5:6]
col_3class <- cols[5:7]
col_10class <- cols <- brewer.pal(name="Set1", n=9)





radar_range <- (-58:58)*2.5

pdf("../../plots/wt_screim_2d_comparison.pdf", width = 5.5, height = 10, bg = "white")
layout(mat=matrix(data=c(1:4, 5, 5, 6, 6), nrow = 4, ncol = 2, byrow = T))
image2D(dbz_vol[, , 6], x=radar_range, y=radar_range, xlab="[Km]", ylab="[Km]",
        main="a) Reflectivity [dBZ], 2Km CAPPI");grid(lty=2)
abline(h=-122.5, col="magenta", lwd=1, lty=3) #cross section at 75 Km

image2D(wt_class_3d[, , 6], x=radar_range, y=radar_range, xlab="[Km]", ylab="[Km]", 
        main="b) WT Class", breaks=c(0, 1, 2), col=col_2class, colkey = F);grid(lty=2)
legend("topleft", fill = rev(col_2class), legend = rev(c("2. Non-convective", "1. Convective")), 
       bty="n", border = T)


image2D(hydro, x=radar_range, y=radar_range, xlab="[Km]", ylab="[Km]", colkey = F,
        main="c) Hydrometeor Class", breaks=seq(0:5), col=col_10class[2:6]);grid(lty=2)
legend("topleft", fill = col_10class[2:6], legend = c("1: Drizzle", "2: Rain", "3: Ice Crystals", "4: Aggregates",
                                                          "5: Wet Snow"),  bty="n", border = T)


image2D(steiner, x=radar_range, y=radar_range, xlab="[Km]", ylab="[Km]", 
        main="d) Steiner Class",breaks=c(0, 1, 2), col=col_2class, colkey = F);grid(lty=2)
legend("topleft", fill = rev(col_2class), legend = rev(c("2. Non-convective", "1. Convective")), 
       bty="n", border = T)

vrange<-(1:30)*0.5
image2D(dbz_vol[,10 , ], x=radar_range, y=vrange, xlab="[Km]", ylab="[Km]", 
        main="e) Vertical Cross-section of Reflectivity");grid(lty=2)
contour2D(wt_class_3d[,10 , ], x=radar_range, y=vrange, colkey = F, col = "black", add=T, 
          lwd=.1, resfac = 2, clab = F)

image2D(dbz_vol[,10 , ], x=radar_range, y=vrange, xlab="[Km]", ylab="[Km]", 
      main="e) Vertical Cross-section of Reflectivity");grid(lty=2)
contour2D(wt_sum_dbz[,10 , ], x=radar_range, y=vrange, colkey = F, col = "black", add=T, 
                lwd=1, breaks=seq(10, 30, 5))

#

dev.off()



