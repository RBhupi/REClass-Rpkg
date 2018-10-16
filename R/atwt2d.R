
#'Computes A trous wavelet transform (ATWT)
#'
#'Computes ATWT of the 2d array up to \code{max_scale}.
#'If \code{max_scale} is outside the boundaries, number of scales will be reduced.
#'Data is mirrored at the boundaries.'Negative WT are removed. Not tested for non-square data.
#'@param data2d 2d image as array or matrix.
#'@param max_scale computes wavelets up to \code{max_scale}. Leave blank for maximum possible scales.
#'@return array containing ATWT of input image with added 3rd dimention for scales.
#'@note Need to break this into smaller functions.
#'@author Bhupendra Raut and Dileep M. Puranik
#'@seealso Press et al. (1992) Numerical Recipes in C.
atwt2d <- function (data2d, max_scale=-1){
    #removed missing and negative values
    data2d <- replace(data2d, is.na(data2d), 0.0)
    
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




