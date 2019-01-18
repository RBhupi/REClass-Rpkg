#' saves clusters of verticle profiles to disk.
#' 
#' @description This function clusters given verticle profiles using K-Modes algorithm 
#' and saves it in RDS file.
#' @param vprof_sample saved samples from the function \code{\link{saveVertProf}}
#' @param nclust Desired number of clusters. (Default 5)
#' @export 
#' @seealso \url{https://cran.r-project.org/web/packages/klaR/index.html}{klaR} package
#' @seealso \code{\link{saveVertProf}}
clusterProfiles <- function(vprof_sample,  nclust =5) {
    dim_vprof <- dim(vprof_sample)
    #use weighed modes
    prof_clust <- klaR::kmodes(vprof_sample, iter.max = 100, modes = nclust, weighted = TRUE) 
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
    wt_class_3d <- getWTClass(dbz_vol, scale_break)
    vprof_samples <- sampleVertProf(wt_class_3d, nsample)
    write.table(x=vprof_samples, file = fout, append = TRUE, col.names = FALSE, row.names = FALSE)
}


#'sample verticle arrangement of Co and non-Co levels in the txt file for clustering.
#'
#'This function is required to get the verticle structure modes using k-modes algorithm.
#'The function is not used in actual classification run.
#'
#'@param wt_class_3d Volume of WT classification from \code{\link{getWTClass}}
#'@param n only few profiles will be saved from this volume for efficiency. default=10. 
#'@seealso \code{link{nearest_kmode}} function.
sampleVertProf <- function(wt_class_3d, n=10, vert_range=1:30){
    dims <- dim(wt_class_3d)
    vstruct <- matrix(data = NA, ncol=length(vert_range), nrow = 0, byrow = TRUE)
    nsample <- n
    #select 1000  random locations and sample vert profiles
    x_select <- sample(seq(dims[1]), size = nsample, replace = TRUE)
    y_select <- sample(seq(dims[2]), size = nsample, replace = TRUE)
    counter <-0
    
    #for each column
    for(i in seq(nsample)){
        vert_column <- wt_class_3d[x_select[i], y_select[i], vert_range]
        
        #if column is neither empty nor full of NAs
        if(sum(!is.na(vert_column), na.rm = TRUE)>3 & !all(vert_column==0)){
            vstruct<-rbind(vstruct, vert_column)
            counter <- counter+1
            if(counter==n) break()
        }
    }
    return(vstruct)
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
