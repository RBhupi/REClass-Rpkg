#' This cript reads txt file containing verticle profiles of convective and non-convective 
#' classes and uses a random sample from it to cluster the profiles using kmodes. 
#' These clusters are saved as an Robject for further classification.
#'  
library(RColorBrewer)

setwd("~/projects/screim/data/Darwin/")

# read prifile data from the tex file written by save_vert_struct() in screim.R
vprof <- read.table("./profiles.txt", sep = " ", header = FALSE)
dims <- dim(vprof)

#Use limited sample for clustering.
x <- sample(1:dims[1], size = 10000, replace = FALSE)


nclust <-5

#use weighed modes
prof_clust <- kmodes(vprof_sample, iter.max = 100, modes = nclust, weighted = TRUE) 
clust<-array(as.numeric(unlist(prof_clust$modes)), dim = c(nclust, 30))

#save file
saveRDS(prof_clust, file = "./kmodes-clust5-1.RDS")

cols <- brewer.pal(name="Set3", n=11)
col_2class <- cols[5:6]
col_3class <- cols[5:7]
col_10class <- cols <- brewer.pal(name="Paired", n=10)

#plot
pdf("~/projects/screim/plots/modes5_weighted.pdf", height=5, width=6)
par(mfrow=c(2,1))
par(oma=c(1, 1, 1, 1), mar=c(0.5, 4, 1, 1), las=1)
par(cex=1.0)
image2D(clust, x=1:nclust, y=1:30, breaks = 0:2, col = col_2class,
        main="Verticle Modes", ylab="Levels [500m]", colkey = FALSE)
#axis(side=1, labels = c("high St", "Cg", "Cb + Anvil", "low St", "nSt"), at = 1:nclust)

grid(lty=1, lwd=3, nx = nclust, ny=0, col="white")

frame()
legend("top", legend = rev(c("Non-convective", "Convective")),bty="n", 
       fill = rev(col_2class), horiz = TRUE)
dev.off()


