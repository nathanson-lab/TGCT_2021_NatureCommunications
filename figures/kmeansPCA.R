rm(list = ls())
setwd("~/Documents/nathansonlab/tecac-manuscript/PCA/")
# pluta 4/7/2017
# determine caucasian v non-caucasian with k-means clustering and plot

library(ggplot2)
library(dplyr)
library(broom)


# function to get distance between two points on a cartesian plane
getDist <- function(x1, y1, x2, y2)
{
  d = sqrt(  ((x2 - x1)^2) + ((y2 - y1)^2)    )
  return(d)
}



# hapmap demographics
HAP.dem = read.table("relationships_w_pops_051208.txt", header=TRUE)

# eigenvectors
ev = paste("EV", seq(1:10), sep="")

evec.dat = read.table("merged.evec", header=FALSE, col.names= c("IID", ev, "Pheno"))
levels(evec.dat$Pheno) <- c("Parent", "Case", "Control")
evec.dat$PRE = substr(evec.dat$IID, 1, 2)
evec.dat$RGN <- "HAPMAP"

# group by geographical region
evec.dat$RGN[evec.dat$PRE %in% c("MD", "YU", "FH", "UP", "CA", "PM")] <- "NA"
evec.dat$RGN[evec.dat$PRE %in% c("GR", "HR", "PV", "TR", "UK")] <- "EUR"
evec.dat$RGN[evec.dat$PRE %in% c("GR", "NJ", "NO", "SE")] <- "SCD"

dat = evec.dat
fam <- read.table("TECAC_ALL_QC8.fam", header = F)
dat$Sex <- 1
dat$Sex[ fam$V2 %in% dat$IID ] <- fam$V5[ match(dat$IID, fam$V2 )]


set.seed(4717)

# from: https://cran.r-project.org/web/packages/broom/vignettes/kmeans.html
kclusts <- data.frame(k=4) %>% group_by(k) %>% do(kclust=kmeans(dat[,2:3], .$k, nstart=20))
clusters <- kclusts %>% group_by(k) %>% do(tidy(.$kclust[[1]]))
assignments <- kclusts %>% group_by(k) %>% do(augment(.$kclust[[1]], dat[,2:3]))
clusterings <- kclusts %>% group_by(k) %>% do(glance(.$kclust[[1]]))
levels(assignments$.cluster) <- as.factor(c( "1", "2", "3", "4", "0"))
set.seed(4717)
ev.cls <- kmeans( dat[,2:3], 4, nstart= 20)
t = table(ev.cls$cluster, dat$RGN)

# largest cluster
l = which(ev.cls$size == max(ev.cls$size))

# the center of each cluster
ev.cls$cluster <- as.factor(ev.cls$cluster)
dat$cluster = ev.cls$cluster

cnt.dat = data.frame(x=ev.cls$centers[,1], y=ev.cls$centers[,2], cluster = c(1:4))

# cluster  center
cls1.cnt = ev.cls$centers[l,]

dat$dist = getDist(cls1.cnt[1], cls1.cnt[2], dat$EV1, dat$EV2)
d.sd   = sd(dat$dist[dat$cluster == l])

# points within the region (caucasian subjects)
ind = which(dat$dist <= 6*d.sd)


# radius of the inclusion region is # of standard deviations of distance from
# the center of the region
r = d.sd * 6
xc = cnt.dat$x[l]
yc = cnt.dat$y[l]


# counts for cluster, in region, not in region
t2 = t[,-which(colnames(t) == "HAPMAP")]
n.cls = sum(t2[l,])
n.cls.not = sum(t2[-l,])

# 162 is the number of caucasian hapmap samples
in.rgn = length(ind) - 162
in.rgn.not = n.cls - in.rgn

cases <- read.table("non-cauc-cases.txt", header = F)
parents <- read.table("non-cauc-parents.txt", header = F)

#levels(ev.cls$cluster) <- c("Caucasian", "Hispanic",  "Asian", "African" )


p3 <- ggplot(dat, aes(x=EV1, y=EV2, color= ev.cls$cluster)) +  
geom_point(alpha=0.3) +
  
  geom_point(data=dat[ind,], aes(x=EV1, y=EV2), col="red") +
  annotate("path", x = xc + r*cos(seq(0,2*pi, length.out=100)),
                   y = yc + r*sin(seq(0,2*pi, length.out=100))) +
  geom_point(data=cnt.dat, aes(x=x, y=y), col = "black") +
#  scale_color_manual(values= c("Caucasian" = "green", 
#                               "African" = "orange", 
#                               "Hispanic" = "purple", 
 #                              "Asian" = "blue"), aesthetics = c("fill", "color")) +
  ggtitle("PCA") + 
  guides(colour = guide_legend(override.aes = list(alpha = 1), title = "Ancestry")) +
  theme_minimal() 

png("k-means-PCA.png", height=800, width=800)
print(p3)
dev.off()




write.table(dat$IID[ind], "caucasians.txt", append=FALSE, quote=FALSE, row.names=FALSE,
            col.names=FALSE)
write.table(dat$IID[dat$cluster == 2], "hispanics.txt", append=FALSE, quote=FALSE, row.names=FALSE,
            col.names=FALSE)
write.table(dat$IID[dat$cluster == 4], "africans.txt", append=FALSE, quote=FALSE, row.names=FALSE,
            col.names=FALSE)
write.table(dat$IID[dat$cluster == 3], "asians.txt", append=FALSE, quote=FALSE, row.names=FALSE,
            col.names=FALSE)
temp <- dat[-ind,]
write.table(temp$IID[temp$cluster == 1], "mixedrace.txt", append=FALSE, quote=FALSE, row.names=FALSE,
            col.names=FALSE)













