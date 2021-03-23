#!/usr/bin/env Rscript

# pluta 3/21/18

# conditional analysis conditioning on > 1 snp


# args come from conditionalAnalysisJoint.sh
args = commandArgs(trailingOnly = TRUE)
if( length(args) < 2)
{
  stop("need to provide 2 argument: snps-pvals  cojofile")
}

library(ggplot2)

refsnpfile = args[1]

ref <- read.table(refsnpfile, header=F, col.names = c("snp", "p"), colClasses = 
c("character", "numeric"))

ref$bp <- as.integer(unlist(lapply(strsplit(ref$snp, ":"), function(x) x[2])))
ref.chr <- as.integer(strsplit(ref$snp[1], ":")[[1]][1]) # this is a constant



# these files are generated from conditionalAnalysis.sh
COJOFILE <- args[2]

# not sure how to incorporate LD just yet
#SNPFILE <- paste(snpname, "snp.ld", sep = ".")
#LDFILE <- paste(snpname, "r.ld", sep = ".")


dat <- read.table(COJOFILE, header=T, 
                  col.names = c("Chr", "SNP", "bp", "refA", "freq", "b", 
                                "se", "p", "n", "freq_geno", "bC", "bC_se", "pC"
), 
                  colClasses = c("character", "character", "integer", "character
", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
                  "numeric", "numeric", "numeric"))

#tmp <- strsplit(COJOFILE, "\\.")[[1]][1]
#tmp <- strsplit(tmp, "_")[[1]]

refsnp.snp <- ref$snp[1]
refsnp.bp <- ref$bp[1]
refsnp.p <- -log10( ref$p[1])


#ld <- t(read.table(LDFILE, header=F))
#snp <- t(read.table(SNPFILE, header=F))

# combine ld to each snp
#ld.mat <- data.frame(snp <- snp[,1], r2 <- (ld[,1])^2)
#colnames(ld.mat) <- c("snp", "r2")

# put r2 values in the main dataframe
#dat$r2 <- 0
#dat$r2 <- ld.mat$r2[match(dat$SNP, ld.mat$snp)]
#dat$r2[is.na(dat$r2)] <- 0

# might need to tune this
dat$indp <- !(dat$p < 0.00001 & dat$pC > 0.00001)

if( any( is.na(dat$pC) ) )
{
	dat$indp[which(is.na(dat$pC))] <- FALSE
}

# top dependent snp
ind <- which(dat$p == min(dat$p[!dat$indp]))
topsnp.bp <- dat$bp[ind]
topsnp.p <- -log10(dat$p[ind])
topsnp.snp <- dat$SNP[ind]

# top independent snp
ind <- which(dat$p == min(dat$p[dat$indp]))
topindsnp.bp <- dat$bp[ind]
topindsnp.p <- -log10(dat$p[ind])
topindsnp.snp <- dat$SNP[ind]



xmin <- min(c(dat$bp, topindsnp.bp, topsnp.bp, ref$bp))
xmax <- max(c(dat$bp, topindsnp.bp, topsnp.bp, ref$bp))
offset <- (xmax - xmin) / 20

ymax <- max(c(-log10(dat$p), -log10(ref$p), topindsnp.p, topsnp.p))

p1 <- ggplot( data = dat, aes(x=bp, y=-log10(p), colour=indp)) +
      theme_light() +
      theme(plot.title = element_text(size = 28, face = "bold"),
	    legend.text = element_text(size = 18, face = "bold"),
            legend.title = element_text(size = 22, face = "bold"), 
            axis.text.x = element_text(size = 16),
            axis.title.x = element_text(size = 20),
            axis.title.y = element_text(size = 20),
            axis.text.y = element_text(size = 16)) +
      ylim(0, ymax + 0.5) +
      geom_point(size=4) +
      geom_hline(yintercept = -log10(1e-5), linetype="dashed", colour="blue") +
      annotate("text", x = round(topsnp.bp + offset), y = topsnp.p, label = topsnp.snp, size = 8) +
      annotate("text", x = round(topindsnp.bp + offset), y = topindsnp.p, label = topindsnp.snp, size = 8, 
hjust = -0.5) +
      geom_point(x = refsnp.bp, y = refsnp.p, size = 4, shape = 24, colour = "black", fill = "blue") +
      annotate("text", x = refsnp.bp + offset, y = refsnp.p, label = refsnp.snp, size= 8) +
      ggtitle(paste(ref.chr, "Joint Conditional Plot", sep = " ")) 

# add each of the conditioned snps
for(i in 1:dim(ref)[1])
{
  p1 <- p1 + geom_point(x = ref$bp[i], y = -log10(ref$p[i]), colour="black", fill="blue", shape = 25, size 
= 4)
} 


# format date for output
x <- c(1,3,5)
dt <- strsplit(date(), " ")[[1]][x]
dt <- paste(dt[1], dt[2], dt[3], sep="_")
# ensure a unique output name
OUTNAME=paste(ref.chr, ref$bp[1], dt, "jt_cond.png", sep="_")

if(file.exists(OUTNAME))
{
	i=2
	repeat 
	{
		OUTNAME = paste(ref.chr, ref$bp[1], dt, "jt_cond", i, ".png", sep="_")
		if( !file.exists(OUTNAME) )
		{
			break
		} else
		  { 
			i = i + 1
		  }
	}	
}

print(paste("Top independent snp: ", topindsnp.snp, sep=""))
print(paste("Top dependent snp: ", topsnp.snp, sep=""))
print("writing image...")


png(OUTNAME, width=1000, height=1200)
print(p1)
dev.off()
