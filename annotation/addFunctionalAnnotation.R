#!/usr/bin/env Rscript

# pluta 11/7/19

# this script was written ad-hoc as analysis developed, a finalized version
# will be much more organized
#rm(list = ls())

#args = commandArgs(trailingOnly = TRUE)
#if( length(args) < 1)
#{
#  stop('need to provide arguments: FILENAME, the reference snp of the credset')
#}

INFILENAME=args[1]
len <- nchar(INFILENAME)
if( substr(INFILENAME, len, len) == "/" )
{
  INFILENAME <- substr(INFILENAME, 1, len - 1)
}
# each locus is computed independently; run script per locus 

# ENCSR418RNI; ENCSR886NTH; ENCSR303XKE; ENCSR898RGU all annotate identically
# drop three of these

# ================== preprocesser =================== #
# ---------------------------- libraries ------------- #
library(data.table)
library(ggplot2)
library(IRanges)
library(BiocManager)

setwd("/Users/johnpluta/Documents/nathansonlab/tecac-manuscript/annotation")
# --------------------------------------------------- #


# ----- header definitions ----- #
methyl.bed.header <- c("Chrom", "Start", "End", "Name", "Score", "strand","Start2", "End2", "rgb", "count", "percentMeth")
narrowpeak.bed.header <- c("Chrom", "Start", "End", "name", "Score", "strand", "signalVal", "-log10(p)", "-log10(q)", "peak")
broadpeak.bed.header <- c("Chrom", "Start", "End", "name", "Score", "strand", "signalVal", "p", "q")
cis.reg.header <- c("Chrom", "Start", "End", "Cis.Reg.Name", "Score", "strand", "Start2", "End2", "color")
chr.in.order <- c(paste("chr", c(seq(1:22), "X"), sep = ""))
# ----------------------------- #

# =================================================== #




# ================== functions =========================== #

# ---------------- getLRT.p ----------------------------- #
# get the p-value for the likelihood ratio test between the baseline
# model and the model with selected covariate(s)
getLRT.p <- function( logBF.Base, logBF.Model )
# input: logBF.Base, logBF.model (numeric): log(baye's factor) from 
#         the two models- this comes from PAINTOR
# output: p-value (numeric), p-value for the LRT between the two models
{
  # LRT ~ x^2 with 1df
  return( 1 - pchisq( -2 * (logBF.Base - logBF.Model), 1))
}
# -------------------------------------------------------- #


# ----------------------- inEqtl ------------------------- #
# check if a crv is in an eqtl
inEqtl <- function( pos, eqtl )
  # input: pos, snp position (hg38)
  #         eqtl: eqtl position
  # output: ind (logical), logical vector of which eqtls are in the crv
{
  return( eqtl$LD.block.end >= pos & eqtl$LD.block.start <= pos )
}
# -------------------------------------------------------- #


# -------------------- snpAnnotSimple -------------------- #
# simple version of annotate snp that returns a logical vector
# no longer need to convert the output to logical
snpAnnotSimple <- function( pos, bed, chr )
  # input: pos (integer), snp poistion
  #        bed (data.frame), data.frame of the bed file that is 
  #         being annotated from
  #        chr (integer), the chromosome of interest
  # output: ind (logical), a vector of which crvs overlap with
  #       annotation
{
  bed <- bed[ bed$V1 == chr,]
  rangeA <- IRanges( pos, pos )
  rangeB <- IRanges(bed$V2, bed$V3)
  ind <- findOverlaps(rangeA, rangeB, type = "within")
  return(ind@from)
}
# -------------------------------------------------------- #

# ----------------------- replaceNA ---------------------- #
# annotation matrix cannot have NA or paintor will crash
# for annotation, it is sufficient to replace NA with 0
replaceNA <- function(x)
{
  x[is.na(x)] <- 0
  return(x)
}
# --------------------------------------------------------- #

# --------------- joinMethylReplicates ------------------- #
# in methylation data, if there are two sets of data, truncate
# data to only those that appear in both sets. returns a single set of
# data
joinMethylReplicates <- function( dat1, dat2, chr = NULL)
  # dat1 (data.frame), bed file of replicate 1
  # dat2 (data.frame), bed file or replicate 2
  # chr (integer), if chr is NULL, run genome wide; else subset
{
  out <- c()
 
  if( is.null(chr))
  {
    tmp1 = dat1
    tmp2 = dat2
  } else
  {
    tmp1 <- dat1[ dat1$Chrom == chr, ]
    tmp2 <- dat2[ dat2$Chrom == chr, ]
  }
  tmp.out <- merge(tmp1, tmp2, by.x = "Start", by.y = "Start")

  tmp.out <- tmp.out[ !is.na(tmp.out$Chrom.x),]
  if( "percentMeth.x" %in% colnames(tmp.out))
  {
    tmp.out <- tmp.out[,colnames(tmp.out) %in% c("Start", "Chrom.x", "End.x", "percentMeth.x")]
    colnames(tmp.out) <- c("Start", "Chrom", "End", "percentMeth")
  } else
  {
    tmp.out <- tmp.out[,colnames(tmp.out) %in% c("Start", "Chrom.x", "End.x", "q.x")]
    colnames(tmp.out) <- c("Start", "Chrom", "End", "q")
  }
  
  
  out <- rbind(tmp.out, out)
  
  
  return(out)
}
# ------------------------------------------------------------#


# ---------------- attachSnpAnnotMethyl ------------------- #
# attach features based on snp id/position
attachSnpAnnotMethyl <- function( pos, dat, chr, varname )
# only slightly different than attachSnpAnnot, to account for different
# formatting in methylation files; probably a smarter way to do this,
# or roll both functions into attachSnpAnnotSimple if we only care about
# binary values
{
  if( is.null(pos) )
  {
    stop("pos has null value, did you pass the right attribute?")
  }
  # subset the functional data so we only pull snps in the same chromosome
  # as start
  dat <- dat[dat$Chrom == chr,]
  if( !(varname %in% colnames(dat)))
  {
    print(paste(varname, " not found in data.", sep = ""))
    stop()
  }
  # find snps in range
  ind <- dat$start %in% pos
  
  if( length(ind) > 0)
  {
    x <- dat[[varname]][ind]
    
    # one snp can map to multiple genes (?); need to concatenate into one var
    if(length(x) > 1)
    {
      print(paste("Start position:", pos, "had", length(x), "matches", sep = " "))
      x <- paste(dat[[varname]][ind], collapse = ";")
    }
    
    return( x )
  }
  
  return(NA)
}
# ----------------------------------------------------------- #




# ---------------- attachSnpAnnot ------------------- #
# attach features based on snp id/position
attachSnpAnnot <- function( pos, dat, chr, varname )
  # input:
  # pos (integer), position of snp from the credible set
  # dat (data.frame), the bed file data
  # chr (integer), chromosome to subset by
  # varname (string), which attribute do we want to return?
  #     this only matters if not using a binary value
  #
  # output: return varname values that are in CRV regions
{
  if( is.null(pos) )
  {
    stop("pos has null value, did you pass the right attribute?")
  }
  # subset the functional data so we only pull snps in the same chromosome
  # as start
  dat <- dat[dat$Chrom == chr,]
  
  
  if( !(varname %in% colnames(dat)))
  {
    print(paste(varname, " not found in data.", sep = ""))
    stop()
  }
  
  
  # find snps in range
  ind = which(dat$Start <= pos & dat$End >= pos)
  
  if( length(ind) > 0)
  {
    x <- dat[[varname]][ind]
    
    # one snp can map to multiple genes (?); need to concatenate into one var
    if(length(x) > 1)
    {
      print(paste("Start position:", pos, "had", length(x), "matches", sep = " "))
      x <- paste(dat[[varname]][ind], collapse = ";")
    }
    
    return( x )
  }
  
  return(NA)
}
# ------------------------------------------------------------- #

# ------------------------ readMethylIdat --------------------- #
# read methylation data in .idat format
readMethylIdat <- function( IDATNAME )
{
  # reads in two idats with the same rootname
  # eg ROOT_Red.idat and ROOT_Grn.idatu
  idats <- c(IDATNAME)
  rgset <- read.metharray(idats, verbose = T)
  mset <- preprocessIllumina(rgset)
  mset <- mapToGenome(mset)
  df <- data.frame( chr = rep(mset@rowRanges@seqnames@values, mset@rowRanges@seqnames@lengths),
                    start = mset@rowRanges@ranges@start,
                    name = rownames(mset) )
  colnames(df) <- c("chr", "start", "name")
  return(df)
}
# ------------------------------------------------------------ #

# --------------------------- readBed ------------------------ #
readBed <- function( BEDFILE, bedheader )
  # read in bedfile, attach header, reorder by chromosome and then position
{
  dat <- fread(BEDFILE, header = F)
  dat <- as.data.frame(dat)
  colnames(dat) <- bedheader
  dat$Chrom <- factor(dat$Chrom, chr.in.order)
  dat <- dat[order(dat$Chrom, dat$Start),]
  
  # qvalue is p-value adjusted for multiple comparisons
  # to get p-value from -log10: 10^-p
  if( "-log10(p)" %in% colnames(dat))
  {
    dat$p <- 10^-(dat$`-log10(p)`)
  }
  
  if( "-log10(q)" %in% colnames(dat))
  {
    dat$q <- 10^-(dat$`-log10(q)`)
  }
  return(dat)
}
# ---------------------------------------------------------- #

# ----------------------- convertToBinIn-------------------- #
convertToBinInt <- function( var )
  # convert data of any type to a binary integer (e.g. is a chromatin feature present or absent)
  # input: var, a column of data from the annotation table
  # output: var, the same data converted to binary int
{
  var[ is.na(var) ] <- 0
  var[ var != 0 ] <- 1
  var <- as.integer(var)
  return(var)
}
# ---------------------------------------------------------- #

# ------------------------ convertToBed -------------------- #
convertToBED <- function( dat )
  # convert idat data to BED
{
  if( all(c("Start", "Chrom", "End") %in% colnames(dat)))
  {
    out <- data.frame( Chrom = dat$Chrom, Start = dat$Start, End = dat$End, percentMeth = dat$percentMeth)
    colnames(out) <- c("Chrom", "Start", "End", "percentMeth")
    return(out)
  }
  if( substr(dat$chr[1], 1, 1) != "c")
  {
    dat$chr <- paste("chr", as.character(dat$chr), sep ="")
  }
  
  out <- data.frame( chrom = dat$chr, chromStart = dat$start - 1, 
                     chromEnd = dat$start, name = dat$name)
  colnames(out) <- c("chrom", "chromStart", "chromEnd", "name")
  return(out)
}
# ---------------------------------------------------------- #

# ======================= end functions ================================ #



newhits <- c("1:212449403",
  "2:111927379",
             "5:1280128",
             "6:32032421",
             "6:33533625",
             "8:120933963",
             "9:779507",
             "9:127190340",
             "9:140073294",
             "10:7534248",
             "11:30351223",
             "12:1051495",
             "12:51301431",
             "12:53793209",
             "17:76691564",
             "18:692095",
             "19:28356614",
             "20:52197366",
             "X:100432681",
             "X:153535143",
             "X:24384181",
             "X:66489986")


# ========================= MAIN ========================== #

print("Parsing credSet file...")

# the CRV file already has hg19 and hg38 positions attached

credSetFile <- paste(INFILENAME, "/", INFILENAME, ".credSet", sep = "")
crv.dat <- read.table(credSetFile, header = TRUE, as.is= TRUE)
crv.dat$h19pos <- as.integer(unlist(lapply(strsplit(crv.dat$SNP.ID.hg19, ":"), function(x) x[2])))
crv.dat$h38pos <- as.integer(unlist(lapply(strsplit(crv.dat$SNP.ID.hg38, ":"), function(x) x[2])))

chr = paste("chr", crv.dat$chr[1], sep = "")
crv.dat$z <- crv.dat$Effect / crv.dat$se
print("done!")

# write locus file
print("writing locus file...")
LOCUSFILE <- paste(INFILENAME, INFILENAME, sep = "/")
write.table(crv.dat[,colnames(crv.dat) %in% c("chr", "SNP.ID.hg19", "rsID", "P", "z", "h19pos" )], 
 LOCUSFILE, col.names = T, row.names = F, quote = F)
print("done!")

print("get features from biomart...")
# get snp and gene mappings
# germline data
grch37.snp = useMart(biomart="ENSEMBL_MART_SNP", host="grch37.ensembl.org",dataset="hsapiens_snp")

#Mart used to map Ensembl Gene IDs to Gene name
grch37 = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")



snp.dat <- getBM(attributes = c("refsnp_id", "ensembl_gene_stable_id", "associated_gene", "consequence_type_tv"), 
                filters = "snp_filter", 
                values = crv.dat$rsID, 
                mart = grch37.snp)

# associated gene doesnt map to every instance of ensembl gene, and some
# annotations are missing- use this list from ENSEMBL and dbsnp
snp.dat$associated_gene <- as.character(snp.dat$associated_gene)

# local db of enseml genes and associated gene name
en.list <- c("ENSG00000152219", "ENSG00000186452", "ENSG00000261824",
             "ENSG00000261770", "ENSG00000267575", "ENSG00000272635",
             "ENSG00000267389", "ENSG00000266977", "ENSG00000267623",
             "ENSG00000267630", "ENSG00000234915", "ENSG00000066027",
             "ENSG00000065600", "ENSG00000234915", "ENSG00000115112",
             "ENSG00000124374", "ENSG00000114850", "ENSG00000108669",
             "ENSG00000267623", "ENSG00000092345", "ENSG00000173890",
             "ENSG00000268220", "ENSG00000173889", "ENSG00000163558",
             "ENSG00000109323", "ENSG00000109332", "ENSG00000246560",
             "ENSG00000145354", "ENSG00000164037", "ENSG00000164038",
             "ENSG00000164039", "ENSG00000138778", "ENSG00000248740",
             "ENSG00000138769", "ENSG00000138757", "ENSG00000246541",
             "ENSG00000163104", "ENSG00000163106", "ENSG00000164362",
             "LRG_343", "ENSG00000204344", "ENSG00000234947",
             "ENSG00000226257", "ENSG00000226033", "ENSG00000206342",
             "ENSG00000272295", "ENSG00000206338", "ENSG00000236250",
             "ENSG00000204344", "ENSG00000234947", "ENSG00000168477",
             "ENSG00000213676", "ENSG00000233323", "ENSG00000234539",
             "ENSG00000231608", "ENSG00000168468", "ENSG00000206258",
             "ENSG00000229353", "ENSG00000228628", "ENSG00000112514",
             "ENSG00000197283", "ENSG00000242014", "ENSG00000226492",
             "ENSG00000227460", "ENSG00000245330", "ENSG00000176058",
             "ENSG00000176101", "ENSG00000176248", "ENSG00000176884",
             "ENSG00000261793", "ENSG00000002016", "ENSG00000170374",
             "ENSG00000185591", "ENSG00000135409", "ENSG00000257379",
             "ENSG00000205352", "ENSG00000197111", "ENSG00000139625",
             "ENSG00000139546", "ENSG00000267281", "ENSG00000170653",
             "ENSG00000135390", "ENSG00000176105", "ENSG00000132199",
             "ENSG00000265490" ,"ENSG00000266171", "ENSG00000261824",
             "ENSG00000261770", "ENSG00000102384", "ENSG00000268013",
             "ENSG00000007350", "ENSG00000269329", "ENSG00000196924",
             "ENSG00000102080", "ENSG00000185254", "ENSG00000223731",
             "ENSG00000226280", "ENSG00000169083", "ENSG00000204581",
             "ENSG00000153093", "ENSG00000153094", "ENSG00000049656",
             "ENSG00000132570", "ENSG00000152705", "ENSG00000224186",
             "ENSG00000069011", "ENSG00000187678", "ENSG00000231185",
             "ENSG00000235168", "ENSG00000055211", "ENSG00000237502",
             "ENSG00000186625", "ENSG00000131023", "ENSG00000120253",
             "ENSG00000120265", "ENSG00000219433", "ENSG00000231760",
             "ENSG00000120256", "ENSG00000268592", "ENSG00000217733",
             "ENSG00000164520", "ENSG00000223701", "ENSG00000216906",
             "ENSG00000203722", "ENSG00000096433", "ENSG00000030110",
             "ENSG00000197251", "ENSG00000204188", "ENSG00000002822",
             "ENSG00000176349", "ENSG00000147596", "ENSG00000137090",
             "ENSG00000259290", "ENSG00000118369", "ENSG00000033327",
             "ENSG00000244573", "ENSG00000121316", "ENSG00000126775",
             "ENSG00000182521", "ENSG00000166450", "ENSG00000259180",
             "ENSG00000166938", "ENSG00000260773", "ENSG00000075131",
             "ENSG00000169032", "ENSG00000261351", "ENSG00000174446",
             "ENSG00000174444", "ENSG00000174442", "ENSG00000188501",
             "ENSG00000262117", "ENSG00000171490", "ENSG00000263307",
             "ENSG00000103342", "ENSG00000261560", "ENSG00000234719",
             "ENSG00000156968", "ENSG00000183793", "ENSG00000205423",
             "ENSG00000259843", "ENSG00000260381", "ENSG00000261170",
             "ENSG00000260539", "ENSG00000090863", "ENSG00000168411",
             "ENSG00000168404", "ENSG00000263456", "ENSG00000108753",
             "ENSG00000259549", "ENSG00000160321", "ENSG00000269615",
             "ENSG00000271095", "ENSG00000269504", "ENSG00000229676",
             "ENSG00000268981", "ENSG00000198153", "ENSG00000268789",
             "ENSG00000213973", "ENSG00000269509", "ENSG00000268696",
             "ENSG00000269067", "ENSG00000271661", "ENSG00000261558",
             "ENSG00000261615", "ENSG00000267886", "ENSG00000183850",
             "ENSG00000213096", "ENSG00000197372", "ENSG00000269289",
             "ENSG00000213967", "ENSG00000229000", "ENSG00000020256",
             "ENSG00000228404", "ENSG00000160285", "ENSG00000223901",
             "ENSG00000215424", "ENSG00000160294", "ENSG00000228137",
             "ENSG00000239415", "ENSG00000182362", "ENSG00000160298",
             "ENSG00000160299", "ENSG00000223692", "ENSG00000160305",
             "ENSG00000099949", "ENSG00000265148", "ENSG00000213246",
             "ENSG00000108375", "ENSG00000176160", "ENSG00000108389",
             "ENSG00000264672", "ENSG00000108387", "ENSG00000181013", 
             "ENSG00000121101", "ENSG00000212195", "ENSG00000108384",
             "ENSG00000175175", "ENSG00000263938", "ENSG00000108395",
             "ENSG00000224738", "ENSG00000182628", "ENSG00000232160",
             "ENSG00000076770", "ENSG00000171004", "ENSG00000123728",
             "ENSG00000079313", "ENSG00000129911", "ENSG00000267141")
gene.list <- c("ARL14EP", "TMPRSS12", "LINC00662",
               "AC006504.1", "CTC-459F4.3", "LLNLF-65H9.1",
               "AC006504.4", "AC006504.2", "not found", 
               "AC005758.1", "AL360091.3", "PPP2R5A",
               "TMEM206", "RP11-384C4.7", "TFCP2L1",
               "PAIP2B", "SSR3", "CYTH1", "AC005357.2", "DAZL",
               "GPR160", "ENSG00000268220","PHC3", "PRKCI",
               "MANBA", "UBE2D3", "UBE2D3-AS1", "CISD2",
               "SLC9B1", "SLC9B2", "BDH2", "CENPE",
               "LINC02428", "CDKL2", "G3BP2", "AC096746.1",
               "SMARCAD1", "HPGDS", "TERT", "TERT",
               "STK19", "STK19", "STK19", "STK19", "STK19",
               "DAQB-331", "CYP21A2", "STK19", "STK19",
               "STK19", "TNXB", "ATF6B", "TNXB", "ATF6B", "TNXB",
               "ATF6B", "TNXB", "TNXB", "ATF6B",
               "CUTA", "SYNGAP1", "RN7SL26P", "CUTA",
               "SYNGAP1", "AP005717.1", "TPRN",
               "SSNA1", "ANAPC2", "GRIN1", "AL929554.1", "RAD52",
               "SP7", "SP1", "AMHR2", "AC023509.1",
               "PRR13", "PCBP2", "MAP3K12", "TARBP2", "ATF7-NPFF",
               "ATF7", "ATP5MC2", "YES1", "ENOSF1", "RP11-806L2.6",
               "AP0010203.5", "LINC00662", "AC006504.1", "CENPI",
               "TKTL1", "TKTL1", "FLNA", "FLNA", "TEX28", "TEX28",
               "SUPT20HL1", "AL049641.1", "AR", "ACOXL-AS1", "ACOXL",
               "BCL2L11", "CLPTM1L", "PCBD2", "CATSPER3", "C5orf66",
               "PITX1", "SPRY4", "SPRY4-AS1", "AL078581.2", "GINM1",
               "RP1-12G14.6", "KATNA1", "LATS1", "NUP43", "PCMT1",
               "BTBD10P2", "AL355312.2", "LRP11", "AL355312.3",
               "CCT7P1", "RAET1E", "RAET1E-AS1", "AL355312.1",
               "RAET1G", "ITPR3", "BAK1", "LINC00336", "GGNBP1", "MAD1L1",
               "AC104129.1", "PRDM14", "DMRT1", "RP11-687M24.7", "USP35",
               "GAB2", "RPL30P11", "PLBD1", "ATG14", "TBPL2", "PRTG", 
               "AC012378.1", "DIS3L", "RP11-352G18.2", "TIPIN",
               "MAP2K1", "AC116913.1", "SNAPC5", "RPL4", "ZWILCH",
               "LCTL", "BCAR4", "RSL1D1", "AC007216.4", "GSPT1",
               "AC007216.3", "NPIPB2", "MPV17L", "NPIPA5", "CNEP1R1",
               "AC007610.1", "RP11-429P3.5", "AC009053.3", "RP11-252A24.7", "GLG1",
               "RFWD3", "MLKL", "MIR5189", "HNF1B", "RP11-115K3.1",
               "ZNF208", "AC003973.1", "BNIP3P28", "AC003973.4", "ZNF492",
               "AC024563.1", "ZNF849P", "VN1R87P", "ZNF99", "BNIP3P34",
               "ZNF723", "ZNF728", "BNIP3P36", "LINC01859", "LINC01858", "AC074135.1",
               "ZNF730", "ZNF254", "ZNF675", "AC011503.1", "ZNF726",
               "SEPT7P8", "ZFP64", "AP001468.1", "LSS", "AP001469.1",
               "MCM3AP-AS1", "MCM3AP", "AP001469.2", "AP001469.9",
               "YBEY", "C21orf58", "PCNT", "DIP2A-IT1", "DIP2A", "LZTR1",
               "TSPOAP1-AS1", "SUPT4H1", "RNF43", "HSF5", "MTMR4",
               "SEPTIN4-AS1", "SEPTIN4", "C17orf47", "TEX14", "U3",
               "RAD51C", "PPM1E", "GC17M058973", "TRIM37",
               "AC099850.1", "SKA2", "RAP2C-AS1", "MBNL3", "HS6ST2", "RAP2C",
               "REXO1", "KLF16", "AC012615.4")
gene.key <- data.frame(ensembl = en.list, gene.id = gene.list)

ind = which(snp.dat$ensembl_gene_stable_id != "" & snp.dat$associated_gene != "")

#if( length(ind) > 0)
#{
#  tmp <- data.frame( ensembl <- unique(snp.dat$ensembl_gene_stable_id[ind]),
#                     gene.id <- unique(snp.dat$associated_gene[ind]) )
#  colnames(tmp) <- c("ensembl", "gene.id")
#  gene.key <- rbind(gene.key, tmp)
#}

if( any(gene.key$gene.id == ""))
{
  print("gene.key is missing some definitions- enter these manually into gene.key and rerun")
  stop()
}

  
# map gene data to snp data
ind = match(snp.dat$ensembl_gene_stable_id, gene.key$ensembl)
snp.dat$associated_gene[!is.na(ind)] <- as.character( gene.key$gene.id[ind[!is.na(ind)]] )
  
  
gene.dat <- getBM(attributes = c("ensembl_gene_id", "5_utr_start", "5_utr_end",
                                   "3_utr_start", "3_utr_end", "exon_chrom_start", "exon_chrom_end"),
                    filters = "ensembl_gene_id", 
                    values =  snp.dat$ensembl_gene_stable_id, 
                    mart = grch37)
  
  
# add exon, 3' utr, 5' utr; snps within range
exonRange <- utr5Range <- utr3Range <- c()
  
for( gene in unique(gene.dat$ensembl_gene_id))
{
   exonRange <- c(IRanges(gene.dat$exon_chrom_start[gene.dat$ensembl_gene_id == gene],
                           gene.dat$exon_chrom_end[gene.dat$ensembl_gene_id == gene]), exonRange)
    
  tmp <- gene.dat[ !is.na(gene.dat$`5_utr_start`),]
  utr5Range <- c(IRanges(tmp$`5_utr_start`[tmp$ensembl_gene_id == gene],
                           tmp$`5_utr_end`[tmp$ensembl_gene_id == gene]), utr5Range)
    
  tmp <- gene.dat[ !is.na(gene.dat$`3_utr_start`),]
  utr3Range <- c(IRanges(tmp$`3_utr_start`[tmp$ensembl_gene_id == gene],
                           tmp$`3_utr_end`[tmp$ensembl_gene_id == gene]), utr3Range)
}
  
crv.dat$exons <- 0
crv.dat$utr5 <- 0
crv.dat$utr3 <- 0
 
if( !any(is.null(c(exonRange, utr5Range, utr3Range) )))
{
      crv.dat$exons[findOverlaps(crv.dat$h19pos, exonRange)@from] <- 1
      crv.dat$utr5[findOverlaps(crv.dat$h19pos, utr5Range)@from] <- 1
      crv.dat$utr3[findOverlaps(crv.dat$h19pos, utr3Range)@from] <- 1
}

  
  
results <- merge(snp.dat, gene.dat, by.x = "ensembl_gene_stable_id", by.y = "ensembl_gene_id", all.x=T)
out <- merge(crv.dat, results, by.x = "rsID", by.y = "refsnp_id", all.x = T)
print("done!")
  
print("adding gene expression data...")
  # gene expression ------
  # column 1: gene_id - gene name of the gene the transcript belongs to (parent gene). If no gene information is provided, gene_id and transcript_id is the same.
  # column 2: transcript_id(s) - transcript name of this transcript
  # column 3: length - the transcript's sequence length (poly(A) tail is not counted)
  # column 4: effective_length - the length containing only the positions that can generate a valid fragment
  # column 5: expected_count - the sum of the posterior probability of each read comes from this transcript over all reads
  # column 6: TPM - transcripts per million, a measure of relative measure of transcript abundance
  # column 7: FPKM - fragments per kilobase of transcript per million mapped reads, another relative measure of transcript abundance
  # column 8: posterior_mean_count - posterior mean estimate calcualted by RSEM's Gibbs sampler
  # column 9: posterior_standard_deviation_of_count - posterior standard deviation of counts
  # column 10: pme_TPM - posterior mean estimate of TPM
  # column 11: pme_FPKM - posterior mean estimate of FPKM
  # column 12: TPM_ci_lower_bound - lower bound of 95% credibility interval for TPM values
  # column 13: TPM_ci_upper_bound - upper bound of 95% credibility interval for TPM values
  # column 14: FPKM_ci_lower_bound - lower bound of 95% credibility interval for FPKM values
  # column 15: FPKM_ci_upper_bound - upper bound of 95% credibility interval for FPKM values
dat4 <- as.data.frame(fread("ENCFF438ZIY.tsv", header = T))
  
# gene ids here have a version number or something, remove this
dat4$gene_id <- unlist(lapply(strsplit(dat4$gene_id, "\\."), function(x) x[1]))
dat4 <- data.frame(dat4$gene_id, dat4$TPM)
colnames(dat4) <- c("gene_id", "ENCFF438ZIY_TPM")
results2 <- merge(dat4, results, by.x = "gene_id", by.y = "ensembl_gene_stable_id", all.x=T)
rm(dat4)
  
# drop any genes that dont pertain to the data
results2 <- results2[!is.na(results2$refsnp_id),]
  
# ENCSR229WIW 2016
# have to go back and look at this again, no idea what the values are here
# dat <- as.data.frame(fread("ENCFF361BBQ.tsv", header = T))
  
# ENCSR755LFM 2013
# transcript quantification (these 2 are identical)
dat <- as.data.frame(fread("ENCFF935XFP.tsv", header = T))
dat <- data.frame(dat$gene_id, dat$TPM)
colnames(dat) <- c("gene_id", "ENCFF935XFP_TPM")
dat$gene_id <- as.character(dat$gene_id)
dat$gene_id <- unlist(lapply(strsplit(dat$gene_id, "\\."), function(x) x[1]))
results3 <- merge(dat, results2, by.x = "gene_id", by.y = "gene_id", all.x = T)
rm(dat)
  
 # not sure how to add in TPPM data, tehre are multiple values per unique snp
results3 <- results3[!is.na(results3$refsnp_id),]

# never actually used this part of the data, leave it in for posterity
# ----------------------------------


# read in each annotation and attach to main crv data.frame

# --------- histone marks ----------
# ENCSR000EXA 2011
# take exp(-q) to get pval
# visualize?
# H3K4me1
print("adding histone marks...")
# originally wrote this to be able to attach any value from the bed file,
# rather than just a binary value; could use attachSnpAnnotSimple instead of
# attachSnpAnnot + convertToBinInt
#ENCSR000EXA <- readBed( "ENCFF450TSY.bed", narrowpeak.bed.header)
#crv.dat$ENCSR000EXA <- convertToBinInt( unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR000EXA, chr, "-log10(q)" )) )
#rm(ENCSR000EXA)


# ENCSR000EWZ 2011
# H3K9me3
#ENCSR000EWZ <- readBed( "ENCFF970JXR.bed", narrowpeak.bed.header)
#crv.dat$ENCSR000EWZ <- convertToBinInt( unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR000EWZ, chr, "-log10(q)" )) )
#rm(ENCSR000EWZ)

# ENCSR000EXC 2011
# H3K9ac
#ENCSR000EXC <- readBed( "ENCFF304WSW.bed", narrowpeak.bed.header)
#crv.dat$ENCSR000EXC <- convertToBinInt( unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR000EXC, chr, "-log10(q)" )) )
#rm(ENCSR000EXC)

# ENCSR494TNM 2016
# CTCF
#ENCSR494TNM <- readBed( "ENCFF146REQ.bed", narrowpeak.bed.header)
#crv.dat$ENCSR494TNM <- convertToBinInt( unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR494TNM, chr, "-log10(q)" )) )
#rm(ENCSR494TNM)

# ENCSR000EXB 2011
# H3K36me3
ENCSR000EXB <- readBed( "ENCFF082VKJ.bed", narrowpeak.bed.header)
crv.dat$ENCSR000EXB <- convertToBinInt( unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR000EXB, chr, "-log10(q)" )) )
rm(ENCSR000EXB)

# ENCSR611DJQ 2017
# H3K4me3
ENCSR611DJQ <- readBed("ENCFF047XWN.bed", narrowpeak.bed.header)
crv.dat$ENCSR611DJQ <- convertToBinInt( unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR611DJQ, chr, "-log10(q)")))
rm(ENCSR611DJQ)

# ENCSR000EXE 2011
# H3K27me3
ENCSR000EXE <- readBed("ENCFF355TTQ.bed", narrowpeak.bed.header)
crv.dat$ENCSR000EXE <- convertToBinInt( unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR000EXE, chr, "-log10(q)")))
rm(ENCSR000EXE)

# ENCSR091MNT 2019
# H3K27me3
ENCSR091MNT <- readBed("ENCFF327DDV.bed", narrowpeak.bed.header)
crv.dat$ENCSR091MNT <- convertToBinInt( unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR091MNT, chr, "-log10(q)")))
rm(ENCSR091MNT)

# ENCSR619EZG 2019
# H3K4me3
ENCSR619EZG <- readBed("ENCFF305TNC.bed", narrowpeak.bed.header)
crv.dat$ENCSR619EZG <- convertToBinInt( unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR619EZG, chr, "-log10(q)")))
rm(ENCSR619EZG)

# ENCSR956VQB 2019
# H3K4me1
ENCSR956VQB <- readBed("ENCFF620AJW.bed", narrowpeak.bed.header)
crv.dat$ENCSR956VQB <- convertToBinInt( unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR956VQB, chr, "-log10(q)")))
rm(ENCSR956VQB)

# ENCSR136ZQZ 2019
# H3K27ac
ENCSR136ZQZ <- readBed("ENCFF567XUE.bed", narrowpeak.bed.header)
crv.dat$ENCSR136ZQZ <- convertToBinInt( unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR136ZQZ, chr, "-log10(q)")))
rm(ENCSR136ZQZ)

# ENCSR954IGQ 2019
# H3K27ac
ENCSR954IGQ <- readBed("ENCFF100QOV.bed", narrowpeak.bed.header)
crv.dat$ENCSR954IGQ <- convertToBinInt( unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR954IGQ, chr, "-log10(q)")))
rm(ENCSR954IGQ)

# ENCSR376JOC
# H3K9me3
ENCSR376JOC <- readBed("ENCFF418NCI.bed", narrowpeak.bed.header)
crv.dat$ENCSR376JOC <- convertToBinInt( unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR376JOC, chr, "-log10(q)")))
rm(ENCSR376JOC)
print("done!")
# ----------------------------------


# ------ open chromatin marks ------


# H3K4me3
#print("adding open chromatin marks...")


# ENCSR303XKE 2017
# 5-group for testis male fetal
ENCSR303XKE <- readBed("ENCFF218UBN.bed", cis.reg.header)
crv.dat$ENCSR303XKE <- convertToBinInt(unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR303XKE, chr, "Cis.Reg.Name" )))
rm(ENCSR303XKE)


# ENCSR000EPS 2011 (exp and replicate)
# UW human NT2-D1 DNase-seq
# signal value, but no p or q val
dat5 <- readBed("ENCFF366XFZ.bed", narrowpeak.bed.header)
dat6 <- readBed("ENCFF506YRM.bed", narrowpeak.bed.header)
ENCSR000EPS <- joinMethylReplicates(dat5, dat6, chr)
crv.dat$ENCSR000EPS <- convertToBinInt(unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR000EPS, chr, "q" )))

# write.table( ENCSR000EPS[,c(2,1,3,4)], "ENCSR000EPS.bed", col.names = F, row.names = F, quote = F)
rm(ENCSR000EPS)
rm(dat5)
rm(dat6)

# ENCSR729DRB 2013
# male embryo testis tissue
ENCSR729DRB <- readBed("ENCFF843ZSC.bed", narrowpeak.bed.header)
crv.dat$ENCSR729DRB <- convertToBinInt(unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR729DRB, chr, '-log10(q)' )))
rm(ENCSR729DRB)

# ENCSR278FHC 2013
# narrowpeak
# male embryo testis tissue
ENCFF012QTD <- readBed("ENCFF012QTD.bed", narrowpeak.bed.header)
crv.dat$ENCFF012QTD <- convertToBinInt(unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCFF012QTD, chr, '-log10(q)' )))
rm(ENCFF012QTD)

# broadpeak - differnet header
ENCFF841TKB <- readBed("ENCFF841TKB.bed", broadpeak.bed.header)
crv.dat$ENCFF841TKB <- convertToBinInt(unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCFF841TKB, chr, 'q' )))
rm(ENCFF841TKB)
print("done!")
# ----------------------------------


# ----- methylation ------ #
print("adding methylation...")
# ENCSR000DED 2011
# RRBS on testis (w/ replicate)
dat1 <- readBed("ENCFF001TKY.bed", methyl.bed.header)
dat2 <- readBed("ENCFF001TKZ.bed", methyl.bed.header)
ENCSR000DED <- joinMethylReplicates( dat1, dat2 )
rm(dat1)
rm(dat2)
write.table( convertToBED(ENCSR000DED), "methyl1.bed", col.names = TRUE, row.names = FALSE, quote = FALSE)

crv.dat$ENCSR000DED <- convertToBinInt(unlist(lapply( crv.dat$h19pos, attachSnpAnnotMethyl, ENCSR000DED, chr, "percentMeth" )))
rm(ENCSR000DED)




# ENCSR080YRO
# RRBS on testis (w/ replicate)
# ---
dat1 <- readBed("ENCFF001TPW.bed", methyl.bed.header)
dat2 <- readBed("ENCFF001TPX.bed", methyl.bed.header)

ENCSR080YRO <- joinMethylReplicates( dat1, dat2, chr )

write.table( convertToBED(ENCSR080YRO), "ENCSR080YRO.bed", col.names = TRUE, row.names = FALSE, quote = FALSE)

crv.dat$ENCSR080YRO <- convertToBinInt(unlist(lapply( crv.dat$h19pos, attachSnpAnnotMethyl, ENCSR080YRO, chr, "percentMeth" )))
rm(ENCSR080YRO)

methyl.bed.header.short <- methyl.bed.header[c(1,2,3,5,6,11)]
# ENCSR011HZJ 2017 - these are in grch38!!
# CHG meth
ENCFF507JBR <- readBed(paste("ENCFF507JBR/", chr, ".bed", sep = ""), methyl.bed.header.short)
crv.dat$ENCFF507JBR <- convertToBinInt(unlist(lapply( crv.dat$h38pos, attachSnpAnnotMethyl, ENCFF507JBR, chr, "percentMeth" )))
rm(ENCFF507JBR)

# CHH meth
ENCFF038JFQ <- readBed(paste("ENCFF038JFQ/", chr, ".bed", sep = ""), methyl.bed.header.short[c(1,2,3,4,6)])
crv.dat$ENCFF038JFQ <- convertToBinInt(unlist(lapply( crv.dat$h38pos, attachSnpAnnotMethyl, ENCFF038JFQ, chr, "percentMeth" )))
rm(ENCFF038JFQ)

# CpG meth
ENCFF715DMX <- readBed(paste("ENCFF715DMX/", chr, ".bed", sep = ""), methyl.bed.header.short[c(1,2,3,4,6)])
crv.dat$ENCFF715DMX <- convertToBinInt(unlist(lapply( crv.dat$h38pos, attachSnpAnnotMethyl, ENCFF715DMX, chr, "percentMeth" )))
rm(ENCFF715DMX)

# ENCSR806NNG 2018
# idats -----
# the pair of idat files needs to have the same basename, and be appened with "_Red" and "_Grn"
# eg, ENCFF001RHW.idat becomes ENCSR000ABD_Red.idat
# and ENCFF001RHX.idat becomes ENCSR000ABD_Grn.idat
# cell.line, adult male testis
# Whole genome seq
ENCSR806NNG <- readMethylIdat("ENCSR000ABD")
write.table( convertToBED(ENCSR806NNG), "ENCSR806NNG", col.names = TRUE, row.names = FALSE, quote = FALSE)
crv.dat$ENCSR806NNG <- convertToBinInt(unlist(lapply( crv.dat$h19pos, attachSnpAnnotMethyl, ENCSR806NNG, chr, "name" )))
rm(ENCSR806NNG)


# ENCSR304AIL 
# tissue, adult male testis
# DNAme
ENCSR304AIL <- readMethylIdat( "ENCSR304AIL" )
write.table( convertToBED(ENCSR304AIL), "ENCSR304AIL", col.names = TRUE, row.names = FALSE, quote = FALSE)
crv.dat$ENCSR304AIL <- convertToBinInt(unlist(lapply( crv.dat$h19pos, attachSnpAnnotMethyl, ENCSR304AIL, chr, "name" )))
rm(ENCSR304AIL)

# ENCSR304PMI
# RRBS cell line
ENCSR304PMI <- readMethylIdat("ENCSR304PMI")
write.table( convertToBED(ENCSR304PMI), "ENCSR304PMI.bed", col.names = TRUE, row.names = FALSE, quote = FALSE)
crv.dat$ENCSR304PMI <- convertToBinInt(unlist(lapply( crv.dat$h19pos, attachSnpAnnotMethyl, ENCSR304PMI, chr, "name" )))
rm(ENCSR304PMI)

# ENCSR962XHD 2017
# DNAme
# cell line
ENCSR962XHD <- readMethylIdat("ENCLB381OUG")
write.table( convertToBED(ENCSR962XHD), "ENCSR962XHD.bed", col.names = TRUE, row.names = FALSE, quote = FALSE)
crv.dat$ENCSR962XHD <- convertToBinInt(unlist(lapply( crv.dat$h19pos, attachSnpAnnotMethyl, ENCSR962XHD, chr, "name" )))
rm(ENCSR962XHD)

# ENCSR942OLI 2017
# tissue adult male testis
# DNAme
ENCSR942OLI <- readMethylIdat("ENCSR942OLI")
write.table( convertToBED(ENCSR942OLI), "ENCSR942OLI.bed", col.names = TRUE, row.names = FALSE, quote = FALSE)
crv.dat$ENCSR942OLI <- convertToBinInt(unlist(lapply( crv.dat$h19pos, attachSnpAnnotMethyl, ENCSR942OLI, chr, "name" )))
rm(ENCSR942OLI)

print("done!")
# ----------------------------------

# ---------- transcription start sites --------- #
print("adding transcription start sites...")
# what are the values here
# adult male testis RAMPAGE
# gene quantifications, or transcription start sites?
# ENCSR866SRG 2016
#dat <- read.table("ENCFF874RBV.tsv", header = FALSE)
#tss.header <- c("chr", "Start", "End", "Name1", "Score", "Strand", "idk", "gene1", "Gene", "gene2")
# colnames(dat) <- tss.header
tss.header <- c("Chrom", "Start", "End", "Name1", "const", "Strand", "Score", "Name2", "Name3", "Name4", "coords")
ENCSR866SRG = readBed("ENCFF648HUU.bed", tss.header)
crv.dat$ENCSR866SRG <- convertToBinInt(unlist(lapply(crv.dat$h19pos, attachSnpAnnot, ENCSR866SRG, chr, "Score")))
rm(ENCSR866SRG)

# ENCSR841EQJ
# rampage of testis
# columns 4, 8, 9, 10 all look like they have identical information
ENCSR841EQJ <- readBed("ENCFF162KMZ.bed", tss.header)
crv.dat$ENCSR841EQJ <- convertToBinInt(unlist(lapply(crv.dat$h19pos, attachSnpAnnot, ENCSR841EQJ,chr, "Score")))
rm(ENCSR841EQJ)
print("done!")
# ---------------------------------------------- #

# ------- transcription factor binding site ---- #
print("adding trnascription factor binding sites...")
# ENCSR000EWY 2011
# ZNF274
ENCSR000EWY <- readBed("ENCFF637LOO.bed", narrowpeak.bed.header)
crv.dat$ENCSR000EWY <- convertToBinInt(unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR000EWY, chr, "q" )))
rm(ENCSR000EWY)


# ENCSR000EXG 2011
# YY1
#hg38
ENCSR000EXG <- readBed("ENCFF293PVG.bed", narrowpeak.bed.header)
crv.dat$ENCSR000EXG <- convertToBinInt(unlist(lapply( crv.dat$h38pos, attachSnpAnnot, ENCSR000EXG, chr, "q" )))
rm(ENCSR000EXG)

# ENCSR981CID 2017
# CTCF adult male 54
ENCSR981CID <- readBed("ENCFF885KKQ.bed", narrowpeak.bed.header)
crv.dat$ENCSR981CID <- convertToBinInt(unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR981CID, chr, "q" )))
rm(ENCSR981CID)

# ENCSR803FAP 2017
# POLR2A
ENCSR803FAP <- readBed("ENCFF940TNN.bed", narrowpeak.bed.header)
crv.dat$ENCSR803FAP <- unlist(lapply( crv.dat$h19pos, attachSnpAnnot, ENCSR803FAP, chr, "q"))
rm(ENCSR803FAP)
print("done!")
# --------------------------------------------- #


# local data from straun; use files with ld > .4
struan.dat1 <- read.table("struan/Testis_Cancer_novel_NTERA2_0.4_snpInOe_ann.csv", sep = ",", header = TRUE)
#struan.dat2 <- read.table("struan/Testis_Cancer_novel_NTERA2_snpInOe_ann_new.csv", sep = ",", header = TRUE)
struan.dat1.rep <- read.table("struan/Testis_Cancer_replicate_NTERA2_0.4_snpInOe_ann.csv", sep = ",", header = TRUE)

struan.dat <- rbind(struan.dat1, struan.dat1.rep)
struan.dat$INT_FRAG <- as.character(struan.dat$INT_FRAG)


tmp <- unlist(lapply(strsplit(struan.dat$INT_FRAG, ":"), function(x) x[2]))
struan <- data.frame( Chrom <- unlist(lapply(strsplit(struan.dat$INT_FRAG, ":"), function(x) x[1])),
                      Start <- unlist(lapply(strsplit(tmp, "-"), function(x) x[1])),
                      End <- unlist(lapply(strsplit(tmp, "-"), function(x) x[2])))
colnames(struan) <- c("Chrom", "Start", "End")

# this line is a duplicate
struan <- struan[-17,]
#write.table(struan, "struan.bed", col.names = F, row.names = F, quote = F)

crv.dat$HIC <- 0
crv.dat$HIC[snpAnnotSimple(crv.dat$h19pos, struan, chr)] <- 1

# local cell line data
EP2102 <- read.table("forChey/2102EP.bed", header = F)
crv.dat$EP2102 <- 0
crv.dat$EP2102[snpAnnotSimple(crv.dat$h19pos, EP2102, chr)] <- 1

TCAM2 <- read.table("forChey/TCAM2.bed", header = F)
crv.dat$TCAM2 <- 0
crv.dat$TCAM2[snpAnnotSimple(crv.dat$h19pos, TCAM2, chr)] <- 1

NTERA2 <- read.table("forChey/NTERA2.bed", header = F)
crv.dat$NTERA2 <- 0
crv.dat$NTERA2[snpAnnotSimple(crv.dat$h19pos, NTERA2, chr)] <- 1

NCCIT <- read.table("forChey/NCCIT.bed", header = F)
crv.dat$NCCIT <- 0
crv.dat$NCCIT[snpAnnotSimple(crv.dat$h19pos, NCCIT, chr)] <- 1
# ----- end annotations


# qc
k <- dim(crv.dat)[2]

print("annotation qc- remove empty columns and NA values...")
crv.dat[,11:k] <- apply(crv.dat[,11:k], 2, replaceNA)

is.emptyAnnot <- function( x )
{
  return(all(x == 0))
}



# ---

# remove any annotations that had no matches in the credset
ind <-  apply(crv.dat, 2, is.emptyAnnot)
crv.dat.bak <- crv.dat
crv.dat <- crv.dat[,!ind]
print("done!")

print(paste("crv.dat has ", dim(crv.dat)[2], " columns", sep = ""))
if(dim(crv.dat)[2] == 10)
{
  print("no valid annotations! exiting")
  stop()
}


# write the output for PAINTOR
print("writing output...")

ANNOTFILE <- paste(INFILENAME, "/", INFILENAME, ".annotations", sep = "")
annot <- crv.dat[,11:dim(crv.dat)[2]]
write.table(annot, ANNOTFILE, 
            col.names = T, row.names = F, quote = F)


# note number of snps and number of redundant snp
annotation.names <- paste(colnames(annot), collapse = ",")
ANNOTNAMESFILE <- paste(INFILENAME, "annotation.names", sep = "/")
write.table(annotation.names, ANNOTNAMESFILE, quote = F, col.names = F, row.names = F, append = F)


INPUTFILES <- paste(INFILENAME, "input.files", sep = "/")
write.table(INFILENAME, INPUTFILES, quote = F, col.names = F, row.names = F, append = F)
print("done!")
print("successfully completed")


xrange <- c(32008931, 320083111)
bed <- readBed( "ENCFF082VKJ.bed", narrowpeak.bed.header)
#p1 <- ggplot(data = bed, aes(x = Start)) + geom_area(stat = "bin", binwidth = 50)
#p1

bed <- bed[ !is.na(bed$Chrom), ]
bed <- bed[ as.character(bed$Chrom) == "chr6" & bed$Start >= xrange[1] & bed$End <= xrange[2],]

# TNXB: 32008931 - 320083111
p1 <- ggplot(data = bed, aes( x = Start, y = -log10(p))) + geom_area(fill = "blue")
             p1


p1 <- ggplot(data = hb) + geom_rect(aes(xmin = BP1, xmax = BP2, ymin = 0.975, ymax = 1.025), fill = "black") +
  geom_hline(yintercept = 1, size = 2)
p1

# -- find eqtls
eqtl <- read.table("significantresults_TECAC.GTEx.csv", sep = ",", header = TRUE)
eqtl$Chr <- as.integer(gsub("chr", "", eqtl$Chr))
#eqtl <- eqtl[which(eqtl$Tissue == "Testis"),]
CHR <- unique(crv.dat$chr)
eqtl <- eqtl[eqtl$Chr == CHR,]

novel.genes <- c("PPP2R5A",
                 "BCL2L11",
                 "TERT",
                   "TNXB",
                 "BAK1",
                   "DEPTOR",
                 "DMRT1",
                   "PSMB7",
                 "ANAPC2",
                 "ITIH5",
                 "ARL14EP",
                 "RAD52",
                 "METTL7A",
                 "SP1",
                 "CYTH1",
                 "ENOSF1",
                 "CTC459",
                 "ZNF217",
                 "FAM48B1",
                 "AR",
                 "CENPI",
                 "TKTL1")

tmp.gene <- c("DEPTOR", "PSMB7", 
              "ITIH5", "METTL7A", 
              "CTC459", "ZNF217", 
              "FAM48B1")
tmp.ensg <- c("ENSG00000155792", "ENSG00000136930",
              "ENSG00000123243", "ENSG00000185432",
              "ENSG00000267575", "ENSG0000171940",
              "ENSG00000223731")
nov.ensg <- gene.key[ gene.key$gene.id %in% novel.genes,]

ensg <- c(tmp.ensg, as.character(nov.ensg$ensembl))

eqtl2 <- eqtl[ eqtl$GeneName %in% c(tmp.gene, as.character(nov.ensg$gene.id)),]

genes <- unique(c(tmp.gene, as.character(nov.ensg$gene.id)))
fuma <- read.table("~/Desktop/gtex_v8_ts_avg_log2TPM_exp.txt", header = TRUE)
fuma.sml <- fuma[ , colnames(fuma) %in% c("ensg", "symbol", "Testis")]
fuma.sml <- fuma.sml[ fuma.sml$ensg %in% ensg,]

# Beta_GWAS is based off of the FIRST snp in ALlele_GWAS
# Beta_eQTL is based off of the SECOND snp in allele_eqtl
# ex: chr8:119914393, BetaGWAS = -0.1259, Allele_GWAS = g_c,
# BetaEqtl = -0.276993, Allele_eqtl = C_G
# these are both based off the C allele

# BCL2L11
# betaGWAS -0.1147, alleGWAS c_t
# betaEQTL -0.129164 alleleEQTL C_T
# eqtls from shweta are in hg38
ind = inEqtl( crv.dat$h38pos[32], eqtl)
sum( !is.na( unique(eqtl$Tissue[ind])))

c("Testis", "testis") %in% eqtl$Tissue[ind]


crv.dat <- crv.dat.bak
k = dim(crv.dat)[2]
# plot annotation matrix
mat <- as.matrix(crv.dat[,11:k])
colnames(mat) <- colnames(crv.dat)[11:k]
rownames(mat) <- crv.dat$SNP.ID.hg19

annot <- read.table("../annotation-table.csv", header = T, sep = ",")
annot$Annotation %in% colnames(mat)

# drop utr3', utr5', exons
k = dim(mat)[2]
mat <- mat[,4:k]

# aligns annotations to whats in the text database
# this is used for ordering the annotations eg by type
ind = match(as.character(annot$Annotation), colnames(mat))
mat <- mat[,ind]



u.stat <- read.table(paste(getwd(), INFILENAME, "univariate.stats", sep = "/"), header = F)
u.stat <-  as.character(u.stat$V2)
u.stat <- data.frame( annot = unlist(lapply(strsplit(u.stat, ","), function(x) return(x[1]))),
                      p = unlist(lapply(strsplit(u.stat, ","), function(x) return(x[2]))) )
u.stat$p <- as.numeric(as.character(u.stat$p))


mat <- melt(mat)
ind <- match(mat$Var2, u.stat$annot  ) 
ind2 <- mat$value != 0
mat$value[!is.na(ind) ] <- u.stat$p[ ind[!is.na(ind)]]
mat$value[!ind2] <- 0
mat$log.value <- -log10(mat$value)
mat$log.value[ mat$log.value == Inf ] <- NA


mat$Var2 <- as.character(mat$Var2)
write.table(mat, paste0(INFILENAME, "_annot_plot.txt"), col.names = TRUE, quote = F, row.names = F)


library(RColorBrewer)
mat$log.value[ mat$log.value == 0 ] <- NA
p1 <- ggplot(data = mat, aes(x = Var2, y = Var1, fill=log.value)) + geom_tile() +
  geom_tile(data = mat[ mat$log.value != 0,], color = "black") +
  xlab("Annotation") +
  ylab("CRV") + 
  theme_minimal() +
  
  theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5), 
        #     legend.position = "none",
        panel.border = element_rect(color = "black", fill = NA)) +
 # ggtitle(gsub("_", ":", INFILENAME)) +
 scale_fill_gradient2( low = "red", mid = "white", high = "green",
                       midpoint = -log10(0.05),
                       breaks = c(-1, -log10(0.05), -log10(0.001), -log10(0.0001)),
                       labels = c(1, 0.05, 0.001, 0.0001),
                       name = "p-value",
                       limits = c(-1,4),
                       na.value = "white") 
  p1
  
 #  scale_fill_gradientn( colors = brewer.pal(n = 9, name = "Blues"), 
#                        breaks = c(0, -log10(0.05), -log10(0.001), -log10(0.00001)),
 #                       labels = c(1, 0.05, 0.001, 0.00001),
  #                      name = "p-value",
   #                     limits = c(0,5),
    #                    na.value = "white") 
#p1
# scale_fill_gradient2( guide = "colourbar", 
#                     low = "white",
#                    high = "royalblue3",
#                   midpoint = -log10(1),
#                 limits = c(0,4),
#                breaks = c(0, -log10(0.05), -log10(0.001), -log10(0.0001)),
#               labels = c(1, 0.05, 0.001, 0.0001),
#              name = "p-value") 
png(paste0(INFILENAME, "_annot.png"), width = 600, height = dim(crv.dat)[1] * 15)
print(p1)
dev.off()


# PAINTOR debugging-
# make sure tabs instead of spaces
# ld matrix is invertable
# same number of snps in all files

# run paintor as follows:
#~/PAINTOR_V3.0/PAINTOR -input input.files -in . -out . 
# -Zhead z -LDname ld -Gname MVTEST -Lname BF_MVTEST -enumerate 4 -annotations H3K36me3,H3K4me3.cisreg
# use getLRT.p to get the pvalue comparing annotation model to baseline


dat <- read.table("hic_snps.txt",header = F )

# empty cells
dat <- dat[-549,]
dat$bp <- unlist(lapply(strsplit(dat$V1, "--"), function(x) x[1]))
dat$refbp <- unlist(lapply(strsplit(dat$V1, "_"), function(x) x[2]))

dat$chr <- unlist(lapply(strsplit(dat$V1, "--"), function(x) x[2]))
dat$chr <- paste0("chr", unlist(lapply(strsplit(dat$chr, "_"), function(x) x[1])))


struan.dat <- rbind(struan.dat1, struan.dat1.rep)
struan.dat$INT_FRAG <- as.character(struan.dat$INT_FRAG)


tmp <- unlist(lapply(strsplit(struan.dat$INT_FRAG, ":"), function(x) x[2]))
struan <- data.frame( Chrom <- unlist(lapply(strsplit(struan.dat$INT_FRAG, ":"), function(x) x[1])),
                      Start <- unlist(lapply(strsplit(tmp, "-"), function(x) x[1])),
                      End <- unlist(lapply(strsplit(tmp, "-"), function(x) x[2])))
colnames(struan) <- c("Chrom", "Start", "End")

# this line is a duplicate
struan <- struan[-17,]

#write.table(struan, "struan.bed", col.names = F, row.names = F, quote = F)
dat$bp <- as.integer(dat$bp)

getOverlap <- function( dat1, dat2)
{
  ind <- rep(0, dim(dat1)[1])
  for( i in 1:dim(dat1)[1])
  {
    if( length( snpAnnotSimple(dat1$bp[i], dat2, dat1$chr[i]) ) > 0 )
    {
      ind[i] <- 1
    }
  }
  return(ind)
}

dat$struan <- getOverlap(dat, struan)
dat$EP2102 <- getOverlap(dat, EP2102)
dat$TCAM2 <- getOverlap(dat, TCAM2)
dat$NNCIT <- getOverlap(dat, NCCIT)
dat$NTERA2 <- getOverlap(dat, NTERA2)

write.table(dat, "hic_ATAC_seq.txt", row.names = F, col.names = T, quote = F)






