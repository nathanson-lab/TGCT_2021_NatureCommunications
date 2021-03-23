# pluta 11/23/20


# functions snp.to.rsid and motifbreak do not work without modification, they are recoded
# as snp.to.rsid2 and motifbreak2.
# this introduces a problem with inheritance, all libraries and functions need to be
# loaded locally. most of this code is directly from motifbreakR
rm(list=ls())
library(motifbreakR)
library(TFMPvalue)
library(matrixStats)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BiocParallel)
library(SNPlocs.Hsapiens.dbSNP142.GRCh37)
library(MotifDb)
library(GenomeInfoDb)
library(Gviz)
library(motifStack)
library(stringr)
data(motifbreakR_motif)
data(hocomoco)

# ghostwriter executable
Sys.setenv(R_GSCMD="/usr/bin/local/gs")



# ======================================================================= #
# ======================================================================= #
# local copies of motifbreakR functions
# ----------------------------------------------------------------------- #
plotMotifLogoStack.2 <- function(pfms, ...) {
  pfms <- rev(pfms)
  n <- length(pfms)
  lapply(pfms, function(.ele) {
    if (!is(.ele, "pfm"))
      stop("pfms must be a list of class pfm")
  })
  opar <- par(mfrow = c(n, 1), mar = c(3.5, 3.5, 1.5, 0.5))
  assign("tmp_motifStack_symbolsCache", list(), pos = ".GlobalEnv")
  motifStack::plotMotifLogo(pfms[[1]], motifName = pfms[[1]]@name, p = rep(0.25, 4))
  for (i in seq.int(n)[-1]) {
    motifStack::plotMotifLogo(pfms[[n]], motifName = pfms[[n]]@name,
                              p=rep(0.25, 4), xlab = NA, newpage = FALSE)
  }
  rm(list = "tmp_motifStack_symbolsCache", pos = ".GlobalEnv")
  par(opar)
}
# ----------------------------------------------------------------------- #


# ----------------------------------------------------------------------- #
#' @importFrom grid grid.newpage pushViewport viewport popViewport
plotMotifLogoStack.3 <- function(pfms, ...) {
  n <- length(pfms)
  lapply(pfms, function(.ele) {
    if (class(.ele) != "pfm")
      stop("pfms must be a list of class pfm")
  })
  assign("tmp_motifStack_symbolsCache", list(), pos = ".GlobalEnv")
  # grid.newpage()
  ht <- 1/n
  y0 <- 0.5 * ht
  for (i in rev(seq.int(n))) {
    pushViewport(viewport(y = y0, height = ht))
    plotMotifLogo(pfms[[i]], motifName = pfms[[i]]@name, ncex = 1,
                  p = pfms[[i]]@background, colset = pfms[[i]]@color,
                  xlab = NA, newpage = FALSE, margins = c(1.5, 4.1,
                                                          1.1, 0.1), ...)
    popViewport()
    y0 <- y0 + ht
  }
  rm(list = "tmp_motifStack_symbolsCache", pos = ".GlobalEnv")
  return()
}
# ----------------------------------------------------------------------- #

# ----------------------------------------------------------------------- #
prepareVariants <- function(fsnplist, genome.bsgenome, max.pwm.width, legacy) {
  k <- max.pwm.width; rm(max.pwm.width)
  ref_len <- nchar(fsnplist$REF)
  alt_len <- nchar(fsnplist$ALT)
  is.indel <- ref_len > 1L | alt_len > 1L
  ## check that reference matches ref genome
  equals.ref <- getSeq(genome.bsgenome, fsnplist) == fsnplist$REF
  if (!all(equals.ref)) {
    stop(paste(names(fsnplist[!equals.ref]), "reference allele does not match value in reference genome",
               sep = " "))
  }
  if (sum(is.indel) < length(is.indel) & sum(is.indel) > 0L) {
    if (legacy) {
      warning("Indels are included in variant input set, but legacy scoring was enabled.\n",
              "Legacy scoring is not availble with indels and they will be dropped from analysis")
      fsnplist.indel <- NULL
      fsnplist.snv <- fsnplist[!is.indel]
    } else {
      fsnplist.indel <- fsnplist[is.indel]
      fsnplist.snv <- fsnplist[!is.indel]
    }
  } else if (sum(is.indel) == length(is.indel)) {
    fsnplist.indel <- fsnplist
    fsnplist.snv <- NULL
    if (legacy) {
      warning("The only variants included in the input set are indels, but legacy scoring was selected.\n",
              "Legacy scoring is not availble for use with indels and will be disabled.")
      legacy <- FALSE
    }
  } else if (sum(is.indel) == 0L) {
    fsnplist.indel <- NULL
    fsnplist.snv <- fsnplist
    if (!legacy) {
      warning("The only variants included in the input set are SNVs, and legacy scoring was not selected.\n",
              "A new scoring algorithm will be used, and may present different scores than previously run.")
      legacy <- FALSE
    }
  }
  if (!is.null(fsnplist.indel)) {
    snp.sequence.ref.indel <- getSeq(genome.bsgenome, promoters(fsnplist.indel, upstream = k - 1,
                                                                downstream = k + max(nchar(fsnplist.indel$REF))))
    at <- as(IRanges(start = k, width = width(fsnplist.indel)), "IRangesList")
    snp.sequence.alt.indel <- DNAStringSet(Map(replaceAt,
                                               x = snp.sequence.ref.indel,
                                               at = at,
                                               fsnplist.indel$ALT))
    
    need.alignment <- !(lengths(fsnplist.indel$REF) == 1 | lengths(fsnplist.indel$ALT) == 1)
    insertion.var <- lengths(fsnplist.indel$REF) < lengths(fsnplist.indel$ALT)
    fsnplist.indel$ALT_loc <- 1L
    fsnplist.indel[insertion.var]$ALT_loc <- Map(seq,
                                                 from = nchar(fsnplist.indel[insertion.var]$REF) + 1L,
                                                 to = nchar(fsnplist.indel[insertion.var]$ALT))
    fsnplist.indel[!insertion.var]$ALT_loc <- Map(seq,
                                                  from = nchar(fsnplist.indel[!insertion.var]$ALT) + 1L,
                                                  to = nchar(fsnplist.indel[!insertion.var]$REF))
    ref.len <- nchar(fsnplist.indel[nchar(fsnplist.indel$REF) == nchar(fsnplist.indel$ALT)]$REF)
    fsnplist.indel[nchar(fsnplist.indel$REF) == nchar(fsnplist.indel$ALT)]$ALT_loc <- lapply(ref.len,
                                                                                             function(x) {
                                                                                               seq(from = 1L,
                                                                                                   to = x)
                                                                                             })
    fsnplist.indel$varType <- "Other"
    if (sum(insertion.var) > 0) {
      fsnplist.indel[insertion.var, ]$varType <- "Insertion"
    }
    if (sum(lengths(fsnplist.indel$REF) > lengths(fsnplist.indel$ALT)) > 0) {
      fsnplist.indel[lengths(fsnplist.indel$REF) > lengths(fsnplist.indel$ALT)]$varType <- "Deletion"
    }
    if (any(need.alignment)) {
      need.del <- any(!insertion.var & need.alignment)
      need.ins <- any(insertion.var & need.alignment)
      if (need.del) {
        pattern.del <- Map(matchPattern,
                           fsnplist.indel[!insertion.var & need.alignment]$ALT,
                           fsnplist.indel[!insertion.var & need.alignment]$REF,
                           with.indels = F, max.mismatch = 0)
        names(pattern.del) <- names(fsnplist.indel[!insertion.var & need.alignment])
        pattern.del <- sapply(pattern.del, function(x) {
          slen <- length(subject(x))
          x <- x[start(x) == 1 | end(x) == slen]
          if (length(x) > 0) {
            x <- x[1]
            x <- gaps(x)
            x <- as(x, "IRanges")
          }
        })
        pattern.del.valid <- sapply(pattern.del, function(x) {length(x) > 0})
        nr.pattern.del <- lapply(pattern.del[pattern.del.valid], function(x) {start(x):end(x)})
        fsnplist.indel[!insertion.var & need.alignment][names(pattern.del[pattern.del.valid])]$ALT_loc <- nr.pattern.del
        if (any((!insertion.var & need.alignment)[!pattern.del.valid])) {
          alignment.del <- Map(pairwiseAlignment,
                               fsnplist.indel[!insertion.var & need.alignment][!pattern.del.valid]$ALT,
                               fsnplist.indel[!insertion.var & need.alignment][!pattern.del.valid]$REF,
                               type = "global")
          names(alignment.del) <- names(fsnplist.indel[!insertion.var & need.alignment][!pattern.del.valid])
          alignment.del <- sapply(sapply(alignment.del, deletion), unlist)
          alignment.del.valid <- sapply(alignment.del, function(x) {length(x) > 0})
          nr.alignment.del <- lapply(alignment.del[alignment.del.valid], function(x) {start(x):end(x)})
          fsnplist.indel[!insertion.var & need.alignment][names(alignment.del[alignment.del.valid])]$ALT_loc <- nr.alignment.del
          rm(alignment.del, nr.alignment.del)
        }
        rm(pattern.del, nr.pattern.del, pattern.del.valid, need.del)
      }
      if (need.ins) {
        pattern.ins <- Map(matchPattern,
                           fsnplist.indel[insertion.var & need.alignment]$REF,
                           fsnplist.indel[insertion.var & need.alignment]$ALT,
                           with.indels = F, max.mismatch = 0)
        names(pattern.ins) <- names(fsnplist.indel[insertion.var & need.alignment])
        pattern.ins <- sapply(pattern.ins, function(x) {
          slen <- length(subject(x))
          x <- x[start(x) == 1 | end(x) == slen]
          if (length(x) > 0) {
            x <- x[1]
            x <- gaps(x)
            x <- as(x, "IRanges")
          }
        })
        pattern.ins.valid <- sapply(pattern.ins, function(x) {length(x) > 0})
        nr.pattern.ins <- lapply(pattern.ins[pattern.ins.valid], function(x) {start(x):end(x)})
        fsnplist.indel[insertion.var & need.alignment][names(pattern.ins[pattern.ins.valid])]$ALT_loc <- nr.pattern.ins
        if (any((insertion.var & need.alignment)[!pattern.ins.valid])) {
          alignment.ins <- Map(pairwiseAlignment,
                               fsnplist.indel[insertion.var & need.alignment][!pattern.ins.valid]$ALT,
                               fsnplist.indel[insertion.var & need.alignment][!pattern.ins.valid]$REF,
                               type = "global")
          names(alignment.ins) <- names(fsnplist.indel[insertion.var & need.alignment][!pattern.ins.valid])
          alignment.ins <- sapply(sapply(alignment.ins, insertion), unlist)
          alignment.ins.valid <- sapply(alignment.ins, function(x) {length(x) > 0})
          nr.alignment.ins <- lapply(alignment.ins[alignment.ins.valid], function(x) {start(x):end(x)})
          fsnplist.indel[insertion.var & need.alignment][names(alignment.ins[alignment.ins.valid])]$ALT_loc <- nr.alignment.ins
          rm(alignment.ins, nr.alignment.ins)
        }
        rm(pattern.ins, nr.pattern.ins, pattern.ins.valid, need.ins)
      }
      rm(ref.len, insertion.var, need.alignment)
    }
  }
  if (!is.null(fsnplist.snv)) {
    snp.sequence.ref.snv <- getSeq(genome.bsgenome, promoters(fsnplist.snv, upstream = k - 1,
                                                              downstream = k + 1))
    at <- matrix(FALSE, nrow = length(snp.sequence.ref.snv), ncol = (k * 2))
    at[, k] <- TRUE
    snp.sequence.alt.snv <- replaceLetterAt(snp.sequence.ref.snv, at, fsnplist.snv$ALT)
    fsnplist.snv$ALT_loc <- 1L
    fsnplist.snv$varType <- "SNV"
  }
  if (sum(is.indel) < length(is.indel) & sum(is.indel) > 0L) {
    fsnplist <- c(fsnplist.indel, fsnplist.snv)
    snp.sequence.alt <- strsplit(as.character(c(snp.sequence.alt.indel,
                                                snp.sequence.alt.snv)), "")
    snp.sequence.ref <- strsplit(as.character(c(snp.sequence.ref.indel,
                                                snp.sequence.ref.snv)), "")
    rm(fsnplist.indel, fsnplist.snv,
       snp.sequence.alt.indel, snp.sequence.ref.indel,
       snp.sequence.alt.snv, snp.sequence.ref.snv)
  } else if (sum(is.indel) == length(is.indel)) {
    fsnplist <- fsnplist.indel
    snp.sequence.alt <- strsplit(as.character(snp.sequence.alt.indel), "")
    snp.sequence.ref <- strsplit(as.character(snp.sequence.ref.indel), "")
    rm(fsnplist.indel, snp.sequence.alt.indel, snp.sequence.ref.indel)
  } else if (sum(is.indel) == 0L) {
    fsnplist <- fsnplist.snv
    snp.sequence.alt <- strsplit(as.character(snp.sequence.alt.snv), "")
    snp.sequence.ref <- strsplit(as.character(snp.sequence.ref.snv), "")
    rm(fsnplist.snv, snp.sequence.alt.snv, snp.sequence.ref.snv)
  }
  rm(at); gc()
  return(list(fsnplist = fsnplist,
              ref.seq = snp.sequence.ref,
              alt.seq = snp.sequence.alt))
}
# ----------------------------------------------------------------------- #



# ----------------------------------------------------------------------- #
varEff <- function(allelR, allelA) {
  score <- allelA - allelR
  if (abs(score) >= 0.7) {
    return(list(score = score, effect = "strong"))
  } else if (abs(score) > 0.4) {
    return(list(score = score, effect = "weak"))
  } else {
    return(list(score = score, effect = "neut"))
  }
}
# ----------------------------------------------------------------------- #


# modified from original
# ----------------------------------------------------------------------- #
snps.from.rsid2 <- function (rsid = NULL, dbSNP = NULL, search.genome = NULL) 
{
  if (is.null(rsid)) {
    stop("no RefSNP IDs have been found, please include RefSNP ID numbers")
  }
  if (!inherits(dbSNP, "SNPlocs")) {
    stop(paste0("dbSNP argument was not provided with a valid SNPlocs object.\n", 
                "Please run availible.SNPs() to check for availible SNPlocs"))
  }
  if (!inherits(search.genome, "BSgenome")) {
    stop(paste0("search.genome argument was not provided with a valid BSgenome object.\n", 
                "Run availible.genomes() and choose the appropriate BSgenome object"))
  }
  if (all(!grepl("rs", rsid))) {
    bad.names <- rsid[!grepl("rs", rsid)]
    stop(paste(paste(bad.names, collapse = " "), "are not rsids, perhaps you want to import your snps from a bed or vcf file with snps.from.file()?"))
  }
  rsid <- unique(rsid)
  rsid.grange <- as(snpsById(dbSNP, rsid, ifnotfound = "warning"), 
                    "GRanges")
  rsid.grange <- change.to.search.genome(rsid.grange, search.genome)
  rsid.grange <- GRanges(rsid.grange)
  rsid.refseq <- getSeq(search.genome, rsid.grange)
  rsid.grange$UCSC.reference <- as.character(rsid.refseq)
  rsid.grange <- sapply(split(rsid.grange, rsid.grange$RefSNP_id), 
                        function(snp) {
                          alt.allele <- determine.allele.from.ambiguous(snp$alleles_as_ambig, 
                                                                        snp$UCSC.reference)
                          if (length(alt.allele) > 1L) {
                            snp <- do.call("c", replicate(length(alt.allele), 
                                                          snp))
                            snp$UCSC.alternate <- alt.allele
                            names(snp) <- paste(snp$RefSNP_id, alt.allele, 
                                                sep = ":")
                          }
                          else {
                            snp$UCSC.alternate <- alt.allele
                            names(snp) <- snp$RefSNP_id
                          }
                          return(snp)
                        })
  rsid.grange <- unlist(do.call("GRangesList", rsid.grange), 
                        use.names = FALSE)
  colnames(mcols(rsid.grange)) <- c("RefSNP_id", "alleles_as_ambig", 
                                    "REF", "ALT")
  rsid.grange$REF <- DNAStringSet(rsid.grange$REF)
  rsid.grange$ALT <- DNAStringSet(rsid.grange$ALT)
  rsid.grange$alleles_as_ambig <- NULL
  colnames(mcols(rsid.grange))[1] <- "SNP_id"
  attributes(rsid.grange)$genome.package <- attributes(search.genome)$pkgname
  return(rsid.grange)
}
# ----------------------------------------------------------------------- #

# modified from original
# ----------------------------------------------------------------------- #
motifbreakR2 <- function (snpList, pwmList, threshold = 0.85, filterp = FALSE, 
                          method = "default", show.neutral = FALSE, verbose = FALSE, 
                          bkg = c(A = 0.25, C = 0.25, G = 0.25, T = 0.25), legacy.score = TRUE, 
                          BPPARAM = bpparam()) 
{
  if (.Platform$OS.type == "windows" && inherits(BPPARAM, "MulticoreParam")) {
    warning(paste0("Serial evaluation under effect, to achive parallel evaluation under\n", 
                   "Windows, please supply an alternative BPPARAM"))
  }
  cores <- bpnworkers(BPPARAM)
  num.snps <- length(snpList)
  if (num.snps < cores) {
    cores <- num.snps
  }
  if (is(BPPARAM, "SnowParam")) {
    # had to change this from bpstart alone, caused a crash
    register(bpstart(MulticoreParam(6)))
    cl <- bpbackend(BPPARAM)
    clusterEvalQ(cl, library("MotifDb"))
  }
  genome.package <- attributes(snpList)$genome.package
  if (requireNamespace(eval(genome.package), quietly = TRUE, 
                       character.only = TRUE)) {
    genome.bsgenome <- eval(parse(text = paste(genome.package, 
                                               genome.package, sep = "::")))
  }
  else {
    stop(paste0(eval(genome.package), " is the genome selected for this snp list and \n", 
                "is not present on your environment. Please load it and try again."))
  }
  pwms <- preparePWM(pwmList = pwmList, filterp = filterp, 
                     scoreThresh = threshold, bkg = bkg, method = method)
  k <- max(sapply(pwms$pwmList, ncol))
  snpList <- prepareVariants(fsnplist = snpList, genome.bsgenome = genome.bsgenome, 
                             max.pwm.width = k, legacy = legacy.score)
  snpList_cores <- split(as.list(rep(names(snpList), times = cores)), 
                         1:cores)
  for (splitr in seq_along(snpList)) {
    splitcores <- sapply(suppressWarnings(split(snpList[[splitr]], 
                                                1:cores)), list)
    for (splitcore in seq_along(snpList_cores)) {
      snpList_cores[[splitcore]][[splitr]] <- splitcores[[splitcore]]
      names(snpList_cores[[splitcore]])[splitr] <- names(snpList)[splitr]
    }
  }
  snpList <- snpList_cores
  rm(snpList_cores)
  x <- lapply(snpList, scoreSnpList, pwmList = pwms$pwmList, 
              threshold = pwms$pwmThreshold, pwmList.pc = pwms$pwmListPseudoCount, 
              pwmRanges = pwms$pwmRange, method = method, bkg = bkg, 
              show.neutral = show.neutral, verbose = ifelse(cores == 
                                                              1, verbose, FALSE), genome.bsgenome = genome.bsgenome, 
              filterp = filterp)
  if (inherits(x, "try-error")) {
    if (is(BPPARAM, "SnowParam")) {
      bpstop(BPPARAM)
    }
    stop(attributes(x)$condition)
  }
  if (is(BPPARAM, "SnowParam")) {
    bpstop(BPPARAM)
  }
  drops <- sapply(x, is.null)
  x <- x[!drops]
  pwmList <- pwms$pwmList
  pwmList@listData <- lapply(pwms$pwmList, function(pwm) {
    pwm <- pwm[c("A", "C", "G", "T"), ]
    return(pwm)
  })
  pwmList.pc <- lapply(pwms$pwmListPseudoCount, function(pwm) {
    pwm <- pwm[c("A", "C", "G", "T"), ]
    return(pwm)
  })
  if (length(x) > 1) {
    x <- unlist(GRangesList(unname(x)))
    snpList <- unlist(GRangesList(lapply(snpList, `[[`, "fsnplist")), 
                      use.names = F)
    x <- x[order(match(names(x), names(snpList)), x$geneSymbol), 
    ]
    attributes(x)$genome.package <- genome.package
    attributes(x)$motifs <- pwmList[mcols(pwmList)$providerId %in% 
                                      unique(x$providerId) & mcols(pwmList)$providerName %in% 
                                      unique(x$providerName), ]
    attributes(x)$scoremotifs <- pwmList.pc[names(attributes(x)$motifs)]
  }
  else {
    if (length(x) == 1L) {
      x <- x[[1]]
      attributes(x)$genome.package <- genome.package
      attributes(x)$motifs <- pwmList[mcols(pwmList)$providerId %in% 
                                        unique(x$providerId) & mcols(pwmList)$providerName %in% 
                                        unique(x$providerName), ]
      attributes(x)$scoremotifs <- pwmList.pc[names(attributes(x)$motifs)]
    }
    else {
      warning("No SNP/Motif Interactions reached threshold")
      x <- NULL
    }
  }
  if (verbose && cores > 1) {
    if (is.null(x)) {
      message(paste("reached end of SNPs list length =", 
                    num.snps, "with 0 potentially disruptive matches to", 
                    length(unique(x$geneSymbol)), "of", length(pwmList), 
                    "motifs."))
    }
    else {
      message(paste("reached end of SNPs list length =", 
                    num.snps, "with", length(x), "potentially disruptive matches to", 
                    length(unique(x$geneSymbol)), "of", length(pwmList), 
                    "motifs."))
    }
  }
  return(x)
}
# ----------------------------------------------------------------------- #

# ----------------------------------------------------------------------- #
change.to.search.genome <- function(granges.object, search.genome) {
  sequence <- seqlevels(granges.object)
  ## sequence is in UCSC format and we want NCBI style
  newStyle <- mapSeqlevels(sequence,seqlevelsStyle(search.genome))
  newStyle <- newStyle[complete.cases(newStyle)] # removing NA cases.
  ## rename the seqlevels
  granges.object <- renameSeqlevels(granges.object,newStyle)
  seqlevels(granges.object) <- seqlevelsInUse(granges.object)
  seqinfo(granges.object) <- keepSeqlevels(seqinfo(search.genome),
                                           value = seqlevelsInUse(granges.object))
  return(granges.object)
}
# ----------------------------------------------------------------------- #


# ----------------------------------------------------------------------- #
determine.allele.from.ambiguous <- function(ambiguous.allele, known.allele) {
  neucleotide.ambiguity.code <- list(Y = c("C", "T"), R = c("A", "G"), W = c("A", "T"),
                                     S = c("G", "C"), K = c("T", "G"), M = c("C", "A"),
                                     D = c("A", "G", "T"), V = c("A", "C", "G"),
                                     H = c("A", "C", "T"), B = c("C", "G", "T"),
                                     N = c("A", "C", "G", "T"))
  specnac <- neucleotide.ambiguity.code[[ambiguous.allele]]
  unknown.allele <- specnac[-grep(known.allele, specnac)]
  return(unknown.allele)
}
# ----------------------------------------------------------------------- #

# ----------------------------------------------------------------------- #
strSort <- function(x) {
  sapply(lapply(strsplit(x, NULL), sort), paste, collapse = "")
}
# ----------------------------------------------------------------------- #

# ----------------------------------------------------------------------- #
DNAmotifAlignment.2snp <- function(pwms, result) {
  from <- min(sapply(result$motifPos, `[`, 1))
  to <- max(sapply(result$motifPos, `[`, 2))
  #  pos <- mcols(result)$motifPos
  #  pos <- pos[as.logical(strand(result) == "+")][1]
  for (pwm.i in seq_along(pwms)) {
    #pwm <- pwms[[pwm.i]]@mat
    ## get pwm info from result data
    pwm.name <- pwms[[pwm.i]]@name
    pwm.name <- str_replace(pwm.name, pattern = "-:rc$", replacement = "")
    pwm.name <- str_replace(pwm.name, pattern = "-:r$", replacement = "")
    pwm.info <- attributes(result)$motifs
    pwm.id <- mcols(pwm.info[pwm.name, ])$providerId
    pwm.name <- mcols(pwm.info[pwm.name, ])$providerName
    mresult <- result[result$providerId == pwm.id & result$providerName == pwm.name, ]
    mstart <- mresult$motifPos[[1]][1]
    mend <- mresult$motifPos[[1]][2]
    # browser()
    if ((mcols(mresult)$varType == "Insertion" & mcols(mresult)$alleleDiff < 0) |
        (mcols(mresult)$varType == "Deletion" & mcols(mresult)$alleleDiff > 0)) {
      new.mat <- cbind(pwms[[pwm.i]]@mat[, 1:abs(mresult$motifPos[[1]][1])],
                       matrix(c(0.25, 0.25, 0.25, 0.25), ncol = length(mcols(mresult)$altPos[[1]]), nrow = 4),
                       pwms[[pwm.i]]@mat[, (abs(mresult$motifPos[[1]][1]) + 1):ncol(pwms[[pwm.i]]@mat)])
      pwms[[pwm.i]]@mat <- new.mat
      start.offset <- mstart - from
      end.offset <- to - mend
    } else {
      if (mstart < 0 | from > 0) {
        start.offset <- mstart - from
      } else {
        start.offset <- (mstart - 1) - from
      }
      if (mend > 0 | to < 0) {
        end.offset <- to - mend
      } else {
        end.offset <- to - (mend + 1)
      }
    }
    if (start.offset > 0) {
      pwms[[pwm.i]] <- addBlank(x = pwms[[pwm.i]], n = start.offset, b = FALSE)
    }
    if (end.offset > 0) {
      pwms[[pwm.i]] <- addBlank(x = pwms[[pwm.i]], n = end.offset, b = TRUE)
    }
  }
  return(pwms)
}
# ----------------------------------------------------------------------- #


# ----------------------------------------------------------------------- #
plotMB2 <- function (results, rsid, reverseMotif = TRUE, effect = c("strong", 
                                                         "weak")) 
{
  motif.starts <- sapply(results$motifPos, `[`, 1)
  motif.starts <- start(results) + motif.starts
  motif.starts <- order(motif.starts)
  results <- results[motif.starts]
  g <- genome(results)[[1]]
  result <- results[names(results) %in% rsid]
  result <- result[order(sapply(result$motifPos, min), sapply(result$motifPos, 
                                                              max)), ]
  result <- result[result$effect %in% effect]
  chromosome <- as.character(seqnames(result))[[1]]
  genome.package <- attributes(result)$genome.package
  genome.bsgenome <- eval(parse(text = genome.package))
  seq.len <- max(length(result$REF[[1]]), length(result$ALT[[1]]))
  distance.to.edge <- max(abs(c(sapply(result$motifPos, min), 
                                sapply(result$motifPos, max)))) + 4
  from <- start(result)[[1]] - distance.to.edge + 1
  to <- end(result)[[1]] + distance.to.edge
  pwmList <- attributes(result)$motifs
  pwm.names <- result$providerId
  results_motifs <- paste0(result$providerId, result$providerName)
  list_motifs <- paste0(mcols(pwmList)$providerId, mcols(pwmList)$providerName)
  pwms <- pwmList <- pwmList[match(results_motifs, list_motifs)]
  if (reverseMotif) {
    for (pwm.i in seq_along(pwms)) {
      pwm.name <- names(pwms[pwm.i])
      pwm.id <- mcols(pwms[pwm.name, ])$providerId
      pwm.name.f <- mcols(pwms[pwm.name, ])$providerName
      doRev <- as.logical(strand(result[result$providerId == 
                                          pwm.id & result$providerName == pwm.name.f, ]) == 
                            "-")
      if (doRev) {
        pwm <- pwms[[pwm.i]]
        pwm <- pwm[, rev(1:ncol(pwm))]
        rownames(pwm) <- c("T", "G", "C", "A")
        pwm <- pwm[c("A", "C", "G", "T"), ]
        pwms[[pwm.i]] <- pwm
        names(pwms)[pwm.i] <- paste0(names(pwms)[pwm.i], 
                                     "-:rc")
      }
    }
  }
  else {
    for (pwm.i in seq_along(pwms)) {
      pwm.name <- names(pwms[pwm.i])
      pwm.id <- mcols(pwms[pwm.name, ])$providerId
      pwm.name.f <- mcols(pwms[pwm.name, ])$providerName
      doRev <- as.logical(strand(result[result$providerId == 
                                          pwm.id & result$providerName == pwm.name.f, ]) == 
                            "-")
      if (doRev) {
        pwm <- pwms[[pwm.i]]
        pwm <- pwm[, rev(1:ncol(pwm))]
        pwms[[pwm.i]] <- pwm
        names(pwms)[pwm.i] <- paste0(names(pwms)[pwm.i], 
                                     "-:r")
      }
    }
  }
  pwms <- lapply(names(pwms), function(x, pwms = pwms) {
    new("pfm", mat = pwms[[x]], name = x)
  }, pwms)
  pwms <- DNAmotifAlignment.2snp(pwms, result)
  pwmwide <- max(sapply(pwms, function(x) {
    ncol(x@mat)
  }))
  markerStart <- result$motifPos[[1]][1]
  if (markerStart > 0) {
    markerEnd <- length(result$altPos[[1]]) + 1
    markerEnd <- markerEnd - markerStart
    markerStart <- 1
  }
  else {
    markerStart <- -1 * markerStart
    markerEnd <- markerStart + length(result$altPos[[1]])
    if (result$varType[[1]] %in% c("Other", "SNV")) {
      markerStart <- markerStart + 1
    }
  }
  varType <- result$varType[[1]]
  varType <- switch(varType, Deletion = "firebrick", Insertion = "springgreen4", 
                    Other = "gray13")
  markerRect <- new("marker", type = "rect", start = markerStart, 
                    stop = markerEnd, gp = gpar(lty = 2, fill = NA, lwd = 3, 
                                                col = varType))
  for (pwm.i in seq_along(pwms)) {
    pwms[[pwm.i]]@markers <- list(markerRect)
  }
  ideoT <- try(IdeogramTrack(genome = g, chromosome = chromosome), 
               silent = TRUE)
  if (inherits(ideoT, "try-error")) {
    backup.band <- data.frame(chrom = chromosome, chromStart = 0, 
                              chromEnd = length(genome.bsgenome[[chromosome]]), 
                              name = chromosome, gieStain = "gneg")
    ideoT <- IdeogramTrack(genome = g, chromosome = chromosome, 
                           bands = backup.band)
  }
  altseq <- genome.bsgenome[[chromosome]]
  at <- IRanges(start = start(result[1]), width = width(result[1]))
  if (result$varType[[1]] == "Deletion") {
    reflen <- length(result$REF[[1]])
    addedN <- DNAString(paste0(rep.int(".", reflen), collapse = ""))
    addedN <- replaceLetterAt(addedN, at = (1:reflen)[-result$altPos[[1]]], 
                              result$ALT[[1]])
    axisT <- GenomeAxisTrack(exponent = 0)
    seqT <- SequenceTrack(genome.bsgenome, fontcolor = colorset("DNA", 
                                                                "auto"))
    altseq <- replaceAt(x = altseq, at = at, addedN)
  }
  else if (result$varType[[1]] == "Insertion") {
    altlen <- length(result$ALT[[1]])
    addedN <- DNAString(paste0(rep.int(".", altlen), collapse = ""))
    addedN <- replaceLetterAt(addedN, at = (1:altlen)[-result$altPos[[1]]], 
                              result$REF[[1]])
    refseq <- genome.bsgenome[[chromosome]]
    refseq <- DNAStringSet(replaceAt(x = refseq, at = at, 
                                     addedN))
    altseq <- replaceAt(x = altseq, at = at, result$ALT[[1]])
    names(refseq) <- chromosome
    seqT <- SequenceTrack(refseq, fontcolor = c(colorset("DNA", 
                                                         "auto"), N = "#FFFFFF", . = "#FFE3E6"), chromosome = chromosome)
  }
  else {
    axisT <- GenomeAxisTrack(exponent = 0)
    altseq <- replaceAt(x = altseq, at = at, result$ALT[[1]])
    seqT <- SequenceTrack(genome.bsgenome, fontcolor = colorset("DNA", 
                                                                "auto"))
  }
  altseq <- DNAStringSet(altseq)
  names(altseq) <- chromosome
  seqAltT <- SequenceTrack(altseq, fontcolor = c(colorset("DNA", 
                                                          "auto"), N = "#FFFFFF", . = "#FFE3E6"), chromosome = chromosome)
  histart <- start(result[1]) + min(result[1]$altPos[[1]]) - 
    2
  histart <- ifelse(result[1]$varType %in% c("Other", "SNV"), 
                    histart + 1, histart)
  hiend <- start(result[1]) + min(result[1]$altPos[[1]]) - 
    2 + length(result[1]$altPos[[1]])
  hiT <- HighlightTrack(trackList = list(seqT, seqAltT), start = histart, 
                        end = hiend, chromosome = chromosome)
  selectingfun <- selcor
  detailfun <- addPWM.stack
  motif_ids <- names(pwmList)
  names(motif_ids) <- mcols(pwmList)$providerName
  for (mymotif_i in seq_along(result)) {
    mymotif <- result[mymotif_i]
    start(mymotif) <- start(mymotif) + min(mymotif$altPos[[1]]) - 
      1
    width(mymotif) <- length(mymotif$altPos[[1]])
    variant.start <- start(mymotif)
    variant.end <- end(mymotif)
    if (mymotif$motifPos[[1]][1] < 0) {
      start(mymotif) <- start(mymotif) + (mymotif$motifPos[[1]][1])
    }
    else {
      start(mymotif) <- start(mymotif) + (mymotif$motifPos[[1]][1] - 
                                            1)
    }
    if (mymotif$motifPos[[1]][2] < 0) {
      end(mymotif) <- end(mymotif) + (mymotif$motifPos[[1]][2] + 
                                        1)
    }
    else {
      end(mymotif) <- end(mymotif) + (mymotif$motifPos[[1]][2])
    }
    if ((result[mymotif_i]$varType == "Deletion" & result[mymotif_i]$alleleDiff > 
         0) | (result[mymotif_i]$varType == "Insertion" & 
               result[mymotif_i]$alleleDiff < 0)) {
      mymotif <- c(mymotif, mymotif)
      end(mymotif)[1] <- variant.start - 1
      start(mymotif)[2] <- variant.end + 1
      mymotif[which.min(width(mymotif))]$motifPos <- NA
    }
    if (exists("mres")) {
      mres <- c(mres, mymotif)
    }
    else {
      mres <- mymotif
    }
  }
  result <- mres
  rm(mres)
  motif_ids <- motif_ids[result$providerName]
  presult <- result
  strand(presult) <- "*"
  pres_cols <- DataFrame(feature = ifelse(!is.na(result$motifPos), 
                                          paste(result$geneSymbol, "motif", sep = "_"), ""), group = result$providerName, 
                         id = motif_ids)
  presult <- GRanges(seqnames = seqnames(result[1]), ranges = ranges(result))
  mcols(presult) <- pres_cols
  motifT <- AnnotationTrack(presult, fun = detailfun, detailsFunArgs = list(pwm_stack = pwms), 
                            name = names(result)[[1]], selectFun = selectingfun, 
                            reverseStacking = FALSE, stacking = "squish")
  if (exists("axisT")) {
    track_list <- list(ideoT, motifT, hiT, axisT)
  }
  else {
    track_list <- list(ideoT, motifT, hiT)
  }
  plotTracks(track_list, from = from, to = to, showBandId = TRUE, 
             cex.main = 0.8, col.main = "darkgrey", add53 = TRUE, 
             labelpos = "below", chromosome = chromosome, fontcolor.item = "black", 
             collapse = FALSE, min.width = 1, featureAnnotation = "feature", 
             cex.feature = 0.8, details.size = 0.85, detailsConnector.pch = NA, 
             detailsConnector.lty = 0, shape = "box", cex.group = 0.8, 
             fonts = c("sans", "Helvetica"))
  return(invisible(NULL))
}
# ----------------------------------------------------------------------- #


# ----------------------------------------------------------------------- #
preparePWM <- function(pwmList = pwmList,
                       filterp = filterp,
                       bkg = bkg,
                       scoreThresh = threshold,
                       method = "default") {
  
  bkg <- bkg[c('A', 'C', 'G', 'T')]
  
  scounts <- as.integer(mcols(pwmList)$sequenceCount)
  scounts[is.na(scounts)] <- 20L
  pwmList.pc <- Map(function(pwm, scount) {
    pwm <- (pwm * scount + 0.25)/(scount + 1)
  }, pwmList, scounts)
  if (method == "ic") {
    pwmOmegas <- lapply(pwmList.pc, function(pwm, b=bkg) {
      omegaic <- colSums(pwm * log2(pwm/b))
    })
  }
  if (method == "default") {
    pwmOmegas <- lapply(pwmList.pc, function(pwm) {
      omegadefault <- colMaxs(pwm) - colMins(pwm)
    })
  }
  if (method == "log") {
    pwmList.pc <- lapply(pwmList.pc, function(pwm, b) {
      pwm <- log(pwm) - log(b)
    }, b = bkg)
    pwmOmegas <- 1
  }
  if (method == "notrans") {
    pwmOmegas <- 1
  }
  pwmList.pc <- Map(function(pwm, omega) {
    if (length(omega) == 1 && omega == 1) {
      return(pwm)
    } else {
      omegamatrix <- matrix(rep(omega, 4), nrow = 4, byrow = TRUE)
      pwm <- pwm * omegamatrix
    }
  }, pwmList.pc, pwmOmegas)
  if (filterp) {
    pwmRanges <- Map(function(pwm, omega) {
      x <- colSums(colRanges(pwm))
      return(x)
    }, pwmList.pc, pwmOmegas)
    pwmList.pc2 <- lapply(pwmList.pc, round, digits = 2)
    pwmThresh <- lapply(pwmList.pc2, TFMpv2sc, pvalue = scoreThresh, bg = bkg, type = "PWM")
    pwmThresh <- Map("+", pwmThresh, -0.02)
  } else {
    pwmRanges <- Map(function(pwm, omega) {
      x <- colSums(colRanges(pwm))
      return(x)
    }, pwmList.pc, pwmOmegas)
    pwmThresh <- rep.int(scoreThresh, times = length(pwmRanges))
  }
  pwmList@listData <- lapply(pwmList, function(pwm) {
    pwm <- rbind(pwm, N = 0)
    colnames(pwm) <- as.character(1:ncol(pwm))
    return(pwm) })
  pwmList.pc <- lapply(pwmList.pc, function(pwm) {
    pwm <- rbind(pwm, N = 0)
    colnames(pwm) <- as.character(1:ncol(pwm))
    return(pwm) })
  return(list(pwmList = pwmList,
              pwmListPseudoCount = pwmList.pc,
              pwmRange = pwmRanges,
              pwmThreshold = pwmThresh))
}
# ----------------------------------------------------------------------- #


# ----------------------------------------------------------------------- #
scoreSnpList <- function(fsnplist, pwmList, method = "default", bkg = NULL,
                         threshold = 1e-3, show.neutral = FALSE, verbose = FALSE, legacy = FALSE,
                         genome.bsgenome=NULL, pwmList.pc = NULL, pwmRanges = NULL, filterp=TRUE) {
  k <- max(sapply(pwmList, ncol))
 
  snp.sequence.alt <- fsnplist$alt.seq
  snp.sequence.ref <- fsnplist$ref.seq
  fsnplist <- fsnplist$fsnplist
  
  res.el.e <- new.env()
  for (snp.map.i in seq_along(snp.sequence.alt)) {
    snp.ref <- snp.sequence.ref[[snp.map.i]]
    snp.alt <- snp.sequence.alt[[snp.map.i]]
    ref.len <- nchar(fsnplist[snp.map.i]$REF)
    alt.len <- nchar(fsnplist[snp.map.i]$ALT)
    alt.loc <- fsnplist[snp.map.i]$ALT_loc[[1]]
    res.el <- rep(fsnplist[snp.map.i], length(pwmList))
    res.el$motifPos <- as.integer(NA)
    res.el$motifID <- mcols(pwmList)$providerID
    res.el$geneSymbol <- mcols(pwmList)$geneSymbol
    res.el$dataSource <- mcols(pwmList)$dataSource
    res.el$providerName <- mcols(pwmList)$providerName
    res.el$providerId <- mcols(pwmList)$providerId
    res.el$seqMatch <- as.character(NA)
    res.el$pctRef <- as.numeric(NA)
    res.el$pctAlt <- as.numeric(NA)
    if (filterp) {
      res.el$scoreRef <- as.numeric(NA)
      res.el$scoreAlt <- as.numeric(NA)
      res.el$Refpvalue <- as.numeric(NA)
      res.el$Altpvalue <- as.numeric(NA)
    }
    if (ref.len > 1 | alt.len > 1 | !legacy) {
      res.el$altPos <- as.numeric(NA)
      res.el$alleleDiff <- as.numeric(NA)
    } else {
      res.el$snpPos <- as.integer(NA)
      res.el$alleleRef <- as.numeric(NA)
      res.el$alleleAlt <- as.numeric(NA)
    }
    res.el$effect <- as.character(NA)
    for (pwm.i in seq_along(pwmList)) {
      pwm.basic <- pwmList[[pwm.i]]
      pwm <- pwmList.pc[[pwm.i]]
      len <- ncol(pwm)
      thresh <- threshold[[pwm.i]]
      seq.start <- min(alt.loc)
      seq.len <- length(alt.loc)
      alt.range <- ref.range <- (k - (ncol(pwm) - seq.start)):(k + ncol(pwm) + seq.start + seq.len - 2)
      if (!show.neutral & identical(snp.ref[ref.range], snp.alt[alt.range])) next()
      seq.remove <- ref.len - alt.len
      if (seq.remove < 0) {
        ref.range <- ref.range[1:(length(ref.range) + seq.remove)]
      } else {
        alt.range <- alt.range[1:(length(alt.range) - seq.remove)]
      }
      ref.windows <- scoreSeqWindows(ppm = pwm, seq = snp.ref[ref.range])
      alt.windows <- scoreSeqWindows(ppm = pwm, seq = snp.alt[alt.range])
      if (!filterp) {
        ref.windows <- (ref.windows - pwmRanges[[pwm.i]][1]) / (pwmRanges[[pwm.i]][2] - pwmRanges[[pwm.i]][1])
        alt.windows <- (alt.windows - pwmRanges[[pwm.i]][1]) / (pwmRanges[[pwm.i]][2] - pwmRanges[[pwm.i]][1])
      }
      if (any(alt.windows > thresh) | any(ref.windows > thresh)) {
        hit.alt <- maxThresholdWindows(alt.windows)
        hit.ref <- maxThresholdWindows(ref.windows)
        bigger <- ref.windows[hit.ref$strand, hit.ref$window] >= alt.windows[hit.alt$strand, hit.alt$window]
        if (bigger) {
          hit <- hit.ref
        } else {
          hit <- hit.alt
        }
      } else {
        hit.alt <- list(window = 0L, strand = 0L)
        hit.ref <- list(window = 0L, strand = 0L)
        hit <- NULL
      }
      if (!show.neutral) {
        if (identical(alt.windows[hit.alt$strand, hit.alt$window],
                      ref.windows[hit.ref$strand, hit.ref$window])) next()
      }
      if (!is.null(hit)) {
        result <- res.el[pwm.i]
        uniquename <- paste(names(result), result$dataSource, result$providerName, result$providerId, sep = "%%")
        if (nchar(result$REF) > 1 | nchar(result$ALT) > 1) {
          allelR <- ref.windows[hit.ref$strand, hit.ref$window]
          allelA <- alt.windows[hit.alt$strand, hit.alt$window]
          scorediff <- varEff(allelR, allelA)
          effect <- scorediff$effect
          score <- scorediff$score
          #ref.pos <- ref.range[(ncol(pwm) - (seq.start - 1)):((ncol(pwm) + ref.len - seq.start))]
          ref.pos <- k:(k + nchar(result$REF) - 1L)
          #alt.pos <- alt.range[(ncol(pwm) - (seq.start - 1)):((ncol(pwm) + alt.len - seq.start))]
          alt.pos <- k:(k + nchar(result$ALT) - 1L)
          if (effect == "neut") {
            if (show.neutral) {
              res.el.e[[uniquename]] <- updateResultsIndel(result,
                                                           snp.ref, snp.alt,
                                                           ref.pos, alt.pos,
                                                           hit.ref, hit.alt,
                                                           ref.windows, alt.windows,
                                                           score, effect, len,
                                                           k, pwm, calcp = filterp)
            }
          } else {
            res.el.e[[uniquename]] <- updateResultsIndel(result,
                                                         snp.ref, snp.alt,
                                                         ref.pos, alt.pos,
                                                         hit.ref, hit.alt,
                                                         ref.windows, alt.windows,
                                                         score, effect, len,
                                                         k, pwm, calcp = filterp)
          }
        } else {
          snp.pos <- len - hit$window + 1L
          if (legacy) {
            if (hit$strand == 1) {
              allelR <- pwm.basic[as.character(result$REF), snp.pos]
              allelA <- pwm.basic[as.character(result$ALT), snp.pos]
            } else {
              allelR <- pwm.basic[as.character(complement(result$REF)), snp.pos]
              allelA <- pwm.basic[as.character(complement(result$ALT)), snp.pos]
            }
          } else {
            allelR <- ref.windows[hit.ref$strand, hit.ref$window]
            allelA <- alt.windows[hit.alt$strand, hit.alt$window]
          }
          scorediff <- varEff(allelR, allelA)
          effect <- scorediff$effect
          score <- scorediff$score
          if (effect == "neut") {
            if (show.neutral) {
              if (legacy) {
                res.el.e[[uniquename]] <- updateResultsSnv(result, snp.ref[ref.range], snp.pos,
                                                           hit, ref.windows, alt.windows,
                                                           allelR, allelA, effect, len,
                                                           k, pwm, calcp = filterp)
              } else {
                res.el.e[[uniquename]] <- updateResultsIndel(result,
                                                             snp.ref, snp.alt,
                                                             21, 21,
                                                             hit.ref, hit.alt,
                                                             ref.windows, alt.windows,
                                                             score, effect, len,
                                                             k, pwm, calcp = filterp)
              }
            }
          } else {
            if (legacy) {
              res.el.e[[uniquename]] <- updateResultsSnv(result, snp.ref[ref.range], snp.pos,
                                                         hit, ref.windows, alt.windows,
                                                         allelR, allelA, effect, len,
                                                         k, pwm, calcp = filterp)
            } else {
              res.el.e[[uniquename]] <- updateResultsIndel(result,
                                                           snp.ref, snp.alt,
                                                           21, 21,
                                                           hit.ref, hit.alt,
                                                           ref.windows, alt.windows,
                                                           score, effect, len,
                                                           k, pwm, calcp = filterp)
            }
          }
        }
      }
    }
  }
  resultSet <- unlist(GRangesList(as.list.environment(res.el.e)), use.names = FALSE)
  if (length(resultSet) < 1) {
    if (verbose) {
      message(paste("reached end of SNPs list length =", length(fsnplist),
                    "with 0 potentially disruptive matches to", length(unique(resultSet$geneSymbol)),
                    "of", length(pwmList), "motifs."))
    }
    return(NULL)
  } else {
    if ("ALT_loc" %in% names(mcols(resultSet))) mcols(resultSet)$ALT_loc <- NULL
    max.match <- max(vapply(str_locate_all(resultSet$seqMatch, "\\w"), max, integer(1)))
    min.match <- min(vapply(str_locate_all(resultSet$seqMatch, "\\w"), min, integer(1)))
    resultSet$seqMatch <- str_sub(resultSet$seqMatch,
                                  start = min.match + 1,
                                  end = max.match + 1)
    if (verbose) {
      message(paste("reached end of SNPs list length =", length(fsnplist),
                    "with", length(resultSet), "potentially disruptive matches to", length(unique(resultSet$geneSymbol)),
                    "of", length(pwmList), "motifs."))
    }
    return(resultSet)
  }
}
# ----------------------------------------------------------------------- #


# ----------------------------------------------------------------------- #
updateResultsSnv <- function(result, snp.seq, snp.pos, hit, ref.windows, alt.windows,
                             allelR, allelA, effect, len, k, pwm, calcp) {
  strand.opt <- c("+", "-")
  strand(result) <- strand.opt[[hit$strand]]
  hit$window <- as.integer(hit$window)
  mresult <- mcols(result)
  mresult[["snpPos"]] <- start(result)
  mresult[["motifPos"]] <- as.integer(snp.pos)
  matchs <- snp.seq
  seq.pos <- snp.pos + hit$window - 1
  matchs[-(seq.pos)] <- tolower(matchs[-(seq.pos)])
  matchs <- paste(matchs, collapse = "")
  mresult[["seqMatch"]] <- str_pad(matchs, width = k * 2, side = "both")
  start(result) <- start(result) - snp.pos + 1
  end(result) <- end(result) - snp.pos + len
  if (calcp) {
    mresult[["scoreRef"]] <- ref.windows[hit$strand, hit$window]
    mresult[["scoreAlt"]] <- alt.windows[hit$strand, hit$window]
    mresult[["Refpvalue"]] <- NA
    mresult[["Altpvalue"]] <- NA
    pwmrange <- colSums(colRanges(pwm))
    mresult[["pctRef"]] <- (mresult[["scoreRef"]] - pwmrange[[1]]) / (pwmrange[[2]] - pwmrange[[1]])
    mresult[["pctAlt"]] <- (mresult[["scoreAlt"]] - pwmrange[[1]]) / (pwmrange[[2]] - pwmrange[[1]])
  } else {
    mresult[["pctRef"]] <- ref.windows[hit$strand, hit$window]
    mresult[["pctAlt"]] <- alt.windows[hit$strand, hit$window]
  }
  mresult[["alleleRef"]] <- allelR
  mresult[["alleleAlt"]] <- allelA
  mresult[["effect"]] <- effect
  mcols(result) <- mresult
  return(result)
}
# ----------------------------------------------------------------------- #


# ----------------------------------------------------------------------- #
updateResultsIndel <- function(result,
                               ref.seq, alt.seq,
                               ref.pos, alt.pos,
                               hit.ref, hit.alt,
                               ref.windows, alt.windows,
                               score, effect, len, k, pwm, calcp) {
  strand.opt <- c("+", "-")
  if (score > 0L) {
    best.hit <- hit.alt
    matchs <- alt.seq
    snp.pos <- alt.pos
  } else {
    best.hit <- hit.ref
    matchs <- ref.seq
    snp.pos <- ref.pos
  }
  strand(result) <- strand.opt[[best.hit$strand]]
  best.hit$window <- as.integer(best.hit$window)
  mresult <- mcols(result)
  alt_loc <- range(mresult$ALT_loc)
  ref_start <- (1 - alt_loc[[1]])
  ref_start <- ifelse(ref_start <= 0, ref_start - 1, ref_start)
  motif.start <- (alt_loc[[1]]) + (-len) + (best.hit$window) + ref_start
  motif.start <- ifelse(motif.start >= 0, motif.start + 1, motif.start)
  if ((mresult$varType == "Insertion" & score < 0) |
      (mresult$varType == "Deletion" & score > 0)) {
    motif.end <- motif.start + len
  } else {
    if (motif.start > 0) {
      motif.end <- len - length(motif.start:length(alt_loc[1]:alt_loc[2]))
    } else {
      motif.end <- motif.start + len - length(alt_loc[1]:alt_loc[2])
    }
  }
  motif.end <- ifelse(motif.end <= 0, motif.end - 1, motif.end)
  mresult$motifPos <- list(c(motif.start, motif.end))
  mresult$altPos <- mresult$ALT_loc
  seq.range <- (k - (len - alt_loc[[1]])):(k + len + alt_loc[[2]] - 2)
  matchs[-(snp.pos)] <- tolower(matchs[-(snp.pos)])
  matchs <- paste(matchs[seq.range], collapse = "")
  mresult[["seqMatch"]] <- str_pad(matchs, width = (k * 2) + alt_loc[[2]], side = "both")
  if (calcp) {
    mresult[["scoreRef"]] <- ref.windows[hit.ref$strand, hit.ref$window]
    mresult[["scoreAlt"]] <- alt.windows[hit.alt$strand, hit.alt$window]
    mresult[["Refpvalue"]] <- NA
    mresult[["Altpvalue"]] <- NA
    pwmrange <- colSums(colRanges(pwm[-5,]))
    mresult[["pctRef"]] <- (mresult[["scoreRef"]] - pwmrange[[1]]) / (pwmrange[[2]] - pwmrange[[1]])
    mresult[["pctAlt"]] <- (mresult[["scoreAlt"]] - pwmrange[[1]]) / (pwmrange[[2]] - pwmrange[[1]])
  } else {
    mresult[["pctRef"]] <- ref.windows[hit.ref$strand, hit.ref$window]
    mresult[["pctAlt"]] <- alt.windows[hit.alt$strand, hit.alt$window]
  }
  mresult[["alleleDiff"]] <- score
  mresult[["effect"]] <- effect
  mcols(result) <- mresult
  return(result)
}
# ----------------------------------------------------------------------- #


# ----------------------------------------------------------------------- #
preparePWM <- function(pwmList = pwmList,
                       filterp = filterp,
                       bkg = bkg,
                       scoreThresh = threshold,
                       method = "default") {
  
  bkg <- bkg[c('A', 'C', 'G', 'T')]
  
  scounts <- as.integer(mcols(pwmList)$sequenceCount)
  scounts[is.na(scounts)] <- 20L
  pwmList.pc <- Map(function(pwm, scount) {
    pwm <- (pwm * scount + 0.25)/(scount + 1)
  }, pwmList, scounts)
  if (method == "ic") {
    pwmOmegas <- lapply(pwmList.pc, function(pwm, b=bkg) {
      omegaic <- colSums(pwm * log2(pwm/b))
    })
  }
  if (method == "default") {
    pwmOmegas <- lapply(pwmList.pc, function(pwm) {
      omegadefault <- colMaxs(pwm) - colMins(pwm)
    })
  }
  if (method == "log") {
    pwmList.pc <- lapply(pwmList.pc, function(pwm, b) {
      pwm <- log(pwm) - log(b)
    }, b = bkg)
    pwmOmegas <- 1
  }
  if (method == "notrans") {
    pwmOmegas <- 1
  }
  pwmList.pc <- Map(function(pwm, omega) {
    if (length(omega) == 1 && omega == 1) {
      return(pwm)
    } else {
      omegamatrix <- matrix(rep(omega, 4), nrow = 4, byrow = TRUE)
      pwm <- pwm * omegamatrix
    }
  }, pwmList.pc, pwmOmegas)
  if (filterp) {
    pwmRanges <- Map(function(pwm, omega) {
      x <- colSums(colRanges(pwm))
      return(x)
    }, pwmList.pc, pwmOmegas)
    pwmList.pc2 <- lapply(pwmList.pc, round, digits = 2)
    pwmThresh <- lapply(pwmList.pc2, TFMpv2sc, pvalue = scoreThresh, bg = bkg, type = "PWM")
    pwmThresh <- Map("+", pwmThresh, -0.02)
  } else {
    pwmRanges <- Map(function(pwm, omega) {
      x <- colSums(colRanges(pwm))
      return(x)
    }, pwmList.pc, pwmOmegas)
    pwmThresh <- rep.int(scoreThresh, times = length(pwmRanges))
  }
  pwmList@listData <- lapply(pwmList, function(pwm) {
    pwm <- rbind(pwm, N = 0)
    colnames(pwm) <- as.character(1:ncol(pwm))
    return(pwm) })
  pwmList.pc <- lapply(pwmList.pc, function(pwm) {
    pwm <- rbind(pwm, N = 0)
    colnames(pwm) <- as.character(1:ncol(pwm))
    return(pwm) })
  return(list(pwmList = pwmList,
              pwmListPseudoCount = pwmList.pc,
              pwmRange = pwmRanges,
              pwmThreshold = pwmThresh))
}
# ----------------------------------------------------------------------- #


# ----------------------------------------------------------------------- #
addPWM.stack <- function(identifier, index, GdObject, pwm_stack, ...) {
  plotMotifLogoStack.3(pwm_stack)
}
# ----------------------------------------------------------------------- #

# ----------------------------------------------------------------------- #
selcor <- function(identifier, index, GdObject, ... ) {
  if (identical(index, 1L)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
# ----------------------------------------------------------------------- #

# ----------------------------------------------------------------------- #
maxThresholdWindows <- function(window.frame) {
  start.ind <- as.integer(colnames(window.frame)[1]) - 1L
  max.win <- arrayInd(which.max(window.frame), dim(window.frame))
  return(list(window = as.integer(colnames(window.frame)[max.win[, 2] + start.ind]),
              strand = c(1, 2)[max.win[, 1]]))
}
# ----------------------------------------------------------------------- #

# ----------------------------------------------------------------------- #
reverseComplementMotif <- function(pwm) {
  rows <- rownames(pwm)
  cols <- colnames(pwm)
  Ns <- pwm["N", ]
  pwm <- pwm[4:1, length(cols):1]
  pwm <- rbind(pwm, Ns)
  rownames(pwm) <- rows
  colnames(pwm) <- cols
  return(pwm)
}
# ----------------------------------------------------------------------- #


# ----------------------------------------------------------------------- #
scoreSeqWindows <- function(ppm, seq) {
  ppm.width <- ncol(ppm)
  seq.len <- length(seq)
  diag.ind <- rep.int(ppm.width, seq.len - ppm.width)
  ranges <- vapply(c(0L, cumsum(diag.ind)),
                   function(x,
                            range = (1L + 0L:(ppm.width - 1L) * (ppm.width + 1L)))
                   {
                     x + range
                   },
                   integer(ppm.width))
  scores <- t(ppm[seq, ])[ranges]
  scores_rc <- t(reverseComplementMotif(ppm)[seq, ])[ranges]
  scores <- split(scores, ceiling(seq_along(scores)/ppm.width))
  scores_rc <- split(scores_rc, ceiling(seq_along(scores_rc)/ppm.width))
  res <- vapply(Map(function(x, y) {matrix(data = c(x, y), nrow = 2,
                                           byrow = TRUE, dimnames = list(c(1, 2)))},
                    scores, scores_rc),
                rowSums,
                numeric(2))
  return(res)
}
# ----------------------------------------------------------------------- #
# ======================================================================= #



# ======================================================================= #
# =================================== MAIN ============================== #
# ======================================================================= #
# --------------------
snps <- read.table("~/Documents/nathansonlab/tecac-manuscript/paintor-res-snps.txt", header = F)

# add two snps that had only 1 causal variant selected
snps <- c(snps$V1,  "rs55873183", "rs17336718")


print("running snps.from.rsid2")
variants <- snps.from.rsid2( rsid = snps, 
                            dbSNP = SNPlocs.Hsapiens.dbSNP142.GRCh37,
                            search.genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19)
print("done")

print("running motifbreakR2")
results <- motifbreakR2(snpList = variants, filterp = TRUE,
                       pwmList = hocomoco,
                       threshold = 1e-4,
                       method = "log",
                       bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                       BPPARAM = BiocParallel::bpparam())
print("done")


found.snps <- unique(names(results))
all.snps <- unique(names(variants))
missing <- all.snps[!(all.snps %in% found.snps)]


results2 <- motifbreakR2(snpList = variants[missing], filterp = TRUE,
                        pwmList = hocomoco,
                        threshold = 1e-4,
                        method = "log",
                        bkg = c(A=0.25, C=0.25, G=0.25, T=0.25),
                        BPPARAM = BiocParallel::bpparam())

# per snp, how many strong effects?
snps.res <- unique(results$SNP_id)
out <- data.frame(snps = snps.res, n.effect = 0)
for( i in snps.res)
{
  tmp <- results[ results$SNP_id == i,]
  out$n.effect[out$snps == i] <- sum(tmp$effect == "strong")
}

write.table(out, "motifbreakR-count.txt", row.names = F, col.names = T, quote = F)

ggplot(data = out, aes(x = n.effect)) + geom_histogram()

write.table(results, "motifbreakR-paintor-res.txt", row.names = F, col.names = T, quote = F, append = F, sep = ";")
print("plotting")
plotMB2(results = results[1:2], rsid = "rs1052053:C", effect = "strong")
print("done")
# ======================================================================= #
# ======================================================================= #
# ======================================================================= #