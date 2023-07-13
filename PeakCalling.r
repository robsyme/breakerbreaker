#!/usr/bin/env Rscript

## ----dependencies, warning=FALSE, message=FALSE--------------------------
#  if (!require("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")

#  BiocManager::install("Rsamtools")
#  BiocManager::install("GenomicRanges")
#  BiocManager::install("GenomicAlignments")
#  install.packages("PeakSegDisk")
library(data.table)
library(Rsamtools)
library(tidyverse)
library(GenomicAlignments)
library(GenomicRanges)
library(PeakSegDisk)
library(writexl)

## ----bamPath-------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
print("args passed in")
BAMs <- sapply(args, BamFile)
print(BAMs)

#Binning genome
gnm <- read_tsv("/data/alistairh/projects/Genome_Comparison/data/Carv2021.genome", colnames(FALSE)) %>%  rename("names" = "X1", "length" = "X2") %>%
  mutate(
    start = 0
  ) %>%
  relocate(start, .after = names)

chr <- c("Ca1_v2.0", "Ca2_v2.0", "Ca3_v2.0", "Ca4_v2.0", "Ca5_v2.0", "Ca6_v2.0", "Ca7_v2.0", "Ca8_v2.0")

#Pileup Parameters
p_params<- PileupParam(max_depth=10e4, min_base_quality= 10, min_mapq= 1,
                       min_nucleotide_depth=0, min_minor_allele_depth= 0,
                       distinguish_strands= TRUE, distinguish_nucleotides= FALSE,
                       ignore_query_Ns= FALSE, include_deletions= TRUE, include_insertions=FALSE)

ProcessBam <- function(bamFile, gnm, chr, p_params) {
  message("Beginning processing of BAM file ")
  print(bamFile)
  sample_id <- gsub(".bam$", "", basename(bamFile))
  message("Processing sample: ", sample_id)
  #Create a list to store the pileup results
  peakResults = list()
  
  #re-write for loop into sapply function for each chromosome with list output
  peakResults <- sapply(chr, PeakRegions, gnm = gnm, p_params = p_params, bamFile = bamFile, simplify = FALSE) 
   
  INV_Results <- lapply(peakResults, `[[`, "INV") %>%
    list() %>%
    flatten() %>%
    reduce(full_join)

  RT_Results <- lapply(peakResults, `[[`, "RT") %>%
    list() %>%
    flatten() %>%
    reduce(full_join)

  message("Writing Peak File to disk ")
  write_csv(INV_Results, path = paste0(sample_id, "_INV_peaks.csv"))
  write_csv(RT_Results, path = paste0(sample_id, "_RT_peaks.csv"))
  return(NULL)
}

# lapply(peakResults, function(x) write.table( data.frame(x), paste0(sample_id, "_peaks.csv"[])  , append= T, sep=',' ))



#rewrite PileUp function to return a list of data frames
PeakRegions <- function(chromosome, gnm,  p_params, bamFile) {
  print("In PeakRegions function")
  genome.stats <- gnm %>% filter(names == chromosome)
  print(genome.stats)
  
  Ca.genome <- GRanges(seqnames = chromosome, ranges = IRanges(start = c(genome.stats$start), end = c(genome.stats$length)))
  
  #Parameters for paired reads
  params <- ScanBamParam(which = Ca.genome, what = scanBamWhat(), flag = Rsamtools::scanBamFlag(isPaired = TRUE))
  #Parameters for discordant reads
  params_dis <- ScanBamParam(which = Ca.genome, what = scanBamWhat(),flag = Rsamtools::scanBamFlag(isPaired=TRUE,isProperPair=FALSE,isDuplicate=FALSE,isNotPassingQualityControls=FALSE,isUnmappedQuery=FALSE,hasUnmappedMate=FALSE))
  
  #pileup of all properly paired reads and discordant reads
  discordant_bins <- pileup(bamFile, scanBamParam = params_dis, pileupParam = p_params) 
  
  message("Discordant bins created")
  
  #Extract the forward strand  discordant read coverage
  discordant_fwd <- discordant_bins %>%
    filter(strand == "+") %>%
    select(c("seqnames", "pos", "count")) %>%
    rename("chromEnd" = "pos", "chrom" = "seqnames") %>%
    mutate(chromStart = as.integer(chromEnd-1)) %>%
    relocate(chromStart, .after = chrom)
  
  #Extract the reverse strand discordant read coverage
  discordant_rev <- discordant_bins %>%
    filter(strand == "-") %>%
    select(c("seqnames", "pos", "count")) %>%
    rename("chromEnd" = "pos", "chrom" = "seqnames") %>%
    mutate(chromStart = as.integer(chromEnd-1)) %>%
    relocate(chromStart, .after = chrom)
  
  rm(discordant_bins)

    all.pos <- data.frame(chrom  = chromosome,
                        chromStart = c(genome.stats$start:(genome.stats$length-1)),
                        chromEnd = ((genome.stats$start +1): (genome.stats$length)))
  
  #Join the discordant read coverage to the genome positions
  dis_fwd <- all.pos %>%
    full_join(discordant_fwd) %>%
    replace(is.na(.),0)

  #Remove the rows with zero coverage for the forward strand

  dis_countRle <- rle(dis_fwd$count)
  dis_rleRows <- cumsum(dis_countRle$lengths)

  final_dis_fwd <- data.frame(
    chrom = dis_fwd$chrom[dis_rleRows],
    chromStart = as.integer(c(0, dis_rleRows) %>% head(-1)),
    chromEnd = as.integer(dis_fwd$chromEnd[dis_rleRows]),
    count = dis_countRle$values
  )

  #Remove the rows with zero coverage for the reverse strand

  dis_rev <- all.pos %>%
    full_join(discordant_rev) %>%
    replace(is.na(.),0)

  disrev_countRle <- rle(dis_rev$count)
  disrev_rleRows <- cumsum(disrev_countRle$lengths)

  final_dis_rev <- data.frame(
    chrom = dis_rev$chrom[disrev_rleRows],
    chromStart = as.integer(c(0, disrev_rleRows) %>% head(-1)),
    chromEnd = as.integer(dis_rev$chromEnd[disrev_rleRows]),
    count = disrev_countRle$values
  )

  message("Starting peak detection")
  
  #Detect the peaks in the discordant read coverage
  peaks<- PeakSegDisk::PeakSegFPOP_df(final_dis_fwd, 1000, base.dir = "/data/alistairh/projects/Genome_Comparison/data/peaks/temp")
  peaks_rev<- PeakSegDisk::PeakSegFPOP_df(final_dis_rev, 1000, base.dir = "/data/alistairh/projects/Genome_Comparison/data/peaks/temp")
  
  PeaksRev <- peaks_rev$segments %>%
  mutate(strand = "-")

  PeaksFwd <- peaks$segments %>%
  mutate(strand = "+")

  rm(peaks)
  rm(peaks_rev)

  message("Finished peak detection")
  
  peak_regions <- data.frame(PeaksRev) %>%
    full_join(PeaksFwd) %>%
    filter(status == "peak") %>%
    select(chromStart, chromEnd, strand) 


  message("peak_regions created")
  
  rm(PeaksRev)
  rm(PeaksFwd)
  
  #Extract the peak starts and ends
  peak_starts <- peak_regions[[1]]
  peak_ends <- peak_regions[[2]]
  peak_strand <- peak_regions[[3]]
  
  #Create a GRanges object of the peak regions
  gr.peaks <- GRanges(seqnames = chromosome,
                      ranges = IRanges(start= c(peak_starts), end = c(peak_ends)), strand = peak_strand)
  
  #ScanBamParam for discordant reads in the peak regions
  params_peak <- ScanBamParam(which = gr.peaks, what = scanBamWhat(),flag = Rsamtools::scanBamFlag(isPaired=TRUE,isProperPair=FALSE,isDuplicate=FALSE,isNotPassingQualityControls=FALSE,isUnmappedQuery=FALSE,hasUnmappedMate=FALSE))
  
  message("creating peak.df...")

rm(peak_regions)

  #Create a data frame of the discordant reads in the peak regions
  peak.scan <- scanBam(bamFile, param = params_peak)

  message("peak.scan created")

  peak.df <- rbindlist(lapply(peak.scan, as.data.frame), idcol = "rn")

  #peak.df <- dplyr::bind_rows(lapply(peak.scan, as.data.frame), .id = 'rn') # nolint
  #re-write line above to not overload memory in bash script
  
  message("peak.df created")

  #Create a data frame of the discordant reads in the peak regions
  peak.split <- peak.df %>%
    group_split(rn) 

    message("peak.split created")

  #Clean up memory by removing the peak.scan object
  rm(peak.scan)

  message("peak.scan removed")

  
  #Detect the proportion of reads with adjacent mate pairs for each peak regioN
  message("Creating results...")
  suppressWarnings({
  results <- peak.split %>% 
    purrr::map(function(peak.df) {
      peak.df %>% 
        dplyr::rowwise() %>% 
        dplyr::summarise(
          rn,
          qname,
          matched.frac = (sum(mrnm == peak.df$mrnm & abs(mpos - peak.df$mpos) < 1000) - 1) / (n() - 1)
        ) %>% 
        dplyr::summarise(
          rn = rn[1],
          peak.metric.mean = mean(matched.frac)
        )
    }) %>% 
    dplyr::bind_rows()
})
  
  message("results created")
  #Filter peak regions for those proprtion of adj. mate pairs above 0.5
  peaks_above <- results %>%
    filter(peak.metric.mean > 0.50) 
  
  rn_above <- as.vector(peaks_above$rn)
  
  message("Filtered peak regions above 50%")
  #Peak regions for intrachromosomal rearrangements greater than 1Mb
  peaks_INV <- peak.df %>%
    filter(rn %in% rn_above) %>%
    filter(rname == mrnm) %>%
    filter(abs(pos - mpos) > 1e6) %>%
    group_by(rn) %>%
    summarise(n = n(), query.mean = mean(qwidth)) %>%
    filter(n > 9) %>%
    separate(col = rn, into = c("chr", "range"), sep = "\\:") %>%
    separate(col = range, into = c("range.start", "range.end"), sep = "\\-") %>%
    mutate(
      range.start = as.numeric(range.start),
      range.end = as.numeric(range.end),
      range.size = range.end - range.start,
    ) 

    # peaks_INV_start <- as.vector(as.numeric(peaks_INV$range.start))  
    # peaks_INV_end <- as.vector(as.numeric(peaks_INV$range.end))

    # #GRange for large inversions (greater than 1Mb)  
    # peaks_INV_GRange <- GRanges(seqnames = chromosome, ranges = IRanges(start = peaks_INV_start, end = peaks_INV_end))

    # #ScanBamParam for large inversions (greater than 1Mb)
    # params_INV <- ScanBamParam(which = peaks_INV_GRange, what = scanBamWhat(),flag = Rsamtools::scanBamFlag(isPaired=TRUE,isProperPair=FALSE,isDuplicate=FALSE,isNotPassingQualityControls=FALSE,isUnmappedQuery=FALSE,hasUnmappedMate=FALSE))

    # #PileupParam for large inversions (greater than 1Mb)
    # p_params<- PileupParam(max_depth=10e4, min_base_quality= 10, min_mapq= 1, min_nucleotide_depth=0, min_minor_allele_depth= 0,
    # distinguish_strands= TRUE, distinguish_nucleotides= FALSE,ignore_query_Ns= FALSE, include_deletions= TRUE, include_insertions= FALSE)

    # #Pileup for large inversions (greater than 1Mb)
    # peaks_cov_INV <- pileup(bamFile, scanBamParam = params_INV, pileupParam = p_params)
  
  peaks_RT <-  peak.df %>%
    filter(rn %in% rn_above) %>%
    filter(rname != mrnm) %>%
    filter(mrnm %in% chr) %>%
    group_by(rn) %>%
    summarise(n = n(), query.mean = mean(qwidth)) %>%
    filter(n > 9) %>%
    separate(col = rn, into = c("chr", "range"), sep = "\\:") %>%
    separate(col = range, into = c("range.start", "range.end"), sep = "\\-") %>%
    mutate(
      range.start = as.numeric(range.start),
      range.end = as.numeric(range.end),
      range.size = range.end - range.start,
    ) 

    # #Peak regions starts and ends
    # peaks_RT_start <- as.vector(as.numeric(peaks_RT$range.start))
    # peaks_RT_end <- as.vector(as.numeric(peaks_RT$range.end))

    # #GRange for reciprocal traslocations
    # peaks_RT_GRange <- GRanges(seqnames = chromosome, ranges = IRanges(start = peaks_RT_start, end = peaks_RT_end))

    # #ScanBamParam for reciprocal traslocations
    # params_RT <- ScanBamParam(which = peaks_RT_GRange, what = scanBamWhat(),flag = Rsamtools::scanBamFlag(isPaired=TRUE,isProperPair=FALSE,isDuplicate=FALSE,isNotPassingQualityControls=FALSE,isUnmappedQuery=FALSE,hasUnmappedMate=FALSE))

    # #PileupParam for reciprocal traslocations
    # p_params<- PileupParam(max_depth=10e4, min_base_quality= 10, min_mapq= 1, min_nucleotide_depth=0, min_minor_allele_depth= 0,
    #     distinguish_strands= TRUE, distinguish_nucleotides= FALSE,ignore_query_Ns= FALSE, include_deletions= TRUE, include_insertions= FALSE)

    # #Pileup for reciprocal translocations
    # peaks_cov_RT <- pileup(bamFile, scanBamParam = params_RT, pileupParam = p_params)
  
    outputs <- list(RT = peaks_RT, INV = peaks_INV)
    # outputs_cov <- list(RT = peaks_cov_RT, INV = peaks_cov_INV)
  
    return(outputs)
 
}
print("Calling sapply...")

sapply(BAMs, ProcessBam, gnm = gnm, chr = chr, p_params = p_params)
