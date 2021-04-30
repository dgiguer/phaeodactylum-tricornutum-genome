#!/usr/bin/env Rscript
# daniel giguere

# ensure on output directory you have trailing backslash
# make sure the directory is already made as well
# genome.fasta may contain multiple sequences 
# usage :
# Rscript ./circlize_gc_information.R -i genome.fasta -o contig_name/contig_name-gc-information.txt -w 1000

suppressPackageStartupMessages(library(optparse))

library(optparse)
library(seqinr)
library(data.table)

option_list <- list(
    make_option(c("-i","--fasta"), type="character", default=NULL,
            help="fasta file of genome",
            dest="genome_filename"),
    make_option(c("-w","--window"), type="integer", default=1000,
            help="window length to calculate GC by",
            dest="window"),
    make_option(c("-o","--output"), type="character", default=NULL,
            help="path output directory [default = ./]",
            dest="output_directory")
            )

parser <- OptionParser(usage = "%prog -i genome.fasta -o contig_name/contig_name-gc-information.txt [options]",option_list=option_list)

opt = parse_args(parser)

fasta <- read.fasta(file = opt$genome_filename)  

window <- opt$window

# if there are multiple sequences in the file, you need to evaluate GC for each one and report it on a per sequence basis. 

fasta.length <- list()
starts <- list()
gc.content <- list()

for (i in seq(fasta)) {
    
    # get length of each fasta file
    fasta.length[i] <- length(fasta[[i]])   
    
    # to change window, change number to desired window size. 
    starts[[i]] <- seq(1, (fasta.length[[i]]), by = window)
    
    # get list of vector of GC content for 1000 base windows
    gc.content[[i]] <- numeric(length(starts[[i]]))
    
    # iterate through each window of the current fasta (i), and calculate GC content
    for (k in seq(length(starts[[i]]))) {

            # last chunk will need special processing

            if (k == length(starts[[i]])) {
              # only get sequence to the end of the fasta, otherwise
              # NAs will mess everything up
              # look in sequence i for window k
              chunk <- fasta[[i]][starts[[i]][k]:fasta.length[[i]]]
              chunkGC <- GC(chunk)
              # add GC content for window to gc vector
              # change the value of list *i* with value *k*
              gc.content[[i]][k] <- chunkGC
            } else {
              # chunks are 0-100,
              chunk <- fasta[[i]][starts[[i]][k]:(starts[[i]][k]+(window -1))]
              chunkGC <- GC(chunk)
              # add GC content for window to gc vector
              gc.content[[i]][k] <- chunkGC
            }
    }
}

# GC Skew = (G - C)/(G + C)
# cumulative gc skew will be summed, should start at the origin.
# fasta sequence needs to imported as a string in this case.
# re-read the fasta file as a string.
fasta.string <- read.fasta(file = opt$genome_filename, as.string=TRUE)  

# to test
# fasta.string <- read.fasta(file = "data/pt_polished.fa", as.string=TRUE)  

# reset variables
fasta.length <- list()
gc.skew <- list()
gc.culm <- list()
starts <- list()

# define function to be used
# taken from https://stat.ethz.ch/pipermail/r-help/2007-November/147069.html
gcskew <- function(x) {
   if (!is.character(x) || length(x) > 1)
   stop("single string expected")
   tmp <- tolower(s2c(x))
   nC <- sum(tmp == "c")
   nG <- sum(tmp == "g")
   if (nC + nG == 0) return(NA)
   return(100 * (nC - nG)/(nC + nG))
}

for (i in seq(fasta.string)) {
    
    # get length of each fasta file
    fasta.length[[i]] <- nchar(fasta.string[[i]])   
    
    # to change window, change number to desired window size. 
    starts[[i]] <- seq(1, (fasta.length[[i]]), by = window)
    
    # get vector of GC content for 100 base windows
    # create vector the length of the starts
    gc.skew[[i]] <- numeric(length(starts[[i]]))

    # for every starting window
    for (k in seq_len(length(starts[[i]]))) {
      # iterate through the entire genome
      # last iteration needs special treatment.
      if (k == length(starts[[i]])) {
        # only go to last base
        gc.skew[[i]][k] <- gcskew(substr(fasta.string[[i]], starts[[i]][k], fasta.length[[i]]))
      } else {
        gc.skew[[i]][k] <- gcskew(substr(fasta.string[[i]], starts[[i]][k], starts[[i]][k] + (window - 1 ) ))
      }
    }
    
    # calculate cumulative gc skew for each window.
    # basically the equation is going to be i + i-1.

    gc.culm[[i]] <- numeric(length(starts[[i]]))

    for (k in seq(gc.skew[[i]])) {

      # first one needs to be treated specially
      if (k == 1) {
        # only take first column
        gc.culm[[i]][k] <- gc.skew[[i]][k]
      } else {
        # this works for last one because it looks back, not ahead.
        gc.culm[[i]][k] <- sum(gc.skew[[i]][1:k])
      }
    }
    
}

# final long file needs to be something parsable:
# chromosome_name\tregion\tgc.content\tgc.skew\tgc_culm 
# make a list of dataframes then basically concatenate them all? 
final_gc_info <- list()

for (i in seq(fasta)) {
    
    final_gc_info[[i]] <- data.frame(chromosome = c(rep(attributes(fasta[[i]])$name, length(gc.content[[i]]) )),   gc_content = gc.content[[i]], gc_skew = gc.skew[[i]], gc_culm = gc.culm[[i]])
    
}

library(data.table)

final_data_frame <- rbindlist(final_gc_info)

write.table(final_data_frame, file = opt$output_directory, col.names=TRUE, row.names=FALSE, quote = FALSE)
