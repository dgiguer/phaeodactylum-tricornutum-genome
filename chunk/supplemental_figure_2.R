# read in the files
d <- read.table("data/final_gc_content.txt", header = TRUE)
reads <- read.table("data/minimum_tiling_path_25kb.paf", sep = "\t")
gff <- read.table("data/final_pt_april_30_final.fasta.mod.EDTA.TEanno.gff3")
methyl <- read.table("data/modified_bases.5mC.bed")
cov <- read.table("data/1000_coverage.regions.bed", sep = "\t")

# for chromosome 3 only
alignments <- read.table("data/genome_alignments_sorted")
pt <- read.table("data/pt_reference_table.txt", sep = "\t", header= TRUE)

# create a PDF output - this does 4 chromosomes per page
pdf("figs/supplemental_figure_2.pdf", width = 12, height = 24)
par(mfrow=c(4, 1), mar=c(5.1,4.1,4.1,4.1))
    

# loop through all the chromosomes 
for (i in unique(d$chromosome)) {
    
    # get chromosome name
    chromosome_name <- i 

    # we don't care about circular elements here so skip them
    if (chromosome_name == "mitochondrian" | chromosome_name == "chloroplast") {
        next()
    }

    # get positions of current chromosome in loop in data files
    chromosome_index <- d[which(d$chromosome == chromosome_name),]

    # find minimum gc content bin
    min_gc_position <- which(chromosome_index$gc_content == min(chromosome_index$gc_content))

    par(mar=c(5.1,4.1,4.1,4.1))
    
    ############################################################################
    # plot gc content

    # generate first plot of gc content
    plot(chromosome_index$gc_content, type = "l", lwd = 0.5, col = rgb(0,0,0,0.2), xlab = "Chromosome position X 100", ylab = "", ylim = c(0.25, 2), yaxt="n")
    axis(2, at = c(0.25, 0.5, 0.75), col.axis = rgb(0,0,0,0.5), col.ticks = rgb(0,0,0,0.5))
    mtext(text = "GC content (%)", side = 2, line = 2.5, at=0.5, col = rgb(0,0,0,0.5))

    ############################################################################
    # plot coverage 
    cov_subset <- cov[which(cov$V1 == chromosome_name),]
    par(new=TRUE)
    maximum <- max(cov_subset$V4)
    plot(cov_subset$V4, ylim = c((-maximum*1.3), (maximum * 3.5)), col = rgb(0,0,0,0.6), type = "l", xlab = "", xaxt="n", yaxt="n", ylab="", main = "", lwd = 2)
    axis(4, at = c(0, round((maximum/2), 0), round(maximum, 0))) 
    mtext(text = "Coverage", side = 4, line = 2.5, at=maximum/2)
    abline(h=(maximum), lty = 2, col = rgb(0,0,0,0.2))
    abline(h=(maximum/2), lty = 2, col = rgb(0,0,0,0.2))
    abline(h=0, lty = 2, col = rgb(0,0,0,0.4))
    
    ############################################################################
    
    # plot reference chromosomes (only for chr 3)
    
    # plot the chromosomes in order of length
    alignments.length <- alignments[order(alignments$V7, decreasing=TRUE),]
    
    alignments.length.large <- alignments.length[which(alignments.length$V11 > 10000),]
    
    current_chromosome = i
    
    subset <- alignments.length.large[which(alignments.length.large$V6 == current_chromosome),]
    
    
    # get reference information to plot 
    ref_chromosomes <- data.frame(name=unique(subset$V1), spot=0, number=0, start=0, end=0)
    
    # if a chromosome, put number, otherwise put S for scaffold
    for (name in seq(length(ref_chromosomes$name))) {
        if (ref_chromosomes$name[name] %in% pt$RefSeq) {
            ref_chromosomes$number[name] <- pt$Name[grep(ref_chromosomes$name[name], pt$RefSeq)]
        } else {
            ref_chromosomes$number[name] <- "S"
        }
    
    }
    
    # now obtain the spot to place the chromosome number on the image. 
    for (w in (seq(length(ref_chromosomes$spot)))) {
        start <- min(subset[which(subset$V1 == ref_chromosomes$name[w]),]$V8)
        end <- max(subset[which(subset$V1 == ref_chromosomes$name[w]),]$V9)
        ref_chromosomes$spot[w] <- ( (start + end) / 2)
        ref_chromosomes$start[w] <- start
        ref_chromosomes$end[w] <- end
    }
    par(new=TRUE)
    plot(chromosome_index$gc_content, xlim=c(0, max(subset$V7)), type = "n", xlab = "", ylab = "", main = "", ylim = c((-maximum*1.3), (maximum * 3.5)), yaxt="n", xaxt="n")
    
    # add colour boxes to show where each chromosome starts and stops
    counter <- 1
    for (v in seq(length(ref_chromosomes$start))) {
        if (counter %% 2 == 0) {
        rect(xleft=ref_chromosomes$start[v], xright=ref_chromosomes$end[v], ybottom=(maximum *1.05),ytop=(maximum*1.15), col = rgb(0,0,0,0.2), border="NA")
        } else {
        rect(xleft=ref_chromosomes$start[v], xright=ref_chromosomes$end[v], ybottom=(maximum *1.05),ytop=(maximum*1.15), col = rgb(0,0,0,0.2), border="NA")
        }
        counter <- counter + 1
    }
    
    
    # add text descriptors each chromosome
    for (v in seq(length(ref_chromosomes$spot))) {
        text(x=ref_chromosomes$spot[v], y=maximum*1.20, labels=ref_chromosomes$number[v], cex = 0.5)
    }


    ############################################################################
    # plot TEs
    ###### replot the without data just to get axes back to correct ones
    par(new=TRUE)
    plot(chromosome_index$gc_content, lwd = 0.5, col = rgb(0,0,0,0.2), xlab = "", ylab = "",ylim = c(0.25, 2), yaxt="n", type = "n")

    # subset the TEs to get only this chromosome
    te_subset <- gff[which(gff$V1 == chromosome_name),]
    te_subset <- te_subset[which(te_subset$V3 == "LTR_retrotransposon" | te_subset$V3 == "Copia_LTR_retrotransposon" | te_subset$V3 == "Gypsy_LTR_retrotransposon"),]

    # plot the regions of each TE by looping through them all
    for (i in seq(nrow(te_subset))) {

        # plot for each one
        # divide by 1000 because we calculated coverage in 1000 base windows
        rect(xleft = te_subset$V4[i] / 100, xright = te_subset$V5[i] /100, ybottom = -0.2, ytop = 2.2, border=NA, col = rgb(156, 199, 255, 50, max = 255))

    }

    # put the minimum gc position in between lines for highlighting
    # abline(v=min_gc_position - 50, lty = 2, col = rgb(1,0,0,0.5))
    # abline(v=min_gc_position + 50, lty = 2, col = rgb(1,0,0,0.5))
    # plot GC content
    rect(xleft = min_gc_position - 50, xright= min_gc_position + 50, ybottom = -0.2, ytop = 2.2, col = rgb(1,0,0,0.1), border=NA)

    # get overlapping reads for current chromosome
    reads_subset <- reads[which(reads$V6 == chromosome_name),]

    # plot
    par(new=TRUE)
    plot(y = c(0, 20), x = c(0, max(reads_subset$V9)), xlab = "", axes=FALSE, ylab = "", type = "n")

    # x0, y0, x1, y1

    # plot overlapping reads for each chromosome
    for (i in 1:dim(reads_subset)[1]) {
        
        # alternative placement based on even or odd number
        # for even numbers, read on top, odd, place on bottom
        if (i %% 2 == 0) {
            
                if (reads_subset$V5[i] == "+") {
                
            # these are all on the negative strand
            arrows(x0=reads_subset$V8[i], y0=0.15, x1=reads_subset$V9[i], y1=0.15, lwd = 2, angle = 20, length = 0.1, col = rgb(156, 199, 255, 200, max = 255))
                
            } else if (reads_subset$V5[i] == "-") {
                
            arrows(x0=reads_subset$V9[i], y0=0.15, x1=reads_subset$V8[i], y1=0.15, lwd = 2, angle = 20, length = 0.1, col = rgb(245, 165, 86, 200, max = 255))
                
            }
        
        } else { 
            
            if (reads_subset$V5[i] == "+") {
        
            # these are all on the negative strand
            arrows(x0=reads_subset$V8[i], y0=1.15, x1=reads_subset$V9[i], y1=1.15, lwd = 2, angle = 20, length = 0.1, col = rgb(156, 199, 255, 200, max = 255))
                
            } else if (reads_subset$V5[i] == "-") {
                
            arrows(x0=reads_subset$V9[i], y0=1.15, x1=reads_subset$V8[i], y1=1.15, lwd = 2, angle = 20, length = 0.1, col = rgb(245, 165, 86, 200, max = 255))
                
            }
        } 
        
    }


    #testing() 
    ### plot the methylation c
    par(new=TRUE)
    #colour tap half light grey 

    ###########################################################################
    
    # methylation

    methyl_subset <- methyl[which(methyl$V1 == chromosome_name),]
    #test_subset <- test[which(test$name == chromosome_name),]
    plot(y = c(-725, 625), x = c(0, max(reads_subset$V7)), xlab = "", axes=FALSE, ylab = "", type = "n")

    rect(xleft = -100000, xright= max(reads_subset$V9)+100000, ybottom = 0, ytop = 300, col = rgb(142,93,176, 20, max = 255), border=NA)

    lines(x=methyl_subset$V3, y=(methyl_subset$V11 * 3), col = rgb(142,93,176, 150, max = 255), xlim = c(0, max(reads_subset$V7)))

    axis(2, at = c(0, 300), labels = c(0, 100), col.ticks=rgb(142,93,176, 225, max = 255), col.axis=rgb(142,93,176, 225, max = 255))
    mtext(text = "% of reads 5mC", side = 2, line = 2.5, at=150, col = rgb(142,93,176, 225, max = 255))

    # try density of regions 
    par(new=TRUE)
    density <- density(te_subset$V5, adjust = 0.01)
    plot(density, col = rgb(0, 0, 0, 0.5), axes=FALSE, xlab = "", ylab="", main = "", ylim = c(- (max(density$y) / 0.3), max(density$y)), xlim = c(0, max(reads_subset$V7)))
    mtext(text = "Density", side = 4, line = 2.5, at= (max(density$y) / 2))
    
}

dev.off()