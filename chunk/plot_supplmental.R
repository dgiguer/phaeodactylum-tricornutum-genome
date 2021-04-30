# read in the files
d <- read.table("data/final_gc_content.txt", header = TRUE)
reads <- read.table("data/april_20_final_path_25kb.paf", sep = "\t")
gff <- read.table("data/final_pt_april_20.fasta.mod.EDTA.TEanno.gff3")
methyl <- read.table("data/modified_bases.5mC.bed")




# create a PDF output - this does 4 chromosomes per page
pdf("final_density.pdf", width = 12, height = 24)
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
    
    # generate first plot of gc content
    plot(chromosome_index$gc_content, type = "l", lwd = 0.5, col = rgb(0,0,0,0.2), xlab = "Chromosome position X 100", ylab = "", main = chromosome_name, ylim = c(0, 2), yaxt="n")
    axis(2, at = c(0, 0.25, 0.5, 0.75))
    mtext(text = "GC content (%)", side = 2, line = 2.5, at=0.5)
    
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
            arrows(x0=reads_subset$V8[i], y0=1, x1=reads_subset$V9[i], y1=1, lwd = 2, angle = 20, length = 0.1, col = rgb(156, 199, 255, 200, max = 255))
                
            } else if (reads_subset$V5[i] == "-") {
                
            arrows(x0=reads_subset$V9[i], y0=1, x1=reads_subset$V8[i], y1=1, lwd = 2, angle = 20, length = 0.1, col = rgb(245, 165, 86, 200, max = 255))
                
            }
        
        } else { 
            
            if (reads_subset$V5[i] == "+") {
        
            # these are all on the negative strand
            arrows(x0=reads_subset$V8[i], y0=2, x1=reads_subset$V9[i], y1=2, lwd = 2, angle = 20, length = 0.1, col = rgb(156, 199, 255, 200, max = 255))
                
            } else if (reads_subset$V5[i] == "-") {
                
            arrows(x0=reads_subset$V9[i], y0=2, x1=reads_subset$V8[i], y1=2, lwd = 2, angle = 20, length = 0.1, col = rgb(245, 165, 86, 200, max = 255))
                
            }
        } 
        
    }
    

    #testing() 
    ### plot the methylation 
    par(new=TRUE)
    #colour tap half light grey 
    
    methyl_subset <- methyl[which(methyl$V1 == chromosome_name),]
    #test_subset <- test[which(test$name == chromosome_name),]
    plot(y = c(-300, 500), x = c(0, max(reads_subset$V7)), xlab = "", axes=FALSE, ylab = "", type = "n")
    
    rect(xleft = -100000, xright= max(reads_subset$V9)+100000, ybottom = 0, ytop = 300, col = rgb(142,93,176, 20, max = 255), border=NA)
    
    lines(x=methyl_subset$V3, y=(methyl_subset$V11 * 3), col = rgb(142,93,176, 150, max = 255), xlim = c(0, max(reads_subset$V7)))
    
    axis(4, at = c(0, 300), labels = c(0, 100))
    mtext(text = "% of reads with 5mC", side = 4, line = 2.5, at=150)
    
    # try density of regions 
    par(new=TRUE)
    density <- density(te_subset$V5, adjust = 0.01)
    plot(density, col = rgb(0, 0, 0, 0.5), axes=FALSE, xlab = "", ylab="", main = "", ylim = c(- (max(density$y) / 0.3), max(density$y)), xlim = c(0, max(reads_subset$V7)))
    
}

plot(chromosome_index$gc_content, lwd = 0.5, col = rgb(0,0,0,0.2), xlab = "", ylab = "", main = "", ylim = c(0, 2), xaxt="n", yaxt="n", type = "n", bty="n")

legend("center", legend = "Supplementary Figure 2. Minimum overlapping tiling path of reads are shown in \n arrows, blue is a forward read,orange is a reverse read. GC content is shown in \n bottom third in grey. Percentage of methlayed bases is plotted in \n the middle third in purple. The density of LTR retrotranposons predicted by \n EDTA is plotted in the top third. Light blue boxes indicates start and stop of \n LTR retrotransposons, while light red shading indicates the window(s) with \n minimum GC content per genome, indicating putative centromeres.", bty="n", cex = 2)

dev.off()



