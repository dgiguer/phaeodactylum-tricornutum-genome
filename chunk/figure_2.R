R

# need to redefine this function from plotrix to accept vectors
smoothColorsNew <- function (..., alpha = NA) 
{
    # accept vector like ("red" "5" "blue" "5")
    if (is.vector(...)) {
            args <- strsplit(..., " ")
            # for every i %%2
            for (i in seq(args)) {
                if (i %% 2 == 0) {
                    args[[i]] <- as.integer(args[[i]])
                }
            }
    } else { 
        args <- list(...)
    }
    
    r <- g <- b <- NULL
    while (length(args) > 0) {
        if (!is.character(args[[1]])) 
            stop("Usage: smoothColors(\"color name\",[n|\"color name\"],...,\"color name\")")
        arglen <- length(args)
        if (arglen > 1) {
            if (is.numeric(args[[2]])) {
                lastarg <- 2
                while (is.numeric(args[[lastarg]])) {
                  lastarg <- lastarg + 1
                  if (lastarg > arglen) 
                    stop("bad argument list")
                }
                from <- col2rgb(args[[1]])
                too <- col2rgb(args[[lastarg]])
                n <- args[[2]] + 2
                r <- c(r, seq(from[1, ], too[1, ], length = n))
                i <- length(r)
                r <- r[-i]
                g <- c(g, seq(from[2, ], too[2, ], length = n))
                g <- g[-i]
                b <- c(b, seq(from[3, ], too[3, ], length = n))
                b <- b[-i]
                args <- args[-(1:(lastarg - 1))]
            }
            else {
                cc <- col2rgb(args[[1]])
                r <- c(r, cc[1, ])
                g <- c(g, cc[2, ])
                b <- c(b, cc[3, ])
                args <- args[-1]
            }
        }
        else {
            cc <- col2rgb(args[[1]])
            r <- c(r, cc[1, ])
            g <- c(g, cc[2, ])
            b <- c(b, cc[3, ])
            args <- args[-1]
        }
    }
    if (is.na(alpha)) 
        rgb(r, g, b, maxColorValue = 255)
    else rgb(r, g, b, alpha = alpha, maxColorValue = 255)
}
d <- read.table("data/50kb_coverage.regions.bed", sep = "\t")
pt <- read.table("data/pt_reference_table.txt", sep = "\t", header= TRUE)

g <- read.table("data/genome_alignments_sorted.paf")

# plot the chromosomes in order of length
d.length <- g[order(g$V7, decreasing=TRUE),]

d.length.large <- d.length[which(d.length$V11 > 10000),]

d.plot <- data.frame(chromosome=unique(d.length$V6), length=unique(d.length$V7), index = 1:25)
# use weird new tools to split column into new column
library(dplyr)
library(tidyr)
library(R.utils)
library(plotrix)

colours <- c("#ff701e", "#a9ddfc",
"#1eacff",
"#1e3cff")

d.new <- d %>% separate(V1, into=c("names", "number"), sep="_", remove=FALSE)
d.new$names <- NULL

d.new$colour <- "nothing"

for (i in seq(rownames(d.new))) {
    
        # if lower or between
        if (d.new$V4[i] < 75) {
                d.new$colour[i] <- as.character(colours[1])
        } else if (d.new$V4[i] >= 75 && d.new$V4[i] < 125) {
            d.new$colour[i] <- as.character(colours[2])
        } else if (d.new$V4[i] >= 125 && d.new$V4[i] < 175) {
            d.new$colour[i] <- as.character(colours[3])
        } else if (d.new$V4[i] >= 175) {
            d.new$colour[i] <- as.character(colours[4])
        } 
        
}


d.sub <- d.new

pdf("figs/figure_2.pdf", height = 8)
# generate blank plot 
plot(d.sub$number, d.sub$V2, col = d.sub$colour, axes = FALSE, xlab="", ylab = "", type = "n")

# loop through each chromosome 
for (i in seq(length(unique(d.sub$V1)))) {
    
    # subset to chromosome of interest
    subset <- d.sub[which(d.sub$V1 == unique(d$V1)[i]),]
    
    chromosome_colours <- insert(subset$colour, ats=2:length(subset$colour), values = 100)
    
    # plot discrete bars for each chromosome as x1=i, y1=0, x2=i+5, y2=length of chromosome
    # coverage plot
    # seperate chromosome
    if (i %%2 == 0) {
        
        # add telomere indicator (make it look like a "T" (bottom of chromosome)
        rect(ybottom=0, ytop=-25000, xleft=i,xright=i, col = rgb(0,0,0,0.5))
        rect(ybottom=-25000, ytop=-25000, xleft=i-0.1,xright=i+0.1, col = rgb(0,0,0,0.5))
        
        # add telomere indicator (make it look like a "T" (top of chromosome)
        rect(ybottom=max(subset$V3), ytop=max(subset$V3)+25000, xleft=i,xright=i, col = rgb(0,0,0,0.5))
        rect(ybottom=max(subset$V3)+25000, ytop=max(subset$V3)+25000, xleft=i-0.1,xright=i+0.1, col = rgb(0,0,0,0.5))
        
        gradient.rect((i-0.15),0, (i+0.15), max(subset$V3), col=smoothColorsNew(chromosome_colours), border="#1e3cff", gradient = "y", nslices=2)
    } else {
        
        # add telomere indicator (make it look like a "T" (bottom of chromosome)
        rect(ybottom=0, ytop=-25000, xleft=i,xright=i, col = rgb(0,0,0,0.5))
        rect(ybottom=-25000, ytop=-25000, xleft=i-0.1,xright=i+0.1, col = rgb(0,0,0,0.5))
        
        # add telomere indicator (make it look like a "T" (top of chromosome)
        rect(ybottom=max(subset$V3), ytop=max(subset$V3)+25000, xleft=i,xright=i, col = rgb(0,0,0,0.5))
        rect(ybottom=max(subset$V3)+25000, ytop=max(subset$V3)+25000, xleft=i-0.1,xright=i+0.1, col = rgb(0,0,0,0.5))
        
        gradient.rect((i-0.15),0, (i+0.15), max(subset$V3), col=smoothColorsNew(chromosome_colours), border="#1e3cff", gradient = "y", nslices=2)
    }
    

}

##### previous assembly

for (i in seq(length(unique(d.plot$chromosome)))) {
    
    current_chromosome <- d.plot$chromosome[i]
    
    # reference mapped
    subset <- d.length.large[which(d.length.large$V6 == current_chromosome),]
    
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
    
    # add colour boxes to show where each chromosome starts and stops
    # x axis is chromosome number, y axis is length
    counter <- 1
    for (v in seq(length(ref_chromosomes$start))) {
        if (counter %% 2 == 0) {
        rect(ybottom=ref_chromosomes$start[v], ytop=ref_chromosomes$end[v], xleft=d.plot$index[i]-0.25,xright=d.plot$index[i]-0.4, col = rgb(0, 0,0,0.3), border = "NA")
        } else {
        rect(ybottom=ref_chromosomes$start[v], ytop=ref_chromosomes$end[v], xleft=d.plot$index[i]-0.25,xright=d.plot$index[i]-0.4, col = rgb(0, 0,0,0.3), border = "NA")
        }
        counter <- counter + 1
    }
    
    # add text descriptors each chromosome
    for (v in seq(length(ref_chromosomes$spot))) {
        text(y=ref_chromosomes$spot[v], x=d.plot$index[i]-0.6, labels=ref_chromosomes$number[v], cex = 0.5, srt=90)
    }
}

axis(1, las =2, cex.axis = 0.5, at = 1:25, lwd = 0.5)
axis(4, cex.axis = 0.8, lwd = 0.5, at = seq(0, 2600000, by = 500000), labels = c("0", "0.5", "1.0", "1.5", "2.0", "2.5"))

dev.off()

pdf("figs/legend.pdf")
plot.new()
legend("topright", legend = c("50-75", "75-125", "125-175", "175+"), col = colours, pch = 15, cex = 1, bty = "n") 
dev.off()

# old colours for old chromosome: rgb(30, 172, 255, 175, maxColorValue = 255) and rgb(255, 112, 30, 175, maxColorValue = 255)

################################################################################

#  recreate chromosome 3
reads <- read.table("data/final_alignments.paf", sep = "\t", fill = NA)

chr3 <- reads[which(reads$V6 == "chromosome_3"),]
chr3_filtered <- chr3[which( ( (chr3$V4-chr3$V3) / chr3$V2) > 0.8),]

# make figure for chromosome three overlapping reads. 
current_chromosome <- "chromosome_3"

pdf("figs/chromosome_3.pdf")
# generate blank plot 
plot(d.plot$length[3], d.plot$index[3], xlab="", ylab = "", xlim = c(-1000, 2600000), axes=FALSE, type="n", main = "Overlapping reads for chromosome 3")

subset <- d.length.large[which(d.length.large$V6 == current_chromosome),]

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

# plot overlapping reads

# plot overlapping reads

# add colour boxes to show where each chromosome starts and stops
counter <- 1
for (v in seq(length(ref_chromosomes$start))) {
    if (counter %% 2 == 0) {
    rect(xleft=ref_chromosomes$start[v], xright=ref_chromosomes$end[v], ybottom=d.plot$index[3]-0.15,ytop=d.plot$index[3] -0.65, col = rgb(0,0,0,0.2))
    } else {
    rect(xleft=ref_chromosomes$start[v], xright=ref_chromosomes$end[v], ybottom=d.plot$index[3]-0.15,ytop=d.plot$index[3] -0.65, col = rgb(0,0,0,0.2))
    }
    counter <- counter + 1
}

# add text descriptors each chromosome
for (v in seq(length(ref_chromosomes$spot))) {
    text(x=ref_chromosomes$spot[v], y=d.plot$index[3] -0.85, labels=ref_chromosomes$number[v], cex = 1, srt=90)
}

#axis(1, cex.axis = 0.8, lwd = 0.5, at = seq(0, 2600000, by = 500000), labels = c("0", "0.5", "1.0", "1.5", "2.0", "2.5"))

reads_subset <- read.table("data/chromosome_3_path.paf")
# remove last read because it's redundant
reads_subset <- reads_subset[1:25,]

# plot overlapping reads 
for (i in 1:dim(reads_subset)[1]) {
    
    # alternative placement based on even or odd number
    # for even numbers, read on top, odd, place on bottom
    if (i %% 2 == 0) {
        
            if (reads_subset$V5[i] == "+") {
            
        # these are all on the negative strand
        arrows(x0=reads_subset$V8[i], y0=3.1, x1=reads_subset$V9[i], y1=3.1, lwd = 2, angle = 20, length = 0.1, col = rgb(156, 199, 255, 200, max = 255))
            
        } else if (reads_subset$V5[i] == "-") {
            
        arrows(x0=reads_subset$V9[i], y0=3.1, x1=reads_subset$V8[i], y1=3.1, lwd = 2, angle = 20, length = 0.1, col = rgb(245, 165, 86, 200, max = 255))
            
        }
    
    } else { 
        
        if (reads_subset$V5[i] == "+") {
    
        # these are all on the negative strand
        arrows(x0=reads_subset$V8[i], y0=3, x1=reads_subset$V9[i], y1=3, lwd = 2, angle = 20, length = 0.1, col = rgb(156, 199, 255, 200, max = 255))
            
        } else if (reads_subset$V5[i] == "-") {
            
        arrows(x0=reads_subset$V9[i], y0=3, x1=reads_subset$V8[i], y1=3, lwd = 2, angle = 20, length = 0.1, col = rgb(245, 165, 86, 200, max = 255))
            
        }
    }} 

# overlapping read 1: ae486c36-e71f-4599-bd01-f6edb622a81c
one <- chr3_filtered[which(chr3_filtered$V1 == "ae486c36-e71f-4599-bd01-f6edb622a81c"),]

arrows(x0=one$V9, y0=2.25, x1=one$V8, y1=2.25, lwd = 2, angle = 20, length = 0.1, col = "black")

# two = ecb2f85c-3420-4ef6-9834-7c5d98861a76
two <- chr3_filtered[which(chr3_filtered$V1 == "ecb2f85c-3420-4ef6-9834-7c5d98861a76"),]

arrows(x0=two$V8, y0=2.25, x1=two$V9, y1=2.25, lwd = 2, angle = 20, length = 0.1, col = "black")

# three = 63328598-9e4f-4ffb-9b2c-5dbf483a7689

two <- chr3_filtered[which(chr3_filtered$V1 == "63328598-9e4f-4ffb-9b2c-5dbf483a7689"),]

arrows(x0=two$V8, y0=2.25, x1=two$V9, y1=2.25, lwd = 2, angle = 20, length = 0.1, col = "black")

dev.off()

# combine the chromosome_3.pdf, legend.pdf, and figure_2.pdf in omnigraffle. 
```