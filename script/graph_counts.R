#!/usr/bin/env Rscript

Args        <- commandArgs()
dataFile    <- ifelse(is.na(Args[6]), "all.tsv",     Args[6])
samplesFile <- ifelse(is.na(Args[7]), "samples.txt", Args[7])
pdfFile     <- ifelse(is.na(Args[8]), "counts.pdf",  Args[8])

# Read data
if (grepl("csv$", dataFile)) {
    data <- read.csv(dataFile, header=TRUE)
} else if (grepl("tsv$", dataFile)) {
    data <- read.delim(dataFile, header=TRUE)
} else {
    stop(sprintf("Can't read %s file", dataFile))
}

# Read samples
samples <- read.table( samplesFile, header=TRUE, row.names=1 )

# Graph parameters
labels <- gsub(".normalised.count$", "",
               names(data)[grepl(".normalised.count$", names(data))])
colours <- as.numeric(samples[labels, "condition"])

pdf(pdfFile)

# Plot each region separately
for (i in 1:nrow(data)) {
    counts <- data[i, grepl(".normalised.count$", names(data)) ]
    par(mar=c(8.1, 4.1, 4.1, 2.1), xpd=TRUE)
    plot(as.numeric(counts), axes=FALSE, ann=FALSE, pch=21, bg=colours)
    axis(1, at=1:length(labels), lab=labels, las=2)
    axis(2)
    title(main=sprintf("%d:%d-%d %s\n%.2f",
                       data[i,"Chr"], data[i,"Region.start"],
                       data[i,"Region.end"], data[i,"Gene.name"],
                       data[i,"Adjusted.p.value"]))
    title(xlab="")
    title(ylab="Normalised Counts")
    legend("topright", inset=c(0, -0.1), levels(samples$condition), pch=21,
           pt.bg=colours)
}

graphics.off()
