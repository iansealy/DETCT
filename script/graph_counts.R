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
    colnames(data)[1] = "Chr" # X.Chr -> Chr
} else {
    stop(sprintf("Can't read %s file", dataFile))
}

# Read samples
samples <- read.table( samplesFile, header=TRUE, row.names=1 )

# Get counts
countData <- data[,grepl(".normalised.count$", names(data))]
names(countData) <- gsub(".normalised.count$", "", names(countData))

# Subset and reorder count data
countData <- countData[, row.names(samples)]

# Graph parameters
colours <- as.numeric(samples$condition)

pdf(pdfFile)

# Plot each region separately
for (i in 1:nrow(data)) {
    par(mar=c(8.1, 4.1, 4.1, 2.1), xpd=TRUE)
    plot(as.numeric(countData[i,]), axes=FALSE, ann=FALSE, pch=21, bg=colours)
    axis(1, at=1:ncol(countData), lab=names(countData), las=2, cex.axis=0.5)
    axis(2)
    title(main=sprintf("%s:%d-%d %s\n%.2f",
                       data[i,"Chr"], data[i,"Region.start"],
                       data[i,"Region.end"], data[i,"Gene.name"],
                       data[i,"Adjusted.p.value"]))
    title(xlab="")
    title(ylab="Normalised Counts")
    legend("topright", inset=c(0, -0.1), levels(samples$condition), pch=21,
           pt.bg=1:length(levels(samples$condition)))
}

graphics.off()
