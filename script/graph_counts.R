#!/usr/bin/env Rscript

library(ggplot2)
library(reshape2)
suppressPackageStartupMessages(library(dplyr))

Args        <- commandArgs()
dataFile    <- ifelse(is.na(Args[6]), "all.tsv",     Args[6])
samplesFile <- ifelse(is.na(Args[7]), "samples.txt", Args[7])
pdfFile     <- ifelse(is.na(Args[8]), "counts.pdf",  Args[8])
plotStyle   <- ifelse(is.na(Args[9]), "default",     Args[9])

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

if (plotStyle == "violin") {
    countData$id <- row.names(countData)
    counts <- melt(countData, id.vars="id", variable.name="condition",
                   value.name="count")
    countData$name <- sprintf("%s:%d-%d %s\n%.2f",
                              data[,"Chr"],
                              data[,"Region.start"],
                              data[,"Region.end"],
                              data[,"Gene.name"],
                              data[,"Adjusted.p.value"])
    levels(counts$condition) <- samples$condition
    for (i in 1:nrow(data)) {
        print(ggplot(counts[counts$id == i,],
               aes(x=condition, y=count, color=condition)) +
            geom_violin() + geom_boxplot(width=0.1, outlier.shape=NA) +
            labs(x="", y="Normalised Counts", title=countData$name[i]) +
            theme_bw() +
            theme(legend.position='none',
                  plot.title = element_text(hjust = 0.5, face="bold")))
    }
} else {
    for (i in 1:nrow(data)) {
        par(mar=c(8.1, 4.1, 4.1, 2.1), xpd=TRUE)
        plot(as.numeric(countData[i,]), axes=FALSE, ann=FALSE, pch=21,
             bg=colours)
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
}

graphics.off()
