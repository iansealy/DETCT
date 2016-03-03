#!/usr/bin/env Rscript

suppressWarnings(library(tcltk))
suppressPackageStartupMessages(library(Mfuzz))

Args           <- commandArgs()
dataFile       <- ifelse(is.na(Args[6]),  "all.tsv",     Args[6])
samplesFile    <- ifelse(is.na(Args[7]),  "samples.txt", Args[7])
outputBase     <- ifelse(is.na(Args[8]),  "all",         Args[8])
numClusters    <- ifelse(is.na(Args[9]),  100,           as.integer(Args[9]))
alphaThreshold <- ifelse(is.na(Args[10]), 0.7,           as.numeric(Args[10]))

# Read data
if (grepl("csv$", dataFile)) {
    data <- read.csv(dataFile, header=TRUE, check.names=FALSE)
} else if (grepl("tsv$", dataFile)) {
    data <- read.delim(dataFile, header=TRUE, check.names=FALSE)
} else {
    stop(sprintf("Can't read %s file", dataFile))
}

# Read samples
samples <- read.table( samplesFile, header=TRUE, row.names=1 )

# Ensure chromosome is a factor, even if all numeric
data[,1] <- factor(data[,1])

# Get counts
countData <- data[,grepl(" normalised count$", names(data))]
names(countData) <- gsub(" normalised count$", "", names(countData))

# Get median counts
medianData <- matrix(nrow=nrow(countData),
                     ncol=length(levels(samples$condition)),
                     dimnames=list(rownames(countData),
                                   levels(samples$condition)))
for (condition in levels(samples$condition)) {
    medianData[,condition] <-
        apply(countData[,samples$condition == condition, drop=FALSE], 1,
              median)
}

# Standardise data
eset <- ExpressionSet(assayData=as.matrix(medianData))
eset <- filter.std(eset, min.std=0, visu=FALSE)
eset.s <- standardise(eset)

# Get fuzzy clusters
cl <- mfuzz(eset.s, c=numClusters, m=mestimate(eset.s))
pdf(paste0(outputBase, '-', numClusters, '-', alphaThreshold, '-mfuzz.pdf'))
mfuzz.plot2(eset.s, cl=cl, mfrow=c(4,4), x11=FALSE, centre=TRUE,
           time.labels=colnames(medianData))
graphics.off()
pdf(paste0(outputBase, '-', numClusters, '-mfuzz.pdf'))
mfuzz.plot2(eset.s, cl=cl, mfrow=c(4,4), x11=FALSE, centre=TRUE,
           time.labels=colnames(medianData), min.mem=alphaThreshold)
graphics.off()

# Output fuzzy clusters
labels <- gsub(".normalised.count$", "",
               names(data)[grepl(".normalised.count$", names(data))])
colours <- as.numeric(samples$condition)
acore.list <- acore(eset.s, cl=cl, min.acore=alphaThreshold)
for (i in 1:numClusters) {
    data.subset <- data[c(as.character(acore.list[[i]][,1])),]
    if (grepl("csv$", dataFile)) {
        write.csv(data.subset, file=paste0(outputBase, '-', numClusters, '-',
                                           alphaThreshold, '-', i, '.csv'),
                  row.names=FALSE)
    } else {
        write.table(data.subset, file=paste0(outputBase, '-', numClusters, '-',
                                             alphaThreshold, '-', i, '.tsv'),
                    quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
    }
    # Plot counts
    pdf(paste0(outputBase, '-', numClusters, '-', alphaThreshold, '-', i,
               '-counts.pdf'))
    for (i in 1:nrow(data.subset)) {
        counts <- data.subset[i, grepl(".normalised.count$",
                                       names(data.subset)) ]
        par(mar=c(8.1, 4.1, 4.1, 2.1), xpd=TRUE)
        plot(as.numeric(counts), axes=FALSE, ann=FALSE, pch=21, bg=colours)
        axis(1, at=1:length(labels), lab=labels, las=2, cex.axis=0.5)
        axis(2)
        title(main=sprintf("%s:%d-%d %s\n%.2f",
                           data.subset[i,"Chr"],
                           data.subset[i,"Region.start"],
                           data.subset[i,"Region.end"],
                           data.subset[i,"Gene.name"],
                           data.subset[i,"Adjusted.p.value"]))
        title(xlab="")
        title(ylab="Normalised Counts")
        legend("topright", inset=c(0, -0.1), levels(samples$condition), pch=21,
               pt.bg=1:length(levels(samples$condition)))
    }
    graphics.off()
}
