#!/usr/bin/env Rscript

suppressWarnings(library(tcltk))
suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(pheatmap))

Args               <- commandArgs()
dataFile           <- ifelse(is.na(Args[6]), "sig.tsv",     Args[6])
samplesFile        <- ifelse(is.na(Args[7]), "samples.txt", Args[7])
pngFile            <- ifelse(is.na(Args[8]), "heatmap.png", Args[8])
transformMethod    <- ifelse(is.na(Args[9]), "rlog",        Args[9])
regionCount        <- ifelse(is.na(Args[10]), 50,           as.integer(Args[10]))

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

# Get counts
countData <- data[,grepl(" count$", names(data)) &
                  !grepl(" normalised count$", names(data)) &
                  !grepl(" end read count$", names(data))]
names(countData) <- gsub(" count$", "", names(countData))

# Transform using DESeq2
dds <- DESeqDataSetFromMatrix(countData, samples, design = ~ condition)
dds <- estimateSizeFactors(dds)
if (transformMethod == "rlog") {
    dds <- rlogTransformation(dds, blind=TRUE)
} else {
    dds <- varianceStabilizingTransformation(dds, blind=TRUE)
}

if (nrow(assay(dds)) < regionCount) {
    regionCount <- nrow(assay(dds))
}
mat <- matrix(nrow=regionCount, ncol=length(levels(samples$condition)),
              dimnames = list(c(), levels(samples$condition)))
for (condition in levels(samples$condition)) {
    if (sum(samples$condition == condition) == 1) {
        mat[,condition] <-
            assay(dds)[1:regionCount,samples$condition == condition]
    } else {
        mat[,condition] <-
            rowMeans(assay(dds)[1:regionCount,samples$condition == condition])
    }
}
rownames(mat) <- data$`Gene name`[1:regionCount]
mat <- mat - rowMeans(mat)
pheatmap(mat, cellheight=10, filename=pngFile)