library(DESeq)
library(RColorBrewer)
library(gplots)

Args            <- commandArgs();
countFile       <- Args[4]
designFile      <- Args[5]
outputFile      <- Args[6]
sizeFactorsFile <- Args[7]
qcPdfFile       <- Args[8]

# Get data and design
countTable    <- read.table(  countFile, header=TRUE, row.names=1 )
design        <- read.table( designFile, header=TRUE, row.names=1 )
numFactors    <- ncol(design)
numConditions <- nlevels(design$condition)

# Can only handle one or two factors
if (numFactors > 2) {
    stop("Too many factors")
}

# Write QC graphs to PDF
pdf(qcPdfFile)

# Create CountDataSets
cdsOneFactFull <- newCountDataSet( countTable, design$condition )
if (numFactors == 2) {
    cdsTwoFactFull <- newCountDataSet( countTable, design )
}

# Remove regions with sum of counts below the 40th quantile
# See "5 Independent filtering and multiple testing" of
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq/inst/doc/DESeq.pdf
rs  <- rowSums ( counts ( cdsOneFactFull ))
use <- (rs > quantile(rs, probs=0.4))
cdsOneFactFilt <- cdsOneFactFull[ use, ]
if (numFactors == 2) {
    cdsTwoFactFilt <- cdsTwoFactFull[ use, ]
}

# Normalise
cdsOneFactFull <- estimateSizeFactors( cdsOneFactFull )
cdsOneFactFilt <- estimateSizeFactors( cdsOneFactFilt )
if (numFactors == 2) {
    cdsTwoFactFull <- estimateSizeFactors( cdsTwoFactFull )
    cdsTwoFactFilt <- estimateSizeFactors( cdsTwoFactFilt )
}
write.table( sizeFactors( cdsOneFactFull ), file=sizeFactorsFile,
    col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t" )

# If not two conditions then write empty output
if (numConditions != 2) {
    write.table( c(), file=outputFile, col.names=FALSE, row.names=FALSE,
        quote=FALSE, sep="\t" )
    stop("Not two conditions")
}

# Estimate variance
cdsOneFactFiltPooled <- tryCatch({
    estimateDispersions( cdsOneFactFilt )
}, error = function(e) {
    estimateDispersions( cdsOneFactFilt, fitType="local" )
})
cdsOneFactFullBlind <- tryCatch({
    estimateDispersions( cdsOneFactFull, method="blind" )
}, error = function(e) {
    estimateDispersions( cdsOneFactFull, method="blind", fitType="local" )
})
if (numFactors == 1) {
    plotDispEsts( cdsOneFactFiltPooled )
} else if (numFactors == 2) {
    cdsTwoFactFiltPooledCR <- tryCatch({
        estimateDispersions( cdsTwoFactFilt, method="pooled-CR" )
    }, error = function(e) {
        estimateDispersions( cdsTwoFactFilt, method="pooled-CR", fitType="local" )
    })
    cdsTwoFactFullBlind <- tryCatch({
        estimateDispersions( cdsTwoFactFull, method="blind" )
    }, error = function(e) {
        estimateDispersions( cdsTwoFactFull, method="blind", fitType="local" )
    })
    plotDispEsts( cdsTwoFactFiltPooledCR )
}

# Compare conditions
conditions <- levels(design$condition)
res <- nbinomTest( cdsOneFactFiltPooled, conditions[1], conditions[2] )
if (numFactors == 2) {
    fit1 <- fitNbinomGLMs( cdsTwoFactFiltPooledCR, count ~ group + condition )
    fit0 <- fitNbinomGLMs( cdsTwoFactFiltPooledCR, count ~ group )
    res$pval <- nbinomGLMTest( fit1, fit0 )
    res$padj <- p.adjust( res$pval, method="BH" )
}
plotMA(res)
hist(res$pval, breaks=100, col="skyblue", border="slateblue",
    main="Histogram of p values")

# Write output
res = data.frame(id=res$id, pval=res$pval, padj=res$padj)
write.table( res, file=outputFile, col.names=FALSE, row.names=FALSE,
    quote=FALSE, sep="\t" )

# Variance stabilising transformation
vsdOneFactFull <- varianceStabilizingTransformation( cdsOneFactFullBlind )
if (numFactors == 2) {
    vsdTwoFactFull <- varianceStabilizingTransformation( cdsTwoFactFullBlind )
}

# Plot heatmap of counts
select <- order(rowMeans(counts(cdsOneFactFull)), decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(exprs(vsdOneFactFull)[select,], col=hmcol, trace="none",
    margin=c(10, 6))

# Plot heatmap of sample to sample distances
dists <- dist( t( exprs(vsdOneFactFull) ) )
mat <- as.matrix( dists )
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))

# Plot PCA of samples
print(plotPCA(vsdOneFactFull, intgroup=c("condition")))
if (numFactors == 2) {
    print(plotPCA(vsdTwoFactFull, intgroup=c("group")))
    print(plotPCA(vsdTwoFactFull, intgroup=c("condition", "group")))
}

dev.off()
