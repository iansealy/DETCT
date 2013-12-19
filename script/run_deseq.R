library(DESeq2)
library(RColorBrewer)
library(gplots)

Args            <- commandArgs();
countFile       <- Args[4]
designFile      <- Args[5]
outputFile      <- Args[6]
sizeFactorsFile <- Args[7]
qcPdfFile       <- Args[8]
outputFile_Wald      <- Args[9]

countTable    <- read.table(  countFile , header=TRUE, row.names=1 )
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

cdsOneFactFull <- DESeqDataSetFromMatrix(countTable, design, formula(~ condition))
if (numFactors == 2) {
    cdsTwoFactFull <- DESeqDataSetFromMatrix(countTable, design, formula(~ group + condition))
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
    estimateDispersions( cdsOneFactFull)
}, error = function(e) {
    estimateDispersions( cdsOneFactFull, fitType="local" )
})
if (numFactors == 1) {
    plotDispEsts( cdsOneFactFiltPooled )
} else if (numFactors == 2) {
    cdsTwoFactFiltPooledCR <- tryCatch({
        estimateDispersions( cdsTwoFactFilt )
    }, error = function(e) {
        estimateDispersions( cdsTwoFactFilt,fitType="local" )
    })
    cdsTwoFactFullBlind <- tryCatch({
        estimateDispersions( cdsTwoFactFull )
    }, error = function(e) {
        estimateDispersions( cdsTwoFactFull, fitType="local" )
    })
    plotDispEsts( cdsTwoFactFiltPooledCR )
}


res <- nbinomLRT(cdsOneFactFiltPooled, full = ~ condition ,reduced = ~ 1)
results = results(res);
regions = rownames(results);
pval = results$pvalue
padj = results$padj	

dds <- nbinomWaldTest(cdsOneFactFiltPooled)
res_wald = results(dds)


if (numFactors == 2) {
    
    res<-nbinomLRT(cdsTwoFactFiltPooledCR,full= ~ group + condition,reduced= ~ group )
    results = results(res);
   	regions = rownames(results);
	pval = results$pvalue
	padj = results$padj
	
	dds <- nbinomWaldTest(cdsTwoFactFiltPooledCR)
	res_wald = results(dds)
}

plotMA(res)
hist(pval, breaks=100, col="skyblue", border="slateblue",
    main="Histogram of p values")

# Write output
res = data.frame(id=regions, pval=pval, padj=padj)
write.table( res, file=outputFile, col.names=FALSE, row.names=FALSE,
    quote=FALSE, sep="\t" )

res_wald = data.frame(id= rownames(res_wald),
					intercept=res_wald$baseMean,
					log2FoldChange=res_wald$log2FoldChange,
					standard_error=res_wald$lfcSE,
					pval=res_wald$pvalue,
					padj=res_wald$padj)
					
write.table( res_wald, file=outputFile_Wald, col.names=FALSE, row.names=FALSE, quote=FALSE, sep="\t" )		
