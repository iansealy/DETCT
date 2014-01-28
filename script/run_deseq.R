library(DESeq2)
library(RColorBrewer)
library(gplots)

estimate_power<-function(design=design, 
						 res_wald =res_wald, 
						 max_sample_size=100){
	
	# the numerator of the z statistic
	diff_vector=as.vector(res_wald$log2FoldChange)
	
	# to get the st. dev. I need to know the number of samples in the non reference condition level (sibling) 
	# in R the reference level is the first level sorted alphabetically, so get the next one (index 2)
	non_ref_sample_size = length( which ( design$condition == sort(levels(design$condition))[2]));
	non_ref_level_name =  sort(levels(design$condition))[2];
	sd = res_wald$lfcSE * sqrt(non_ref_sample_size )
	
	
	potential_regions<-c();
	for(sample_size in 1:max_sample_size){
		z = diff_vector / (sd / sqrt(sample_size ))
		pvalue = (1-pnorm(z))+pnorm(-z)
		padj=p.adjust(pvalue,"BH")
		potential_regions<-c(potential_regions,length(which(padj <= 0.05)))
	}
	
	increase_vector<-c(0);
	for(i in 2:length(potential_regions)){
		increase = potential_regions[i] - potential_regions[i-1];
		increase_vector<-c(increase_vector,increase);
	}
	
	asterisk_vector = rep('',length(potential_regions))
	max_increase_index = which.max(increase_vector)	
	asterisk_vector[max_increase_index]='*'
	
	power_data = data.frame(c(1:max_sample_size),potential_regions,increase_vector,asterisk_vector)
	colnames(power_data)=c(paste("sample_size",non_ref_level_name,sep="__"),"potential_regions","increase","best")
	
	
	write.table( power_data, file=outputFilePowerMatrix, col.names=colnames(power_data), row.names=FALSE, quote=FALSE, sep="\t" )
	pdf(file = outputPdfPowerCurve);
	plot(power_data$sample_size,power_data$potential_regions,
		type="l", main="z stat power curve",xlab=colnames(power_data)[1],ylab=colnames(power_data)[2])
	abline(v=max_increase_index)

	dev.off();
	
	return(power_data);
}

     

Args            <- commandArgs();
countFile       <- Args[4]
designFile      <- Args[5]
outputFile      <- Args[6]
sizeFactorsFile <- Args[7]
qcPdfFile       <- Args[8]
outputFilePowerMatrix      <- Args[9]
outputPdfPowerCurve      <- Args[10]

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
rs  <- rowSums ( DESeq2::counts ( cdsOneFactFull ))
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
    estimateDispersions( cdsOneFactFilt, fitType="parametric" )
})
cdsOneFactFullBlind <- tryCatch({
    estimateDispersions( cdsOneFactFull)
}, error = function(e) {
    estimateDispersions( cdsOneFactFull, fitType="parametric" )
})
if (numFactors == 1) {
    plotDispEsts( cdsOneFactFiltPooled )
} else if (numFactors == 2) {
    cdsTwoFactFiltPooledCR <- tryCatch({
        estimateDispersions( cdsTwoFactFilt )
    }, error = function(e) {
        estimateDispersions( cdsTwoFactFilt,fitType="parametric" )
    })
    cdsTwoFactFullBlind <- tryCatch({
        estimateDispersions( cdsTwoFactFull )
    }, error = function(e) {
        estimateDispersions( cdsTwoFactFull, fitType="parametric" )
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


# Variance stabilising transformation
vsdOneFactFull <- varianceStabilizingTransformation( cdsOneFactFullBlind )
if (numFactors == 2) {
    vsdTwoFactFull <- varianceStabilizingTransformation( cdsTwoFactFullBlind )
}

## Plot heatmap of counts
select <- order(rowMeans(counts(cdsOneFactFull)), decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(assay(vsdOneFactFull)[select,], col=hmcol, trace="none",margin=c(10, 6))

## Plot heatmap of sample to sample distances
dists <- dist(t(assay(vsdTwoFactFull)))
mat <- as.matrix( dists )
heatmap.2(mat, trace="none", col = rev(hmcol), margin=c(13, 13))

# Plot PCA of samples
print(plotPCA(vsdOneFactFull, intgroup=c("condition")))
if (numFactors == 2) {
    print(plotPCA(vsdTwoFactFull, intgroup=c("group")))
    print(plotPCA(vsdTwoFactFull, intgroup=c("condition", "group")))
}
dev.off()


power_data=estimate_power(design,res_wald,20);


