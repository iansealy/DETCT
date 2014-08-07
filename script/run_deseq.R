library(DESeq2)
library(RColorBrewer)
library(gplots)

Args             <- commandArgs()
countFile        <- Args[4]
samplesFile      <- Args[5]
outputFile       <- Args[6]
sizeFactorsFile  <- Args[7]
qcPdfFile        <- Args[8]
filterPercentile <- as.numeric( Args[9] )


# Get data and samples
countData     <- read.table(   countFile, header=TRUE, row.names=1 )
samples       <- read.table( samplesFile, header=TRUE, row.names=1 )
numFactors    <- ncol( samples )
numConditions <- nlevels( samples$condition )

# Can only handle one or two factors
if (numFactors > 2) {
    stop("Too many factors")
}

# If not two conditions then just write empty output
if (numConditions != 2) {
    dds <- DESeqDataSetFromMatrix(countData, samples, design = ~ 1)
    write.table( c(), file=outputFile, col.names=FALSE, row.names=FALSE,
        quote=FALSE, sep="\t" )
    dds <- estimateSizeFactors(dds)
    write.table( sizeFactors( dds ), file=sizeFactorsFile, col.names=FALSE,
        row.names=FALSE, quote=FALSE, sep="\t" )
    stop("Not two conditions")
}

# Optionally remove regions with sum of counts below the specified percentile
if (filterPercentile) {
    rs  <- rowSums( countData )
    use <- (rs > quantile(rs, probs=filterPercentile/100))
    countData <- countData[ use, ]
}

# Create DESeqDataSet (with design according to number of factors)
dds <- DESeqDataSetFromMatrix(countData, samples, design = ~ condition)
if (numFactors == 2) {
    design(dds) <- formula(~ group * condition)
}

# Ensure control level (usually "sibling") is first level (i.e. before "mutant")
colData(dds)$condition <- factor(colData(dds)$condition,
    levels=rev(levels(colData(dds)$condition)))

# Differential expression analysis
dds <- DESeq(dds)

# Find the condition term (does not contain '.' and constains 'condition')
results_names=resultsNames(dds)
non_interaction_terms = which(grepl("\\.",results_names,perl =T) == FALSE)
all_condition_terms_index = which(grepl("condition",results_names) == TRUE)
condition_terms_index = intersect(non_interaction_terms,all_condition_terms_index)

# create matrix with pvalues of all terms in the model
pvalues_matrix = matrix(ncol=0,nrow=nrow(dds))
col_names<-c();
for (i  in 1:length(results_names) )
{
	equation_term = results_names[i]; 
	results = results(dds, name= equation_term)
	palues = as.matrix(results$pvalue)
	padj =   as.matrix(results$padj)
	pvalues_matrix = cbind(pvalues_matrix,palues)
	pvalues_matrix = cbind(pvalues_matrix,padj)
	col_names<-c(col_names,paste(equation_term,"_pvalue",sep=""),paste(equation_term,"_pvalue_adjusted",sep=""))
}
colnames(pvalues_matrix)=col_names;



# Write output
out <- data.frame(pvalues_matrix, row.names=rownames(res))
write.table( out, file=outputFile, col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t" )
write.table( sizeFactors( dds ), file=sizeFactorsFile, col.names=FALSE,row.names=FALSE, quote=FALSE, sep="\t" )

# Data transformations for QC
rld <- rlogTransformation(dds, blind=TRUE)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

# Write QC graphs to PDF
pdf(qcPdfFile)

# MA-plot
plotMA(dds)

# Heatmaps of counts
select <- order(rowMeans(counts(dds, normalized=TRUE)), decreasing=TRUE)[1:30]
hmcol  <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(assay(rld)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE,
    scale="none", dendrogram="none", trace="none", margin=c(10, 6))
heatmap.2(assay(vsd)[select,], col = hmcol, Rowv = FALSE, Colv = FALSE,
    scale="none", dendrogram="none", trace="none", margin=c(10, 6))

# Heatmaps of sample to sample distances
heatmap.2(as.matrix(dist(t(assay(rld)))), trace="none", col = rev(hmcol),
    margin=c(13, 13))
heatmap.2(as.matrix(dist(t(assay(vsd)))), trace="none", col = rev(hmcol),
    margin=c(13, 13))

# PCA of samples
print(plotPCA(rld, intgroup=c("condition")))
if (numFactors == 2) {
    print(plotPCA(rld, intgroup=c("group")))
    print(plotPCA(rld, intgroup=c("condition", "group")))
}
print(plotPCA(vsd, intgroup=c("condition")))
if (numFactors == 2) {
    print(plotPCA(vsd, intgroup=c("group")))
    print(plotPCA(vsd, intgroup=c("condition", "group")))
}

# Dispersion plot
plotDispEsts(dds)

# Plot optimisation of independent filtering using mean of normalised counts
# See "3.6 Independent filtering of results" of DESeq2 vignette
plot(attr(res, "filterNumRej"), type="b", xlab="quantile",
    ylab="number of rejections")

# Plot p-value over mean of counts
plot(res$baseMean+1, -log10(res$pvalue), log="x",
    xlab="mean of normalised counts", ylab=expression(-log[10](pvalue)),
    ylim=c(0,30), cex=.4, col=rgb(0,0,0,.3))

# p-value histogram
use <- res$baseMean > attr(res,"filterThreshold")
h1 <- hist(res$pvalue[!use], breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use],  breaks=0:50/50, plot=FALSE)
colori <- c("filtered"="khaki", "unfiltered"="powderblue")
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE, col = colori,
    space = 0, main = "", xlab="p-value", ylab="frequency")
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0,1)),
    adj = c(0.5,1.7), xpd=NA)
legend("topright", fill=rev(colori), legend=rev(names(colori)))

# Calculate power

# Constants
max_pairs <- 20
sig_level <- 0.05
num_pairs <- nrow(samples) / 2 # Assume even

# Estimate number of significant regions for different numbers of pairs
sd = res$lfcSE * sqrt(num_pairs)
potential_sig <- c()
for (pairs in 1:max_pairs) {
    z <- as.vector(res$log2FoldChange) / ( sd / sqrt(pairs) )
    pvalue <- ( 1 - pnorm(z) ) + pnorm(-z)
    padj   <- p.adjust(pvalue, "BH")
    potential_sig <- c(potential_sig, length(which(padj <= sig_level)))
}

# Calculate increase in potentially significant regions
increase <- c(0)
for (i in 2:length(potential_sig)) {
    increase <- c(increase, potential_sig[i] - potential_sig[i-1])
}

# Identify optimal number of pairs
optimal <- rep('', length(potential_sig))
optimal[which.max(increase)] <- '*'

# Output power data
power <- data.frame(c(1:max_pairs), potential_sig, increase, optimal)
colnames(power) <- c("Pairs", "Significant", "Increase", "Optimal")
textplot(power, show.rownames=FALSE)

# Plot power curve
plot(power$Pairs,power$Significant, type="l", main="Z-statistic power curve",
    xlab="Number of pairs", ylab="Number of significant regions")
abline(v=which.max(increase))

dev.off()
