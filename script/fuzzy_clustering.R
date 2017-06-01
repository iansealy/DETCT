#!/usr/bin/env Rscript

suppressWarnings(library(tcltk))
suppressPackageStartupMessages(library(Mfuzz))
library(naturalsort)

Args           <- commandArgs()
dataFile       <- ifelse(is.na(Args[6]),  "all.tsv",     Args[6])
samplesFile    <- ifelse(is.na(Args[7]),  "samples.txt", Args[7])
outputBase     <- ifelse(is.na(Args[8]),  "all",         Args[8])
numClusters    <- ifelse(is.na(Args[9]),  100,           as.integer(Args[9]))
alphaThreshold <- ifelse(is.na(Args[10]), 0.7,           as.numeric(Args[10]))

# Temporarily include modified version of mfuzz.plot2
mfuzz.plot2.tmp <- function(eset, cl, mfrow = c(1, 1), colo, min.mem = 0, time.labels,
    time.points, ylim.set = c(0, 0), xlab = "Time", ylab = "Expression changes",
    x11 = TRUE, ax.col = "black", bg = "white", col.axis = "black", col.lab = "black",
    col.main = "black", col.sub = "black", col = "black", centre = FALSE, centre.col = "black",
    centre.lwd = 2, Xwidth = 5, Xheight = 5, single = FALSE, ...) {
    clusterindex <- cl[[3]]
    memship <- cl[[4]]
    memship[memship < min.mem] <- -1
    colorindex <- integer(dim(exprs(eset))[[1]])

    if (missing(colo)) {
        colo <- c("#FF0000", "#FF1800", "#FF3000", "#FF4800", "#FF6000", "#FF7800",
            "#FF8F00", "#FFA700", "#FFBF00", "#FFD700", "#FFEF00", "#F7FF00", "#DFFF00",
            "#C7FF00", "#AFFF00", "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00",
            "#20FF00", "#08FF00", "#00FF10", "#00FF28", "#00FF40", "#00FF58", "#00FF70",
            "#00FF87", "#00FF9F", "#00FFB7", "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF",
            "#00CFFF", "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF", "#0040FF",
            "#0028FF", "#0010FF", "#0800FF", "#2000FF", "#3800FF", "#5000FF", "#6800FF",
            "#8000FF", "#9700FF", "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF",
            "#FF00D7", "#FF00BF", "#FF00A7", "#FF008F", "#FF0078", "#FF0060", "#FF0048",
            "#FF0030", "#FF0018")
    } else {
        if (colo == "fancy") {
            fancy.blue <- c(c(255:0), rep(0, length(c(255:0))), rep(0, length(c(255:150))))
            fancy.green <- c(c(0:255), c(255:0), rep(0, length(c(255:150))))
            fancy.red <- c(c(0:255), rep(255, length(c(255:0))), c(255:150))
            colo <- rgb(b = fancy.blue/255, g = fancy.green/255, r = fancy.red/255)
        }
    }
    colorseq <- seq(0, 1, length = length(colo))

    for (j in 1:dim(cl[[1]])[[1]]) {
        if (single)
            j <- single
        tmp <- exprs(eset)[clusterindex == j, , drop = FALSE]
        tmpmem <- memship[clusterindex == j, j]
        if (((j - 1)%%(mfrow[1] * mfrow[2])) == 0 | single) {
            if (x11)
                X11(width = Xwidth, height = Xheight)
            if (sum(clusterindex == j) == 0) {
                ymin <- -1
                ymax <- +1
            } else {
                ymin <- min(tmp)
                ymax <- max(tmp)
            }
            if (sum(ylim.set == c(0, 0)) == 2) {
                ylim <- c(ymin, ymax)
            } else {
                ylim <- ylim.set
            }
            if (!is.na(sum(mfrow))) {
                par(mfrow = mfrow, bg = bg, col.axis = col.axis, col.lab = col.lab,
                  col.main = col.main, col.sub = col.sub, col = col)
            } else {
                par(bg = bg, col.axis = col.axis, col.lab = col.lab, col.main = col.main,
                  col.sub = col.sub, col = col)
            }
            xlim.tmp <- c(1, dim(exprs(eset))[[2]])
            if (!(missing(time.points)))
                xlim.tmp <- c(min(time.points), max(time.points))
            plot.default(x = NA, xlim = xlim.tmp, ylim = ylim, xlab = xlab, ylab = ylab,
                main = paste("Cluster", j), axes = FALSE, ...)
            if (missing(time.labels) && missing(time.points)) {
                axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]), col = ax.col,
                  ...)
                axis(2, col = ax.col, ...)
            }
            if (missing(time.labels) && !(missing(time.points))) {
                axis(1, time.points, 1:length(time.points), time.points, col = ax.col,
                  ...)
                axis(2, col = ax.col, ...)
            }
            if (missing(time.points) & !(missing(time.labels))) {
                axis(1, 1:dim(exprs(eset))[[2]], time.labels, col = ax.col, ...)
                axis(2, col = ax.col, ...)
            }
            if (!(missing(time.points)) & !(missing(time.labels))) {
                axis(1, time.points, time.labels, col = ax.col, ...)
                axis(2, col = ax.col, ...)
            }
        } else {
            if (sum(clusterindex == j) == 0) {
                ymin <- -1
                ymax <- +1
            } else {
                ymin <- min(tmp)
                ymax <- max(tmp)
            }
            if (sum(ylim.set == c(0, 0)) == 2) {
                ylim <- c(ymin, ymax)
            } else {
                ylim <- ylim.set
            }
            xlim.tmp <- c(1, dim(exprs(eset))[[2]])
            if (!(missing(time.points)))
                xlim.tmp <- c(min(time.points), max(time.points))
            plot.default(x = NA, xlim = xlim.tmp, ylim = ylim, xlab = xlab, ylab = ylab,
                main = paste("Cluster", j), axes = FALSE, ...)
            if (missing(time.labels) && missing(time.points)) {
                axis(1, 1:dim(exprs(eset))[[2]], c(1:dim(exprs(eset))[[2]]), col = ax.col,
                  ...)
                axis(2, col = ax.col, ...)
            }
            if (missing(time.labels) && !(missing(time.points))) {
                axis(1, time.points, 1:length(time.points), time.points, col = ax.col,
                  ...)
                axis(2, col = ax.col, ...)
            }
            if (missing(time.points) & !(missing(time.labels))) {
                axis(1, 1:dim(exprs(eset))[[2]], time.labels, col = ax.col, ...)
                axis(2, col = ax.col, ...)
            }
            if (!(missing(time.points)) & !(missing(time.labels))) {
                axis(1, time.points, time.labels, col = ax.col, ...)
                axis(2, col = ax.col, ...)
            }
        }

        if (length(tmpmem) > 0) {
            for (jj in 1:(length(colorseq) - 1)) {
                tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <= colorseq[jj + 1])
                if (sum(tmpcol) > 0) {
                  tmpind <- which(tmpcol)
                  for (k in 1:length(tmpind)) {
                    if (missing(time.points)) {
                      lines(tmp[tmpind[k], ], col = colo[jj])
                    } else lines(time.points, tmp[tmpind[k], ], col = colo[jj])
                  }
                }
            }
        }
        if (centre) {
            lines(cl[[1]][j, ], col = centre.col, lwd = centre.lwd)
        }
        if (single)
            return()
    }
}

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

# Subset and reorder count data
countData <- countData[, row.names(samples)]

# Get median counts
medianData <- matrix(nrow=nrow(countData),
                     ncol=length(levels(samples$condition)),
                     dimnames=list(rownames(countData),
                                   naturalsort(levels(samples$condition))))
for (condition in naturalsort(levels(samples$condition))) {
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
pdf(paste0(outputBase, '-', numClusters, '-mfuzz.pdf'))
mfuzz.plot2.tmp(eset.s, cl=cl, mfrow=c(4,4), x11=FALSE, centre=TRUE,
           time.labels=colnames(medianData))
graphics.off()
pdf(paste0(outputBase, '-', numClusters, '-', alphaThreshold, '-mfuzz.pdf'))
mfuzz.plot2.tmp(eset.s, cl=cl, mfrow=c(4,4), x11=FALSE, centre=TRUE,
           time.labels=colnames(medianData), min.mem=alphaThreshold)
graphics.off()

# Output fuzzy clusters
colours <- as.numeric(samples$condition)
acore.list <- acore(eset.s, cl=cl, min.acore=alphaThreshold)
for (i in 1:numClusters) {
    data.subset <- data[c(as.character(acore.list[[i]][,1])),]
    data.subset$membership <- acore.list[[i]][,2]
    if (grepl("csv$", dataFile)) {
        write.csv(data.subset, file=paste0(outputBase, '-', numClusters, '-',
                                           alphaThreshold, '-', i, '.csv'),
                  row.names=FALSE)
    } else {
        write.table(data.subset, file=paste0(outputBase, '-', numClusters, '-',
                                             alphaThreshold, '-', i, '.tsv'),
                    quote=FALSE, sep='\t', row.names=FALSE, col.names=TRUE)
    }
    data.subset <- data.subset[,grepl(" normalised count$", names(data.subset))]
    names(data.subset) <- gsub(" normalised count$", "", names(data.subset))
    data.subset <- data.subset[, row.names(samples)]
    if (nrow(data.subset) == 0) {
        next
    }
    # Plot counts
    pdf(paste0(outputBase, '-', numClusters, '-', alphaThreshold, '-', i,
               '-counts.pdf'))
    for (j in 1:nrow(data.subset)) {
        counts <- data.subset[j,]
        par(mar=c(8.1, 4.1, 4.1, 2.1), xpd=TRUE)
        plot(as.numeric(counts), axes=FALSE, ann=FALSE, pch=21, bg=colours)
        axis(1, at=1:ncol(data.subset), lab=names(data.subset), las=2,
             cex.axis=0.5)
        axis(2)
        title(main=sprintf("%s:%d-%d %s\n%.2f",
                           data.subset[j,"Chr"],
                           data.subset[j,"Region.start"],
                           data.subset[j,"Region.end"],
                           data.subset[j,"Gene.name"],
                           data.subset[j,"Adjusted.p.value"]))
        title(xlab="")
        title(ylab="Normalised Counts")
        legend("topright", inset=c(0, -0.1), levels(samples$condition), pch=21,
               pt.bg=1:length(levels(samples$condition)))
    }
    graphics.off()
}
