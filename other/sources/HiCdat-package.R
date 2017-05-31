# This file is part of HiCdat.
# 
# HiCdat is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# HiCdat is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# See <http://www.gnu.org/licenses/> for a a copy of the GNU General Public License.

cat("##############################################################################################
###                                                                                        ###
### IMPORTAN NOTE:                                                                         ###
### this package requires organism specific functions, load them with:                     ###
### f.source.organism.specific.code('HiCdat-A-thaliana-TAIR10.R')                          ###
###                                                                                        ###
### these files can be automatically created using HiCdatPre                               ###
### download the binary from https://github.com/MWSchmid/HiCdat                            ###
###                                                                                        ###
##############################################################################################
")

#'@title a global variable for plotting control.
#'@note set to true if you wish all plots to be saved as svg and converted to png (requires the command-line tool rsvg-convert).
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
GLOBAL_VARIABLE_USE_SVG_AND_RSVG_CONVERT <- FALSE

#'@title another global variable for plotting control.
#'@note set to true if you wish some plots to be saved as svg instead of tiff/png.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
GLOBAL_VARIABLE_USE_SVG <- TRUE; if (Sys.info()['sysname'] == "Darwin") { GLOBAL_VARIABLE_USE_SVG <- FALSE }

#'@title Add organism-specific R code.
#'@param rScript the file with the organism-specific R code.
#'@note The file required can be built using HiCdatPre (\url{https://github.com/MWSchmid/HiCdat}). 
#'@export
f.source.organism.specific.code <- function(rScript) {
	source(rScript)
	assignInMyNamespace("f.internal.get.chrom.sizes", f.get.chrom.sizes)
	assignInMyNamespace("f.internal.get.frag.list", f.get.frag.list)
	assignInMyNamespace("f.internal.get.relevant.chromosomes", f.get.relevant.chromosomes)
}

#'@title Definition of the chromosome sizes. 
#'@return a list with chromosome sizes.
#'@note this is an empty function warning the user if it was not overwritten.
#'Create the file with organism-specific code using HiCdatPre (\url{https://github.com/MWSchmid/HiCdat})
#'and load it after loading the HiCdatR package with \code{\link{f.source.organism.specific.code}}. 
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.get.chrom.sizes <- function() {
	stop("If you see this, you likely forgot to run\nf.source.organism.specific.code('/path/to/the/file')\nThis file can be created with HiCdatPre available on:\ngithub.com/MWSchmid/HiCdat")
}

#'@title Restriction fragment indices of the chromosomes. 
#'@return a list containing the start and end restriction fragments for each chromosome.
#'@note this is an empty function warning the user if it was not overwritten.
#'Create the file with organism-specific code using HiCdatPre (\url{https://github.com/MWSchmid/HiCdat})
#'and load it after loading the HiCdatR package with \code{\link{f.source.organism.specific.code}}. 
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.get.frag.list <- function() {
	stop("If you see this, you likely forgot to run\nf.source.organism.specific.code('/path/to/the/file')\nThis file can be created with HiCdatPre available on:\ngithub.com/MWSchmid/HiCdat")
}

#'@title Definition of chromosomes of interest.
#'@return a vector containing the chromosome names of interest.
#'@note this is an empty function warning the user if it was not overwritten.
#'Create the file with organism-specific code using HiCdatPre (\url{https://github.com/MWSchmid/HiCdat})
#'and load it after loading the HiCdatR package with \code{\link{f.source.organism.specific.code}}. 
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.get.relevant.chromosomes <- function() {
	stop("If you see this, you likely forgot to run\nf.source.organism.specific.code('/path/to/the/file')\nThis file can be created with HiCdatPre available on:\ngithub.com/MWSchmid/HiCdat")
}

#'@title Get the size of relevant chromosomes.
#'@return a list of chromosome sizes, accessible via chromosome name.
#'@note requires an organism-specific R-script, see \url{https://github.com/MWSchmid/HiCdat}.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.get.relevant.chrom.sizes <- function() {
	chromToSize <- f.internal.get.chrom.sizes()
	chroms <- f.internal.get.relevant.chromosomes()
	out <- list()
	for (chrom in chroms) { out[[chrom]] <- chromToSize[[chrom]] }
	return(out)
}

#'@title Get a list containing the start and end bins of ALL chromosomes.
#'@param binSize size of the genomic bins in bp, if set to 0, restriction fragments are used instead.
#'@return a list containing the start and end bins for each chromosome.
#'@note requires an organism-specific R-script, see \url{https://github.com/MWSchmid/HiCdat}.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.get.se.list}}.
f.get.se.list.with.irrelevant.chromosomes <- function(binSize) {
	if (binSize == 0) {
		out <- f.internal.get.frag.list()
		out[["ALL"]] <- c(1, max(do.call("rbind", out)[,2]))
	} else {
		chromToSize <- f.internal.get.chrom.sizes()
		out <- list()
		prevEnd <- 0
		for (chrom in names(chromToSize)) {
			s <- prevEnd + 1
			e <- s + floor(chromToSize[[chrom]]/binSize)
			out[[chrom]] <- c(s, e)
			prevEnd <- e
		}
		out[["ALL"]] <- c(1, prevEnd)
	}
	return(out)
}

#'@title Get a list containing the start and end bins of the RELEVANT chromosomes.
#'@param binSize size of the genomic bins in bp, if set to 0, restriction fragments are used instead.
#'@return a list containing the start and end bins for each chromosome.
#'@note requires an organism-specific R-script, see \url{https://github.com/MWSchmid/HiCdat}.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@export
f.get.se.list <- function(binSize) {
	if (binSize == 0) {
		out <- f.internal.get.frag.list()
		out[["ALL"]] <- c(1, max(do.call("rbind", out)[,2]))
	} else {
		chromToSize <- f.internal.get.chrom.sizes()
		chroms <- f.internal.get.relevant.chromosomes()
		out <- list()
		prevEnd <- 0
		for (chrom in chroms) {
			s <- prevEnd + 1
			e <- s + floor(chromToSize[[chrom]]/binSize)
			out[[chrom]] <- c(s, e)
			prevEnd <- e
		}
		out[["ALL"]] <- c(1, prevEnd)
	}
	return(out)
}

#'@title Convert a genomic coordinate (chrom, position) into its binned coordinate (index).
#'@param chrom a chromosome identifier.
#'@param position a position within the chromosome.
#'@param binSize size of genomic bins in bp (0 will not produce correct results).
#'@param seList a list with first and last bins per chromosome (\code{\link{f.get.se.list}}).
#'@return the binned coordinate.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.translate.chrom.pos.vector.to.index}}, \code{\link{f.translate.index.to.chrom.pos}}, and
#'\code{\link{f.translate.index.vector.to.chrom.pos}}.
#'@export
f.translate.chrom.pos.to.index <- function(chrom, position, seList, binSize) {
	out <- floor(position/binSize)+seList[[chrom]][1]
	return(out)
}

#'@title Convert multiple genomic coordinates (chrom, position) into their binned coordinates (index).
#'@param chromVec a vector of chromosomes in the genome.
#'@param positionVec a vector of positions within the chromosomes.
#'@param binSize size of genomic bins in bp (0 will not produce correct results).
#'@param seList a list with first and last bins per chromosome (\code{\link{f.get.se.list}}).
#'@return a vector with binned coordinates.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.translate.chrom.pos.to.index}}, \code{\link{f.translate.index.to.chrom.pos}}, and
#'\code{\link{f.translate.index.vector.to.chrom.pos}}.
#'@export
f.translate.chrom.pos.vector.to.index <- function(chromVec, positionVec, seList, binSize) {
	endVec <- do.call("rbind",seList[chromVec])[,1]
	out <- floor(positionVec/binSize)+endVec
	return(out)
}

#'@title Convert a binned coordinate (index) into its genomic coordinate (chrom, start, end).
#'@param index a genomic bin number.
#'@param binSize size of genomic bins in bp (0 will not produce correct results).
#'@param seList a list with first and last bins per chromosome (\code{\link{f.get.se.list}}).
#'@return the genomic coordinate corresponding to the index.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.translate.chrom.pos.to.index}}, \code{\link{f.translate.chrom.pos.vector.to.index}}, and
#'\code{\link{f.translate.index.vector.to.chrom.pos}}.
#'@export
f.translate.index.to.chrom.pos <- function(index, seList, binSize) {
	chromToSize <- f.get.relevant.chrom.sizes()
	chromIntervals <- do.call("rbind", seList)[,2][names(chromToSize)]+1
	chrom <- names(chromIntervals)[findInterval(index, chromIntervals, rightmost.closed = TRUE)+1]
	s <- (index - seList[[chrom]][1])*binSize
	e <- s + binSize
	e[e>chromToSize[chrom]] <- chromToSize[[chrom]]
	out <- data.frame(chrom = chrom, start = s, end = e, stringsAsFactors = FALSE)
	return(out)
}

#'@title Convert a binned coordinates (index) into their genomic coordinates (chrom, start, end).
#'@param indexVec a vector of genomic bin numbers.
#'@param binSize size of genomic bins in bp (0 will not produce correct results).
#'@param seList a list with first and last bins per chromosome (\code{\link{f.get.se.list}}).
#'@return the genomic coordinates corresponding to the indices.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.translate.chrom.pos.to.index}}, \code{\link{f.translate.chrom.pos.vector.to.index}}, and
#'\code{\link{f.translate.index.to.chrom.pos}}.
#'@export
f.translate.index.vector.to.chrom.pos <- function(indexVec, seList, binSize) {
	chromToSize <- f.get.relevant.chrom.sizes()
	chromIntervals <- do.call("rbind", seList)[,2][names(chromToSize)]+1
	chrom <- names(chromIntervals)[findInterval(indexVec, chromIntervals, rightmost.closed = TRUE)+1]
	s <- (indexVec - do.call("rbind", seList)[chrom,1])*binSize
	e <- s + binSize
	csv <- unlist(chromToSize[chrom])
	toReplace <- e>csv
	e[toReplace] <- csv[toReplace]
	out <- data.frame(chrom = chrom, start = s, end = e, stringsAsFactors = FALSE)
	return(out)
}

#'@title Summarize the diagonal and the neighboring fields of a matrix into a new vector.
#'@param dataMatrix a HiC interaction matrix (see \code{\link{f.load.one.sample}}).
#'@param fieldsToAdd how many fields to add on both sides (e.g. 1 means that 3x3 fields are summarized).
#'@param summaryFunction the function used to summarize the fields - e.g. mean or sum.
#'@param useLog if given, the resulting vector will be logarithmized.
#'@return a vector of the summarized diagonal.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@export
f.summary.along.diagonal <- function(dataMatrix, fieldsToAdd = 1, summaryFunction = sum, useLog = FALSE) {
	from <- (1:ncol(dataMatrix))-fieldsToAdd
	to <- (1:ncol(dataMatrix))+fieldsToAdd
	from[from < 1] <- 1
	to[to > ncol(dataMatrix)] <- ncol(dataMatrix)
	winNums <- cbind(from, to, from, to)
	out <- apply(winNums, 1, function(x) summaryFunction(dataMatrix[x[1]:x[2],x[3]:x[4]], na.rm = TRUE))
	if (useLog) {out <- log2(out+1)}
	return(out)
}

#'@title Extract a subset of the HiC interaction matrix using genomic coordinates (chrom, start, end).
#'@param dataMatrix a HiC interaction matrix (see \code{\link{f.load.one.sample}}).
#'@param binSize size of genomic bins in bp (0 will not produce correct results).
#'@param xChrom first chromosome.
#'@param xStart start on first chromosome.
#'@param xEnd end on first chromosome.
#'@param yChrom second chromosome.
#'@param yStart start on second chromosome.
#'@param yEnd end on second chromosome.
#'@return a HiC interaction matrix (a subset of dataMatrix).
#'@note to extract entire chromosomes, set start to 0 and end to 1e100.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@export
f.extract.subset <- function(dataMatrix, binSize, xChrom, yChrom, xStart = 1, yStart = 1, xEnd = 1e100, yEnd = 1e100) {
	seList <- f.get.se.list(binSize)
	chromSizes <- f.internal.get.chrom.sizes()
	if (xEnd > chromSizes[[xChrom]]) {xEnd <- chromSizes[[xChrom]]}
	xS <- f.translate.chrom.pos.to.index(xChrom, xStart, seList, binSize)
	xE <- f.translate.chrom.pos.to.index(xChrom, xEnd, seList, binSize)
	if (yEnd > chromSizes[[yChrom]]) {yEnd <- chromSizes[[yChrom]]}
	yS <- f.translate.chrom.pos.to.index(yChrom, yStart, seList, binSize)
	yE <- f.translate.chrom.pos.to.index(yChrom, yEnd, seList, binSize)
	return(dataMatrix[xS:xE, yS:yE])
}

#'@title Normalize HiC data iteratively.
#'@param dataMatrix a HiC interaction matrix (see \code{\link{f.load.one.sample}}).
#'@param repetitions number of iterations.
#'@return a normalized HiC interaction matrix.
#'@references Zhang, Y. and McCord, R. P. and Ho, Y. J. and Lajoie, B. R. and Hildebrand, D. G. and
#'Simon A. C. and Becker, M. S. and Alt, F. W. and Dekker, J. (2012) Spatial Organization of the 
#'Mouse Genome and Its Role in Recurrent Chromosomal Translocations. \emph{Cell} \bold{148}, 908--921.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@export
f.normalize.like.zhang <- function(dataMatrix, repetitions) {
	cat(paste("normalizing with", repetitions, "iterations\n"))
	for (i in 1:repetitions) {
		cat('#')
		rs <- apply(dataMatrix, 1, sum) / sum(dataMatrix)
		winMatCor <- (rs %*% t(rs)) * sum(dataMatrix)
		dataMatrix <- dataMatrix / winMatCor
		dataMatrix[is.na(dataMatrix)] <- 0
	}
	cat("\n")
	return(dataMatrix)
}

#'@title Normalize HiC data using linear regression (HiCNorm).
#'@param dataMatrix a HiC interaction matrix (see \code{\link{f.load.one.sample}}).
#'@param binSize size of genomic bins in bp, if set to 0, restriction fragments are used instead.
#'@param annotation a dataframe holding the annotation of the bins (see \code{\link{f.read.annotation}}).
#'@param lenCol name of the column holding the length of the bins.
#'@param gccCol name of the column holding the GC-content of the bins.
#'@param mapCol name of the column holding the alignability of the bins.
#'@param useNegativeBinomial use negative binomial model (default is Poisson).
#'@return a normalized HiC interaction matrix.
#'@examples \dontrun{
#'normalizedDataMatrix <- f.normalize.like.hu(
#'  dataMatrix = dataMatrixSampleX,
#'  binSize = 1e5,
#'  annotation = annotationTable,
#'  lenCol = "length",
#'  gccCol = "gcContent",
#'  mapCol = "mappability",
#'  useNegativeBinomial = FALSE
#')}
#'@details
#'Note that the three parameters fragment length, GC-content, and mappability are not defined per default in the annotation tables.
#'However, you can use HiCdatPre \url{https://github.com/MWSchmid/HiCdat} to obtain them. For example:
#'\itemize{
#'  \item{}{fragment length can be calculated directly in R: annotation$length <- annotation$end - annotation$start.}
#'  \item{}{GC-content can be imported as a density-feature. Instead of using a regular DNA-methylation table, 
#'  one can supply a table where the CG-positions are marked as methylated and the non-CG positions are marked as unmethylated.
#'  An example for an artificial chromosome ``Chr1'' starting with the sequence ACGTA:
#' \tabular{ccc}{
#'  Chr1 \tab 0 \tab u\cr
#'  Chr1 \tab 0 \tab u\cr
#'  Chr1 \tab 1 \tab m\cr
#'  Chr1 \tab 2 \tab m\cr
#'  Chr1 \tab 3 \tab u\cr
#'  Chr1 \tab 4 \tab u
#'  }
#' }
#'  \item{}{for mappability, one can align either artificial reads (from a chopped genome) or real genomic sequencing reads and 
#'  import them as a density-feature using.}
#' }
#'@references Hu, M. and Deng, K. and Selvaraj, S. and Qin, Z. S. and Ren, B. and Liu, J. S. (2012)
#'HiCNorm: removing biases in Hi-C data via Poisson regression. \emph{Bioinformatics} \bold{28}, 3131--3133.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@export
f.normalize.like.hu <- function(dataMatrix, binSize, annotation, lenCol, gccCol, mapCol, useNegativeBinomial = FALSE) {
	# set minimal GC and mappability to 0.0001
	annotation[annotation[,mapCol] < 0.0001,mapCol] <- 0.0001
	annotation[annotation[,gccCol] < 0.0001,gccCol] <- 0.0001
	columnsToUse <- c(lenCol, gccCol, mapCol)
	seList <- f.get.se.list(binSize)
	chroms <- f.internal.get.relevant.chromosomes()
	out <- matrix(0, nrow = nrow(dataMatrix), ncol = ncol(dataMatrix))
	for (i in 1:length(chroms)) {
		chromA <- chroms[i]
		sA <- seList[[chromA]][1]
		eA <- seList[[chromA]][2]
		for (j in i:length(chroms)) {
			chromB <- chroms[j]
			sB <- seList[[chromB]][1]
			eB <- seList[[chromB]][2]
			if (chromA == chromB) { 
				subAnno <- annotation[sA:eA, columnsToUse]
				out[sA:eA, sB:eB] <- f.internal.normalize.intra.like.hu(dataMatrix[sA:eA, sB:eB], subAnno, useNegativeBinomial)
			} else {
				subAnnoA <- annotation[sA:eA, columnsToUse]
				subAnnoB <- annotation[sB:eB, columnsToUse]
				out[sA:eA, sB:eB] <- f.internal.normalize.inter.like.hu(dataMatrix[sA:eA, sB:eB], subAnnoA, subAnnoB, useNegativeBinomial)
				out[sB:eB, sA:eA] <- f.internal.normalize.inter.like.hu(dataMatrix[sB:eB, sA:eA], subAnnoB, subAnnoA, useNegativeBinomial)
			}
		}
	}
	return(out)
}

#'@title Normalize HiC data using linear regression (HiCNorm).
#'@param subMatrix a subset of the HiC matrix (generally one chromosome).
#'@param subAnnotation a subset of the annotation corresponding to subMatrix.
#'@param useNegativeBinomial TRUE if using negative binomial model instead of Poisson.
#'@note internal function called by \code{\link{f.normalize.like.hu}} for intra-chromosomal normalization.
#'@references Hu, M. and Deng, K. and Selvaraj, S. and Qin, Z. S. and Ren, B. and Liu, J. S. (2012)
#'HiCNorm: removing biases in Hi-C data via Poisson regression. \emph{Bioinformatics} \bold{28}, 3131--3133.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.normalize.intra.like.hu <- function(subMatrix, subAnnotation, useNegativeBinomial) {
	# NOTE that this code is largely taken from HiCnorm: http://www.people.fas.harvard.edu/~junliu/HiCNorm/
	#change matrix into vector
	subVec <- subMatrix[upper.tri(subMatrix,diag=F)]
	#get cov matrix
	len_m<-as.matrix(log(subAnnotation[,1]%o%subAnnotation[,1]))
	gcc_m<-as.matrix(log(subAnnotation[,2]%o%subAnnotation[,2]))
	map_m<-as.matrix(log(subAnnotation[,3]%o%subAnnotation[,3]))
	#centralize cov matrix of enz, gcc
	len_m<-(len_m-mean(c(len_m)))/sd(c(len_m), na.rm = FALSE)
	gcc_m<-(gcc_m-mean(c(gcc_m)))/sd(c(gcc_m), na.rm = FALSE)
	#change matrix into vector
	len_vec<-len_m[upper.tri(len_m,diag=F)]
	gcc_vec<-gcc_m[upper.tri(gcc_m,diag=F)]
	map_vec<-map_m[upper.tri(map_m,diag=F)]
	# fit model
	if (useNegativeBinomial) {
		#require("MASS")
		fit<-MASS::glm.nb(subVec~len_vec+gcc_vec+offset(map_vec))
	} else {
		fit<-glm(subVec~len_vec+gcc_vec+offset(map_vec),family="poisson")
	}
	#summary(fit)
	coeff<-round(fit$coeff,4)
	out <- round(subMatrix/exp(coeff[1]+coeff[2]*len_m+coeff[3]*gcc_m+map_m), 4)
	return(out)
}

#'@title Normalize HiC data using linear regression (HiCNorm).
#'@param subMatrix a subset of the HiC matrix (generally two different chromosome).
#'@param subAnnotationA a subset of the annotation corresponding to chromosome A of subMatrix.
#'@param subAnnotationB a subset of the annotation corresponding to chromosome B of subMatrix.
#'@param useNegativeBinomial TRUE if using negative binomial model instead of Poisson.
#'@note internal function called by \code{\link{f.normalize.like.hu}} for inter-chromosomal normalization.
#'@references Hu, M. and Deng, K. and Selvaraj, S. and Qin, Z. S. and Ren, B. and Liu, J. S. (2012)
#'HiCNorm: removing biases in Hi-C data via Poisson regression. \emph{Bioinformatics} \bold{28}, 3131--3133.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.normalize.inter.like.hu <- function(subMatrix, subAnnotationA, subAnnotationB, useNegativeBinomial) {
	# NOTE that this code is largely taken from HiCnorm: http://www.people.fas.harvard.edu/~junliu/HiCNorm/
	#change matrix into vector
	subVec<-c(subMatrix)
	#get cov matrix
	len_m<-as.matrix(log(subAnnotationA[,1]%o%subAnnotationB[,1]))
	gcc_m<-as.matrix(log(subAnnotationA[,2]%o%subAnnotationB[,2]))
	map_m<-as.matrix(log(subAnnotationA[,3]%o%subAnnotationB[,3]))
	#centralize cov matrix of enz, gcc
	len_m<-(len_m-mean(c(len_m)))/sd(c(len_m))
	gcc_m<-(gcc_m-mean(c(gcc_m)))/sd(c(gcc_m))
	#change matrix into vector
	len_vec<-c(len_m)
	gcc_vec<-c(gcc_m)
	map_vec<-c(map_m)
	# fit model
	if (useNegativeBinomial) {
		#require("MASS")
		fit<-MASS::glm.nb(subVec~len_vec+gcc_vec+offset(map_vec))
	} else {
		fit<-glm(subVec~len_vec+gcc_vec+offset(map_vec),family="poisson")
	}
	#summary(fit)
	coeff<-round(fit$coeff,4)
	out<-round(subMatrix/exp(coeff[1]+coeff[2]*len_m+coeff[3]*gcc_m+map_m), 4)
	return(out)
}

#'@title Normalize HiC data for distance and coverage.
#'@param subMatrix a subset of the HiC matrix (generally one chromosome).
#'@note internal function called by \code{\link{f.normalize.data.matrix.like.LA}} and \code{\link{f.correlate.data.matrix}}.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.normalize.distance <- function(subMatrix) {
	m <- matrix(0, nrow = nrow(subMatrix), ncol = ncol(subMatrix))
	for (i in 1:(nrow(subMatrix)-1)) {
		band <- (row(m)==col(m)+i)
		m[band] <- mean(subMatrix[band])
	}
	m <- m + t(m)
	diag(m) <- mean(diag(subMatrix))
	out <- subMatrix - m
	return(out)
}

#'@title Normalize HiC data for distance and coverage.
#'@param dataMatrix a HiC interaction matrix (see \code{\link{f.load.one.sample}}).
#'@param binSize size of genomic bins in bp, if set to 0, restriction fragments are used instead.
#'@return a normalized HiC interaction matrix.
#'@references Liebermann-Aiden, E. and van Berkum, N. L. and Williams, L. and Imakaev, M. and Ragoczy, T. and Telling, A. and 
#'Amit, I. and Lajoie, B. R. and Sabo, P. J. and Dorschner, M. O. and Sandstrom, R. and Bernstein, B. and Bender, M. A. and 
#'Groudine, M. and Gnirke, A. and Stamatoyannopoulos, J. and Mirny, L. A. and Lander, E. S. and Dekker, J. (2009)
#'Comprehensive mapping of long-range interactions reveals folding principles of the human genome. \emph{Science} \bold{326}, 289--293.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@export
f.normalize.data.matrix.like.LA <- function(dataMatrix, binSize) {
	if (binSize == 0) { cat("this function does not work on individual fragments"); return(0);}
	dataMatrix <- log2(dataMatrix+1)
	seList <- f.get.se.list(binSize)
	chroms <- f.internal.get.relevant.chromosomes()
	for (i in 1:length(chroms)) {
		chromA <- chroms[i]
		sA <- seList[[chromA]][1]
		eA <- seList[[chromA]][2]
		for (j in i:length(chroms)) {
			chromB <- chroms[j]
			sB <- seList[[chromB]][1]
			eB <- seList[[chromB]][2]
			if (chromA == chromB) { 
				dataMatrix[sA:eA, sB:eB] <- f.internal.normalize.distance(dataMatrix[sA:eA, sB:eB])
			} else {
				dataMatrix[sA:eA, sB:eB] <- dataMatrix[sA:eA, sB:eB] - mean(dataMatrix[sA:eA, sB:eB])
				dataMatrix[sB:eB, sA:eA] <- dataMatrix[sB:eB, sA:eA] - mean(dataMatrix[sB:eB, sA:eA])
			}
		}
	}
	return(dataMatrix)
}

#'@title Correlate inter-chromosomal HiC interaction matrices.
#'@param subMatrix a subset of the HiC matrix (generally two different chromosome).
#'@note internal function called by \code{\link{f.correlate.data.matrix}}.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.correlate.inter <- function(subMatrix) {
	xMeans <- apply(subMatrix,1,mean)
	yMeans <- apply(subMatrix,2,mean)
	xMeans[which(xMeans==0)] <- mean(xMeans)
	yMeans[which(yMeans==0)] <- mean(yMeans)
	fhatNorm <- subMatrix/sqrt(xMeans %*% t(yMeans))
	totMean <- mean(as.vector(fhatNorm^(1/5)))
	devi <- (fhatNorm^(1/5)) - totMean
	divi <- max(abs(devi))
	out <- devi / divi
	return(out)
}

#'@title Correlate HiC interaction data.
#'@param dataMatrix a HiC interaction matrix (see \code{\link{f.load.one.sample}}).
#'@param binSize size of genomic bins in bp, if set to 0, restriction fragments are used instead.
#'@param useLog use logarithm instead of counts (default).
#'@return a matrix with correlated HiC interaction data.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@export
f.correlate.data.matrix <- function(dataMatrix, binSize, useLog = TRUE) {
	if (binSize == 0) { cat("this function does not work on individual fragments"); return(0);}
	if (useLog) {dataMatrix <- log2(dataMatrix+1)}
	seList <- f.get.se.list(binSize)
	chroms <- f.internal.get.relevant.chromosomes()
	for (i in 1:length(chroms)) {
		chromA <- chroms[i]
		sA <- seList[[chromA]][1]
		eA <- seList[[chromA]][2]
		for (j in i:length(chroms)) {
			chromB <- chroms[j]
			sB <- seList[[chromB]][1]
			eB <- seList[[chromB]][2]
			if (chromA == chromB) { 
				cm <- cor(f.internal.normalize.distance(dataMatrix[sA:eA, sB:eB]))
				diag(cm) <- mean(cm[(row(cm)==col(cm)+1)])
				dataMatrix[sA:eA, sB:eB] <- cm
			} else {
				dataMatrix[sA:eA, sB:eB] <- f.internal.correlate.inter(dataMatrix[sA:eA, sB:eB])
				dataMatrix[sB:eB, sA:eA] <- f.internal.correlate.inter(dataMatrix[sB:eB, sA:eA])
			}
		}
	}
	return(dataMatrix)
}

#'@title Correlate different HiC samples based on the pairwise correlation between their virtual 4C tracks.
#'@param dataMatrixA a HiC interaction matrix of sample A (see \code{\link{f.load.one.sample}}).
#'@param dataMatrixB a HiC interaction matrix of sample B (see \code{\link{f.load.one.sample}}).
#'@param corMethod method used to calculate correlation (pearson, spearman, kendall).
#'@param summaryFunction function used to summarize the individual 4C correlations (default is median).
#'@param useOnlyHighVar use only the virtual 4Cs with high variability (the ones which are in the top half in both samples).
#'Note that if this is turned off, the function removes 5 percent of the bins with the highest numbers of zeroes.
#'@return the correlation between the two samples (a numeric value).
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@export
f.HiC.correlation.via.fourC <- function(dataMatrixA, dataMatrixB, corMethod = "pearson", summaryFunction = median, useOnlyHighVar = TRUE) {
	if (useOnlyHighVar) {
		varA <- apply(dataMatrixA, 1, var)
		varB <- apply(dataMatrixB, 1, var)
		lbA <- quantile(varA, 0.5, na.rm = TRUE)
		lbB <- quantile(varB, 0.5, na.rm = TRUE)
		lvA <- varA < lbA
		lvB <- varB < lbB
		lvAB <- lvA | lvB
		indices <- 1:nrow(dataMatrixA)
		indices <- indices[!lvAB]
		cat(paste("removed", sum(lvAB, na.rm = TRUE), "values due to low variance\n"))
	} else {
		nzA <- apply(dataMatrixA == 0,1,sum)
		nzB <- apply(dataMatrixB == 0,1,sum)
		ubA <- quantile(nzA, 0.95, na.rm = TRUE)
		ubB <- quantile(nzB, 0.95, na.rm = TRUE)
		rmA <- nzA > ubA
		rmB <- nzB > ubB
		rmAB <- rmA | rmB
		indices <- 1:nrow(dataMatrixA)
		indices <- indices[!rmAB]
		cat(paste("removed", sum(rmAB, na.rm = TRUE), "values due to zeroes\n"))
	}
	temp <- rep(0, length(indices))
	for (i in 1:length(indices)) {
		ind <- indices[i]
		temp[i] <- cor(dataMatrixA[ind,], dataMatrixB[ind,], method = corMethod)
	}
	out <- summaryFunction(temp)
	return(out)
}

#'@title Get a set of colors for correlation heatmaps.
#'@param x number of colors.
#'@note internal function called by \code{\link{f.internal.blueredyellow}}.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.yellowredblue <- function(x) {
	r <- approx(c(0, 0.5, 1), c(1, 1, 0), n = x)$y
	g <- approx(c(0, 0.5, 1), c(1, 0, 0), n = x)$y
	b <- approx(c(0, 0.5, 1), c(0, 0, 1), n = x)$y
	return(rgb(r, g, b))
}

#'@title Get a set of colors for correlation heatmaps.
#'@param x number of colors.
#'@note internal function called by \code{\link{f.HiC.correlation.matrix}}.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.blueredyellow <- function(x) {
	return(rev(f.internal.yellowredblue(x)))
}

#'@title Draw a correlation matrix for multiple samples (correlation between multiple samples)
#'@param dataMatrixList a list of HiC matrices with the samplename as key to the entry (see \code{\link{f.load.samples}}).
#'@param rDir a directory where the figure is stored.
#'@param outfile name of the figure and the table (without file extension, _cor_mat.txt/.svg will be added).
#'@param corMethod method used to calculate correlation (pearson, spearman, kendall).
#'@param summaryFunction function used to summarize the individual 4C correlations (default is median).
#'@param useOnlyHighVar use only the virtual 4Cs with high variability (the ones which are in the top half in both samples).
#'Note that if this is turned off, the function removes 5 percent of the bins with the highest numbers of zeroes.
#'@return nothing (the figure is stored on the HD).
#'@examples \dontrun{
#'f.HiC.correlation.matrix(
#'  dataMatrixList = dataMatrices,
#'  rDir = "/path/to/where/the/figure/is/stored",
#'  outfile = "aNameForTheFigureWithoutExtension",
#'  corMethod = "pearson",
#'  summaryFunction = median,
#'  useOnlyHighVar = TRUE
#') }
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.HiC.correlation.via.fourC}}
#'@export
f.HiC.correlation.matrix <- function(dataMatrixList, rDir, outfile, corMethod = "pearson", summaryFunction = median, useOnlyHighVar = TRUE) {
	#require("gplots")
	samplesForCor <- names(dataMatrixList)
	sampleCorMat <- matrix(0, nrow = length(samplesForCor), ncol = length(samplesForCor), dimnames = list(samplesForCor, samplesForCor))
	for (i in 1:length(samplesForCor)) {
		sampleA <- samplesForCor[i]
		for (j in i:length(samplesForCor)) {
			sampleB <- samplesForCor[j]
			pCor <- f.HiC.correlation.via.fourC(log2(dataMatrixList[[sampleA]]+1), log2(dataMatrixList[[sampleB]]+1), corMethod, summaryFunction, useOnlyHighVar)
			sampleCorMat[i,j] <- pCor
			sampleCorMat[j,i] <- pCor
		}
	}
	write.csv(sampleCorMat, file.path(rDir, paste(outfile, '_cor_mat.txt', sep = '')))
	# some plotting parameters
	minCor <- min(sampleCorMat, na.rm = TRUE)
	maxCor <- max(sampleCorMat[sampleCorMat != 1], na.rm = TRUE)
	if (minCor > 0.8) { minCor <- 0.8 } else { if (minCor > 0.5) { minCor <- 0.5 }}
	if (maxCor > 0.8) { maxCor <- 1 }   else { if (maxCor < 0.8) { maxCor <- 0.8 } else { if (maxCor < 0.5) { maxCor <- 0.5 }}}
	numCols <- ceiling((maxCor-minCor)*50)
	if (GLOBAL_VARIABLE_USE_SVG) {
		svg(file.path(rDir, paste(outfile, '_cor_mat.svg', sep = '')), onefile = TRUE, bg = FALSE, antialias = "default", pointsize = 10, width = 12, height = 12)
	} else {
		tiff(file.path(rDir, paste(outfile, "_cor_mat.tiff", sep = '')), width = 1600, height = 1600)
	}
	if (nrow(sampleCorMat) > 6) {
		gplots::heatmap.2(sampleCorMat, col = f.internal.blueredyellow(numCols), trace="none", scale = "none", margins = c(15,15), cellnote = round(sampleCorMat,2), notecol = "black", notecex = 1, breaks = seq(minCor, maxCor, length.out = numCols+1))
	} else {
		gplots::heatmap.2(sampleCorMat, col = f.internal.blueredyellow(numCols), trace="none", scale = "none", margins = c(15,15), cellnote = round(sampleCorMat,2), notecol = "black", notecex = 3, breaks = seq(minCor, maxCor, length.out = numCols+1))
	}
	dev.off()
}

#'@title Find bins which are not mainly filled with zeroes.
#'@param dataMatrix a HiC interaction matrix (see \code{\link{f.load.one.sample}}).
#'@param threshold the percentile above which the bins are defined as "zero-bins".
#'@note internal function of HiCdatR.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.internal.find.zero.indices}}
f.internal.find.non.zero.indices <- function(dataMatrix, threshold = 0.95) {
	# make sure that it is a full dataMatrix and not just a subset
	zeros <- apply(dataMatrix == 0, 1, sum)
	upperBound <- quantile(zeros, threshold)
	valid <- which(zeros <= upperBound)
	return(valid)
}

#'@title Find bins which are mainly filled with zeroes.
#'@param dataMatrix a HiC interaction matrix (see \code{\link{f.load.one.sample}}).
#'@param threshold the percentile above which the bins are defined as "zero-bins".
#'@note internal function of HiCdatR.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.internal.find.non.zero.indices}}
f.internal.find.zero.indices <- function(dataMatrix, threshold = 0.95) {
	# make sure that it is a full dataMatrix and not just a subset
	zeros <- apply(dataMatrix == 0, 1, sum)
	upperBound <- quantile(zeros, threshold)
	invalid <- which(zeros > upperBound)
	return(invalid)
}

#'@title Check the annotation and remove non-relevat chromosomes and NAs.
#'@param annotation a dataframe holding the annotation of the bins (see \code{\link{f.read.annotation}}).
#'@param binSize size of genomic bins in bp, if set to 0, restriction fragments are used instead.
#'@note internal function called by \code{\link{f.read.annotation}} and \code{\link{f.read.annotation.via.fragment.annotation}}
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.strip.annotation <- function(annotation, binSize) {
	# keep only relevant chromosomes
	seList <- f.get.se.list.with.irrelevant.chromosomes(binSize)
	seListRelevant <- f.get.se.list(binSize)
	relevantChromosomes <- names(seListRelevant)[1:(length(seListRelevant)-1)]
	relevantBins <- c()
	for (rc in relevantChromosomes) { relevantBins <- c(relevantBins, seList[[rc]][1]:seList[[rc]][2]) }
	annotation <- annotation[relevantBins,]
	# check for na entires (low coverage DNA methylation can cause that)
	naSums <- apply(annotation, 2, function(x) sum(is.na(x)))
	if (sum(naSums) > 0) {
		haveNAs <- naSums[naSums > 0]
		cat("following columns have the displayed number of NAs which will be replaced by zero:\n")
		cat(paste(names(haveNAs), haveNAs, sep = ': ', collapse = '\n'))
		cat('\n')
		for (hasNA in names(haveNAs)) { annotation[is.na(annotation[,hasNA]),hasNA] <- 0 }
	}
	return(annotation)
}

#'@title Read the annotation file produced by HiCdatPre.
#'@param annotationFile a file holding restriction fragments or genomic bins and their annotation, see \url{https://github.com/MWSchmid/HiCdat}.
#'@param binSize size of genomic bins in bp, if set to 0, restriction fragments are used instead.
#'@param useLog if given, the sum_ and ann_ features are logged after the summarization/loading (default).
#'@return a table holding the annotation.
#'@examples \dontrun{
#'annotation <- f.read.annotation(
#'  annotationFile = "/path/to/file/created/with/HiCdatPre", 
#'  binSize = 1e5, 
#'  useLog = TRUE
#') }
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso for bins with fixed size, the annotation can also be obtained from an annotation obtained for restriction fragments or smaller bins
#'\code{\link{f.read.annotation.via.fragment.annotation}}.
#'@export
f.read.annotation <- function(annotationFile, binSize, useLog = TRUE) {
	out <- read.table(annotationFile, sep="\t", header=TRUE, stringsAsFactors = FALSE)
	out <- f.internal.strip.annotation(out, binSize)
	dataColsForSum <- grep("^sum_|^ann_", colnames(out), value = TRUE);
	if (useLog & (length(dataColsForSum) > 0)) { out[,dataColsForSum] <- log2(out[,dataColsForSum]+1); }
	return(out)
}

#'@title Read the annotation file produced by HiCdatPre.
#'@param annotationFile a file holding restriction fragments or genomic bins and their annotation, see \url{https://github.com/MWSchmid/HiCdat}.
#'@param binSize size of genomic bins in bp, if set to 0, restriction fragments are used instead.
#'@param useLog if given, the sum_ and ann_ features are logged after the summarization/loading (default).
#'@return a table holding the annotation.
#'@examples \dontrun{
#'annotation <- f.read.annotation.via.fragment.annotation(
#'  annotationFile = "/path/to/file/created/with/HiCdatPre", 
#'  binSize = 1e5, 
#'  useLog = TRUE
#') }
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.read.annotation}}
#'@export
f.read.annotation.via.fragment.annotation <- function(annotationFile, binSize, useLog = TRUE) {
	temp <- read.table(annotationFile, sep="\t", header=TRUE, stringsAsFactors = FALSE)
	relevantChromosomes <- f.internal.get.relevant.chromosomes()
	temp <- temp[temp$chrom %in% relevantChromosomes,] ## REMOVE IRRELEVANT CHROMOSOMES!
	# the annotation data is the table with the fragment wise annotations 
	# separate count and den features - the first are summed up, the latter are averaged
	dataColsForSum <- grep("^sum_|^ann_", colnames(temp), value = TRUE)
	dataColsForMean <- grep("^den_", colnames(temp), value = TRUE)
	# test if correct entries exist
	if ((length(dataColsForSum) == 0) & (length(dataColsForMean) == 0)) { cat("ERROR: there are no annotation features at all.\nMake sure that the columns holding annotation start with ann_, sum_, or den_.\n") }
	if (length(dataColsForSum) == 0) { cat("there are no features to sum up (i.e. no sum_* or ann_* features)\n") }
	if (length(dataColsForMean) == 0) { cat("there are no features to average (i.e. no den_* features)\n") }
	# translate chrom pos in matrix indices - use the start
	binNums <- f.translate.chrom.pos.vector.to.index(temp$chrom, temp$start, f.get.se.list.with.irrelevant.chromosomes(binSize), binSize)
	if (length(dataColsForSum) > 0) { winDataSum <- aggregate(temp[, dataColsForSum], by = list(winNum = binNums),  function(x) sum(x, na.rm = TRUE)) }
	if (length(dataColsForMean) > 0) { winDataMean <- aggregate(temp[, dataColsForMean], by = list(winNum = binNums),  function(x) mean(x, na.rm = TRUE)) }
	if ((length(dataColsForSum) > 0) & (length(dataColsForMean) > 0)) {
		out <- cbind(winDataSum[,dataColsForSum], winDataMean[,dataColsForMean])
	} else {
		if (length(dataColsForSum) > 0) {out <- winDataSum[,dataColsForSum]}
		else if (length(dataColsForMean) > 0) {out <- winDataMean[,dataColsForMean]}
		else {out <- NULL}
	}
	if (useLog & (length(dataColsForSum) > 0)) { out[,dataColsForSum] <- log2(out[,dataColsForSum]+1)}
	# add the columns fragmentNumber, chrom, start, end (to make it identical to the version where one loads the annotation directly)
	temp <- f.translate.index.vector.to.chrom.pos(1:nrow(out), f.get.se.list.with.irrelevant.chromosomes(binSize), binSize)
	out <- data.frame(fragmentNumber = 1:nrow(out), temp, out, stringsAsFactors = FALSE)
	out <- f.internal.strip.annotation(out, binSize)
	return(out)
}

#'@title Load one HiC interaction data sample.
#'@param dataDir the directory where the files are stored.
#'@param files a vector of files with HiC interaction tables produced by HiCdatPre, columns are: binA, binB, count.
#'All files will be read into one table and normalized together.
#'@param binSize size of genomic bins in bp, if set to 0, restriction fragments are used instead.
#'@param repetitions number of iterations for the normalization (see \code{\link{f.normalize.like.zhang}}).
#'If set to 0, normalization is omitted.
#'@param useLog if set to true, counts are logged (base 2, default off).
#'@return a HiC interaction matrix.
#'@examples \dontrun{
#'dataMatrix <- f.load.one.sample(
#' dataDir = "/path/to/files",
#' files = c("sampleA_run1.txt", "sampleA_run2.txt"),
#' binSize = 1e5,
#' repetitions = 50
#')}
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.load.samples}} and \code{\link{f.normalize.like.zhang}}.
#'@export
f.load.one.sample <- function(dataDir, files, binSize, repetitions = 0, useLog = FALSE) {
	seList <- f.get.se.list.with.irrelevant.chromosomes(binSize)
	numberInteractions <- 0
	out <- matrix(0, nrow = seList[["ALL"]][2], ncol = seList[["ALL"]][2])
	for (file in files) {
		cat(paste("loading ", file, "\n", sep = ''))
		temp <- read.table(file.path(dataDir, file), sep = '\t', col.names = c("binA", "binB", "count"))
		if ((min(temp$binA) < 1) | (min(temp$binB) < 1)) {
		  cat("WARNING: detected bins with an ID below 1 - assuming 0 and adding 1!\n")
		  temp$binA <- temp$binA+1
		  temp$binB <- temp$binB+1
		}
		if (useLog) { toAdd <- out[cbind(temp$binA,temp$binB)] + log2(temp$count+1) }
		else { toAdd <- out[cbind(temp$binA,temp$binB)] + temp$count }
		out[cbind(temp$binA,temp$binB)] <- toAdd
		out[cbind(temp$binB,temp$binA)] <- toAdd
		numberInteractions <- numberInteractions+sum(temp$count)
	}
	seListRelevant <- f.get.se.list(binSize) # remove now the irrelevant chromosomes
	relevantChromosomes <- names(seListRelevant)[1:(length(seListRelevant)-1)]
	relevantBins <- c()
	for (rc in relevantChromosomes) { relevantBins <- c(relevantBins, seList[[rc]][1]:seList[[rc]][2]) }
	out <- out[relevantBins, relevantBins]
	if (repetitions > 0) { out <- f.normalize.like.zhang(out, repetitions) }
	rm(temp)
	gc()
	return(out)
}

#'@title Load multiple HiC interaction data samples.
#'@param dataDir the directory where the files are stored.
#'@param sampleToFiles a list with sampleNames as keys to vectors holding the corresponding files.
#'eg: list(mySample = c("/home/myName/mySample_partA.txt", "/home/myName/mySample_partB.txt"))
#'@param binSize size of genomic bins in bp, if set to 0, restriction fragments are used instead.
#'@param repetitions number of iterations for the normalization (see \code{\link{f.normalize.like.zhang}}).
#'If set to 0, normalization is omitted.
#'@param useLog if set to true, counts are logged (base 2, default off).
#'@return a list holding HiC interaction matrix.
#'@examples \dontrun{
#'dataMatrices <- f.load.samples(
#'dataDir = "/path/to/files",
#'  sampleToFiles = list(
#'    sampleA = c("sampleA_run1.txt", "sampleA_run2.txt"), 
#'    sampleB = c("sampleB_run1.txt", "sampleB_run2.txt", "sampleB_run3.txt")
#'  ),
#'  binSize = 1e5,
#'  repetitions = 50
#') }
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.load.samples}} and \code{\link{f.normalize.like.zhang}}.
#'@export
f.load.samples <- function(dataDir, sampleToFiles, binSize, repetitions = 0, useLog = FALSE) {
	out <- list()
	for (sample in names(sampleToFiles)) {
		cat(paste("loading ", sample, "\n", sep = ''))
		files <- sampleToFiles[[sample]]
		out[[sample]] <- f.load.one.sample(dataDir, files, binSize, repetitions, useLog)
	}
	gc()
	return(out)
}

#'@title Create axis for the HiC plots.
#'@param binSize size of genomic bins in bp. 0 (restriction fragments) does not work.
#'@param axStep the size of one step on the axis in bp.
#'@param seList a list with first and last bins per chromosome (\code{\link{f.get.se.list}}).
#'@note internal function of HiCdatR.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.internal.axis.maker.on.index}}
f.internal.axis.maker <- function(binSize, axStep, seList){
	chroms <- names(f.get.se.list(binSize))
	out <- list()
	for (i in 1:length(chroms)) {
		if (i == 1) {
			axisVec <- seq(0,(max(seList[[chroms[i]]])*binSize), by=axStep)
			out[[chroms[i]]] <- axisVec
		} else {
			axisVec <- seq(0,(max(seList[[chroms[i]]])-max(seList[[chroms[i-1]]]))*binSize, by=axStep)
			out[[chroms[i]]] <- axisVec
		}
	}
	return(out)
}

#'@title Create axis for the HiC plots.
#'@param binSize size of genomic bins in bp. 0 (restriction fragments) does not work.
#'@param axStep the size of one step on the axis in bp.
#'@param seList a list with first and last bins per chromosome (\code{\link{f.get.se.list}}).
#'@note internal function of HiCdatR.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.internal.axis.maker}}
f.internal.axis.maker.on.index <- function(binSize, axStep, seList) {
	axisList <- f.internal.axis.maker(binSize, axStep, seList)
	out <- list()
	for (chrom in names(axisList)) {
		out[[chrom]] <- sapply(axisList[[chrom]], function(x) f.translate.chrom.pos.to.index(chrom, x, seList, binSize))
	}
	return(out)
}

#'@title Visualize HiC interaction data.
#'@param matrixToPlot a HiC interaction matrix (see \code{\link{f.load.one.sample}}).
#'@param binSize size of genomic bins in bp. 0 (restriction fragments) does not work.
#'@param axStep the size of one step on the axis in bp.
#'@param rDir a directory where the figure is stored.
#'@param outfile name of the figure and the table (without file extension).
#'@param chromA see \code{\link{f.extract.subset}}.
#'@param startA see \code{\link{f.extract.subset}}.
#'@param endA see \code{\link{f.extract.subset}}.
#'@param chromB see \code{\link{f.extract.subset}}.
#'@param startB see \code{\link{f.extract.subset}}.
#'@param endB see \code{\link{f.extract.subset}}.
#'@param useLog use logarithm instead of counts (default).
#'@param drawGrid draw a grid onto the plot (mainly useful for subregions)
#'@param doNorm normalize the data for distance and coverage (see \code{\link{f.normalize.data.matrix.like.LA}}).
#'@param doCor correlate the HiC matrix before drawing it (see \code{\link{f.correlate.data.matrix}}).
#'@param useSplineInterPol a plotting parameter (color interpolation)
#'@return noting, the figure is stored on the HD.
#'@examples \dontrun{
#'f.plot.XY.matrix(
#'  matrixToPlot = dataMatrix,
#'  binSize = 1e5,
#'  axStep = 10e6,
#'  rDir = "/path/to/where/the/figure/is/stored",
#'  outfile = "aNameForTheFigureWithoutExtension",
#'  chromA = "ALL",
#'  startA = 0,
#'  endA = 0,
#'  chromB = "ALL",
#'  startB = 0,
#'  endB = 0,
#'  useLog = TRUE,
#'  drawGrid = FALSE,
#'  doNorm = FALSE,
#'  doCor = FALSE, # or TRUE to draw a distance-normalized, correlated Hi-C-matrix
#'  useSplineInterPol = TRUE
#') }
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.load.one.sample}} and \code{\link{f.normalize.like.zhang}}.
#'@export
f.plot.XY.matrix <- function(matrixToPlot, binSize, axStep, rDir, outfile, chromA = "ALL", startA = 0, endA = 0, chromB = "ALL", startB = 0, endB = 0, useLog = TRUE, drawGrid = FALSE, doNorm = FALSE, doCor = FALSE, useSplineInterPol = TRUE) {
	# get seList and chromSizes
	seList <- f.get.se.list(binSize)
	chromSizes <- f.internal.get.chrom.sizes()
	if (startA < axStep) { startA <- 0 }
	if (startB < axStep) { startB <- 0 }
	# get the axes
	if (chromA == "ALL") {
		sa <- 1
		ea <- nrow(matrixToPlot)
		xAxAt <- as.vector(unlist(f.internal.axis.maker.on.index(binSize, axStep, seList)))
		xAxLab <- as.vector(unlist(f.internal.axis.maker(binSize,axStep,seList)))/1e6
		xGrid <- do.call("rbind", seList)[,1]
	} else {
		if (endA > chromSizes[[chromA]]) {endA <- chromSizes[[chromA]]}
		sa <- f.translate.chrom.pos.to.index(chromA, startA, seList, binSize)
		ea <- f.translate.chrom.pos.to.index(chromA, endA, seList, binSize)
		xAxAt <- seq(0, (ea-sa), by = axStep/binSize)
		xAxLab <- seq(startA, endA, by = axStep)/1e6
		xGrid <- xAxAt
	}
	if (chromB == "ALL") {
		sb <- 1
		eb <- ncol(matrixToPlot)
		yAxAt <- as.vector(unlist(f.internal.axis.maker.on.index(binSize, axStep, seList)))
		yAxLab <- as.vector(unlist(f.internal.axis.maker(binSize,axStep,seList)))/1e6
		yGrid <- do.call("rbind", seList)[,1]
	} else {
		if (endB > chromSizes[[chromB]]) {endB <- chromSizes[[chromB]]}
		sb <- f.translate.chrom.pos.to.index(chromB, startB, seList, binSize)
		eb <- f.translate.chrom.pos.to.index(chromB, endB, seList, binSize)
		yAxAt <- seq(0,(eb-sb), by = axStep/binSize)
		yAxLab <- seq(startB, endB, by = axStep)/1e6
		yGrid <- yAxAt
	}
	if (doCor) { matrixToPlot <- f.correlate.data.matrix(matrixToPlot, binSize, useLog) }
	else {
		if (doNorm) {matrixToPlot <- f.normalize.data.matrix.like.LA(matrixToPlot, binSize)}
		else { if (useLog) {matrixToPlot <- log2(matrixToPlot + 1)} }
	}
	matrixToPlot <- matrixToPlot[sa:ea,sb:eb]
	if (GLOBAL_VARIABLE_USE_SVG_AND_RSVG_CONVERT) {
		svg(file.path(rDir, paste(outfile, ".svg", sep = '')), height = 19, width = 19)
	} else {
		tiff(file.path(rDir, paste(outfile, ".tiff", sep = '')), width = 2400, height = 2400)
	}
	par(oma=c(5,5,0,0), mar = c(5,5,0,0))
	if (doCor) { 
		if (useSplineInterPol) {colorSet <- colorRampPalette(c("#ffeda0", "#feb24c","#f03b20"), interpolate="spline")(64)}
		else {colorSet <- colorRampPalette(c("#ffeda0", "#feb24c","#f03b20"))(64)}
		image(1:nrow(matrixToPlot),1:ncol(matrixToPlot),matrixToPlot, col = colorSet, useRaster=TRUE, yaxt = "n", xaxt = "n", xlab = "", ylab = "")
	}
	else { 
		if (useSplineInterPol) {colorSet <- colorRampPalette(c("black", "yellow","red"), interpolate="spline")(64)}
		else {colorSet <- colorRampPalette(c("black", "yellow","red"))(64)}
		image(1:nrow(matrixToPlot),1:ncol(matrixToPlot),matrixToPlot, col = colorSet, useRaster=TRUE, yaxt = "n", xaxt = "n", xlab = "", ylab = "")
	}
	axis(1, at = xAxAt, labels = xAxLab, outer=FALSE, line=2, lwd=2, cex.axis=1.5, las=2)
	axis(2, at = yAxAt, labels = yAxLab, outer=FALSE, line=2, lwd=2, cex.axis=1.5, las=1)
	if (drawGrid) {abline(h = yGrid, v = xGrid, lwd= 2, col= "white")}
	dev.off()
	if (GLOBAL_VARIABLE_USE_SVG_AND_RSVG_CONVERT) { system(paste("rsvg-convert -a -d 100 -p 100 ", file.path(rDir, paste(outfile, ".svg", sep = ''))," > ", file.path(rDir, paste(outfile, ".png", sep = '')), sep = '')) }
}

#'@title Randomize a HiC interaction matrix (for intrachromosomal map, the sampling is done only within a given distance).
#'@param dataMatrix a HiC interaction matrix (see \code{\link{f.load.one.sample}}).
#'@param binSize size of genomic bins in bp, if set to 0, restriction fragments are used instead.
#'@return a randomized HiC interaction matrix.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.randomize.matrix <- function(dataMatrix, binSize) {
	out <- matrix(0, nrow = nrow(dataMatrix), ncol = ncol(dataMatrix))
	chromosomes <- f.internal.get.relevant.chromosomes()
	seList <- f.get.se.list(binSize)
	for (i in 1:(length(chromosomes)-1)) {
		chromA <- chromosomes[i]
		for (j in (i+1):length(chromosomes)) {
			chromB <- chromosomes[j]
			if (chromA == chromB) {
				fromTo <- seList[[chromA]][1]:seList[[chromA]][2]
				temp <- dataMatrix[fromTo, fromTo]
				m <- matrix(0, nrow = nrow(temp), ncol = ncol(temp))
				for(i in 1:(nrow(temp)-1)) {
					band <- (row(m)==col(m)+i)
					m[band] <- sample(temp[band], length(temp[band]))
				}
				m <- m+t(m)
				diag(m) <- sample(diag(temp), length(diag(temp)))
				out[fromTo, fromTo] <- m
			} else {
				fromToA <- seList[[chromA]][1]:seList[[chromA]][2]
				fromToB <- seList[[chromB]][1]:seList[[chromB]][2]
				temp <- dataMatrix[fromToA, fromToB]
				randTemp <- sample(as.vector(temp))
				m <- matrix(randTemp, nrow = nrow(temp), ncol = ncol(temp))
				out[fromToA, fromToB] <- m
				out[fromToB, fromToA] <- t(m)
			}
		}
	}
	return(out)
}

#'@title Compare and plot two HiC samples using signed difference matrices (SDM).
#'@param dataMatrixA a HiC interaction matrix of sample A (see \code{\link{f.load.one.sample}}).
#'@param dataMatrixB a HiC interaction matrix of sample B (see \code{\link{f.load.one.sample}}).
#'@param binSize size of genomic bins in bp, if set to 0, restriction fragments are used instead.
#'@param rDir a directory where the figure is stored.
#'@param outfile name of the figure and the table (without file extension).
#'@param figureTitle the title of the figure.
#'@param filterZero whether or not to filter the x percent of bins with the highest number of 0 entries.
#'@param filterThreshold if filterZero is given, this specifies the fraction of bins which shall be kept.
#'@param pValueThreshold the threshold for the bin-wise pValues.
#'@param randomizeDiff randomize the matrices to test if the result is an artifact.
#'@return a list with the over all P-value and a vector with the significant columns.
#'@examples \dontrun{
#'SDMresult <- f.plot.signed.difference(
#'  dataMatrixA = dataMatrixSampleA,
#'  dataMatrixB = dataMatrixSampleB,
#'  binSize = 1e5,
#'  rDir = "/path/to/where/the/figure/is/stored",
#'  outfile = "aNameForTheFigureWithoutExtension",
#'  filterZero = TRUE,
#'  filterThreshold = 0.95,
#'  pValueThreshold = 0.01
#') }
#'@references Grob, S. and Schmid, M. W. and Grossniklaus, U. (2014) Hi-C analysis in Arabidopsis identifies
#'the KNOT, a structure with similarities to the flamenco locus of Drosophila. \emph{Molecular Cell} \bold{55}, 678--93.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.compare.samples.signed.difference}}
#'@export
f.plot.signed.difference <- function(dataMatrixA, dataMatrixB, binSize, rDir, outfile, figureTitle = "", filterZero = TRUE, filterThreshold = 0.95, pValueThreshold = 0.01, randomizeDiff = FALSE){
	#require("randomizeBE")
	# create the matrix with the signed differences
	if (randomizeDiff)	{ binDiff <- sign(f.internal.randomize.matrix(dataMatrixA - dataMatrixB, binSize)) }
	else 				{ binDiff <- sign(dataMatrixA - dataMatrixB) }
	# get the bins which are ok
	if (filterZero) 	{ toKeep <- intersect(f.internal.find.non.zero.indices(dataMatrixA, filterThreshold), f.internal.find.non.zero.indices(dataMatrixB, filterThreshold)) }
	else 				{ toKeep <- 1:nrow(binDiff) }
	# set the zero-filtered to zero (this is only for drawing):
	if (filterZero)		 { toRemove <- setdiff(1:nrow(binDiff), toKeep); binDiff[toRemove, ] <- 0; binDiff[, toRemove] <- 0;}
	# calculate p-values within rows
	pValues <- apply(binDiff, 1, function(x) randomizeBE::runs.pvalue(x[x!=0], pmethod="normal")) # take only the signs
	isSignificant <- pValues <= pValueThreshold
	pValuesSignificant <- intersect(which(isSignificant), toKeep)
	# calculate overall P-value - with runs-test (only - removed the empirical p-Value)
	overallPvalue <- randomizeBE::runs.pvalue(isSignificant[toKeep])
	# plot
	if (GLOBAL_VARIABLE_USE_SVG_AND_RSVG_CONVERT & (nrow(binDiff) < 1400)) {
		svg(file.path(rDir, paste(outfile, ".svg", sep = '')), width = 10, height = 11, pointsize = 12)
		cexTextTop <- 0.75
	} else {
		tiff(file.path(rDir, paste(outfile, ".tiff", sep = '')), width = 2400, height = 2500)
		cexTextTop <- 3
	}
	par(mar=c(0,0,0,0))
	layout(mat = c(1:2), widths = rep(1,2), heights = c(1,10))
	plot(pValuesSignificant, rep(1, length(pValuesSignificant)), xlim = c(0,length(pValues)), ylim = c(0,3), type = "h", xaxs = "i", xaxt = "n", yaxt = "n")
	legend("topleft", legend = paste("Number of significant stretches: ", length(pValuesSignificant), sep = ""), bty = "n", cex = cexTextTop)
	legend("topright", legend = figureTitle, bty = "n", cex = cexTextTop)
	legend("top", legend = c(paste("WW P: ",overallPvalue, sep="")), bty = "n", cex = cexTextTop)
	gridPoints <- do.call("rbind", f.get.se.list(binSize))[,1]
	image(x = c(0:nrow(binDiff)), y = c(0:nrow(binDiff)), z = binDiff, col = c("blue","white","red"), useRaster=TRUE) #col=c("#fc8d59","#ffffbf","#91bfdb")
	abline(h = gridPoints, v = gridPoints, lwd = 2, col = "black")
	dev.off()
	out <- list(overallPvalue = overallPvalue, significantRows = pValuesSignificant)
	return(out)
	if (GLOBAL_VARIABLE_USE_SVG_AND_RSVG_CONVERT) { system(paste("rsvg-convert -a -d 100 -p 100 ", file.path(rDir, paste(outfile, ".svg", sep = ''))," > ", file.path(rDir, paste(outfile, ".png", sep = '')), sep = '')) }
}

#'@title Perform multiple pairwise comparisons between HiC samples using signed difference matrices (SDM).
#'@param dataMatrixList a list of HiC matrices with the samplename as key to the entry (see \code{\link{f.load.samples}}).
#'@param binSize size of genomic bins in bp, if set to 0, restriction fragments are used instead.
#'@param rDir a directory where the figure is stored.
#'@param outfilePrefix a prefix for the figure and table files.
#'@param filterZero whether or not to filter the x percent of bins with the highest number of 0 entries.
#'@param filterThreshold if filterZero is given, this specifies the fraction of bins which shall be kept.
#'@param pValueThreshold the threshold for the bin-wise pValues.
#'@param randomizeDiff randomize the matrices to test if the result is an artifact.
#'@return for each pair a list with the over all P-value and a vector with the significant columns (see \code{\link{f.plot.signed.difference}}).
#'@examples \dontrun{
#'SDMresultList <- f.compare.samples.signed.difference(
#'  dataMatrixList = dataMatrices,
#'  binSize = 1e5,
#'  rDir = "/path/to/where/the/figures/are/stored",
#'  outfilePrefix = "aPrefixForTheFigureNames",
#'  filterZero = TRUE,
#'  filterThreshold = 0.95,
#'  pValueThreshold = 0.01
#'  ) }
#'@references Grob, S. and Schmid, M. W. and Grossniklaus, U. (2014) Hi-C analysis in Arabidopsis identifies
#'the KNOT, a structure with similarities to the flamenco locus of Drosophila. \emph{Molecular Cell} \bold{55}, 678--93.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.plot.signed.difference}}
#'@export
f.compare.samples.signed.difference <- function(dataMatrixList, binSize, rDir, outfilePrefix = "binDiff_", filterZero = TRUE, filterThreshold = 0.95, pValueThreshold = 0.01, randomizeDiff = FALSE) {
	samplesForDiff <-  names(dataMatrixList)
	out <- list()
	for (i in 1:(length(samplesForDiff)-1)) {
		sampleA <- samplesForDiff[i]
		for (j in (i+1):length(samplesForDiff)) {
			sampleB <- samplesForDiff[j]
			compStr <- paste(sampleA, sampleB, sep = "_vs_")
			cat(paste(compStr, '\n', sep = ''))
			out[[compStr]] <- f.plot.signed.difference(dataMatrixList[[sampleA]], dataMatrixList[[sampleB]], binSize, rDir, paste(outfilePrefix, compStr, sep = ''), compStr, filterZero, filterThreshold, pValueThreshold, randomizeDiff)
		}
	}
	return(out)
}

#'@title Compare and plot two HiC samples using relative differences.
#'@param dataMatrixA a HiC interaction matrix of sample A (see \code{\link{f.load.one.sample}}).
#'@param dataMatrixB a HiC interaction matrix of sample B (see \code{\link{f.load.one.sample}}).
#'@param binSize size of genomic bins in bp, if set to 0, restriction fragments are used instead.
#'@param rDir a directory where the figure is stored.
#'@param outfile name of the figure and the table (without file extension).
#'@param filterZero whether or not to filter the x percent of bins with the highest number of 0 entries.
#'@param filterThreshold if filterZero is given, this specifies the fraction of bins which shall be kept.
#'@param randomizeDiff randomize the matrices to test if the result is an artifact.
#'@return nothing, figures are stored on the HD.
#'@examples \dontrun{
#'f.plot.relative.difference(
#'  dataMatrixA = dataMatrixSampleA,
#'  dataMatrixB = dataMatrixSampleB,
#'  binSize = 1e5,
#'  rDir = "/path/to/where/the/figure/is/stored",
#'  outfile = "aNameForTheFigureWithoutExtension",
#'  filterZero = TRUE,
#'  filterThreshold = 0.95
#') }
#'@references Moissiard, G. and Cokus, S. J. and Cary, J. and Feng, S. and Billi, A. C. and Stroud, H. and
#'Husmann, D. and Zhan, Y. and Lajoie, B. R. and McCord, R. P. and Hale, C. J. and Feng, W. and Michaels, S. D. and
#'Frand, A. R. and Pellegrini, M. and Dekker, J. and Kim, J. K. and Jacobsen, S. E. (2012)
#'MORC family ATPases required for heterochromatin condensation and gene silencing. \emph{Science} \bold{336}, 1448--1451.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.compare.samples.relative.difference}}
#'@export
f.plot.relative.difference <- function(dataMatrixA, dataMatrixB, binSize, rDir, outfile, filterZero = TRUE, filterThreshold = 0.95, randomizeDiff = FALSE) {
	# create the matrix with the raw differences
	if (randomizeDiff)	{ rawDiff <- f.internal.randomize.matrix(dataMatrixA-dataMatrixB, binSize) }
	else 				{ rawDiff <- dataMatrixA-dataMatrixB }
	# get the bins which are ok
	if (filterZero) 	{ toKeep <- intersect(f.internal.find.non.zero.indices(dataMatrixA, filterThreshold), f.internal.find.non.zero.indices(dataMatrixB, filterThreshold)) }
	else 				{ toKeep <- 1:nrow(rawDiff) }
	# set the zero-filtered to zero (this is only for drawing):
	if (filterZero)		{ toRemove <- setdiff(1:nrow(rawDiff), toKeep); rawDiff[toRemove, ] <- 0; rawDiff[, toRemove] <- 0;}
	# normalize the difference
	normDiff <- rawDiff/((dataMatrixA+dataMatrixB)/2)
	normDiff[is.na(normDiff)] <- 0
	# plotting parameters
	seList <- f.get.se.list(binSize)
	axAt <- as.vector(unlist(f.internal.axis.maker.on.index(binSize, 100*binSize, seList)))
	axLab <- as.vector(unlist(f.internal.axis.maker(binSize, 100*binSize, seList)))/1e6
	grid <- do.call("rbind", seList)[,1]
	if (GLOBAL_VARIABLE_USE_SVG_AND_RSVG_CONVERT) {
		svg(file.path(rDir, paste(outfile, ".svg", sep = '')), height = 19, width = 19)
	} else {
		tiff(file.path(rDir, paste(outfile, ".tiff", sep = '')), width = 2400, height = 2400)
	}
	par(oma=c(5,5,0,0), mar = c(5,5,0,0))
	image(1:nrow(normDiff),1:ncol(normDiff), normDiff, col = colorRampPalette(c("blue", "white","red"))(64), useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = "", ylab = "")
	axis(1, at = axAt, labels = axLab, outer = FALSE, line = 2, lwd = 2, cex.axis = 1.5, las = 2)
	axis(2, at = axAt, labels = axLab, outer = FALSE, line = 2, lwd = 2, cex.axis = 1.5, las = 1)
	abline(h = grid, v = grid, lwd = 2, col = "black")
	dev.off()
	if (GLOBAL_VARIABLE_USE_SVG_AND_RSVG_CONVERT) { system(paste("rsvg-convert -a -d 100 -p 100 ", file.path(rDir, paste(outfile, ".svg", sep = ''))," > ", file.path(rDir, paste(outfile, ".png", sep = '')), sep = '')) }
}

#'@title Perform multiple pairwise comparisons between HiC samples using relative differences.
#'@param dataMatrixList a list of HiC matrices with the samplename as key to the entry (see \code{\link{f.load.samples}}).
#'@param binSize size of genomic bins in bp, if set to 0, restriction fragments are used instead.
#'@param rDir a directory where the figure is stored.
#'@param outfilePrefix a prefix for the figure and table files.
#'@param filterZero whether or not to filter the x percent of bins with the highest number of 0 entries.
#'@param filterThreshold if filterZero is given, this specifies the fraction of bins which shall be kept.
#'@param randomizeDiff randomize the matrices to test if the result is an artifact.
#'@return nothing, figures are stored on the HD.
#'@examples \dontrun{
#'f.compare.samples.relative.difference(
#'  dataMatrixList = dataMatrices,
#'  binSize = 1e5,
#'  rDir = "/path/to/where/the/figures/are/stored",
#'  outfilePrefix = "aPrefixForTheFileNames",
#'  filterZero = TRUE,
#'  filterThreshold = 0.95
#') }
#'@references Moissiard, G. and Cokus, S. J. and Cary, J. and Feng, S. and Billi, A. C. and Stroud, H. and
#'Husmann, D. and Zhan, Y. and Lajoie, B. R. and McCord, R. P. and Hale, C. J. and Feng, W. and Michaels, S. D. and
#'Frand, A. R. and Pellegrini, M. and Dekker, J. and Kim, J. K. and Jacobsen, S. E. (2012)
#'MORC family ATPases required for heterochromatin condensation and gene silencing. \emph{Science} \bold{336}, 1448--1451.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.plot.relative.difference}}
#'@export
f.compare.samples.relative.difference <- function(dataMatrixList, binSize, rDir, outfilePrefix = "relDiff_", filterZero = TRUE, filterThreshold = 0.95, randomizeDiff = FALSE) {
	samplesForDiff <-  names(dataMatrixList)
	out <- list()
	for (i in 1:(length(samplesForDiff)-1)) {
		sampleA <- samplesForDiff[i]
		for (j in (i+1):length(samplesForDiff)) {
			sampleB <- samplesForDiff[j]
			compStr <- paste(sampleA, sampleB, sep = "_vs_")
			cat(paste(compStr, '\n', sep = ''))
			out[[compStr]] <- f.plot.relative.difference(dataMatrixList[[sampleA]], dataMatrixList[[sampleB]], binSize, rDir, paste(outfilePrefix, compStr, sep = ''), filterZero, filterThreshold, randomizeDiff)
		}
	}
	#return(out)
}

#'@title Compare and plot two HiC samples using correlated differences.
#'@param dataMatrixA a HiC interaction matrix of sample A (see \code{\link{f.load.one.sample}}).
#'@param dataMatrixB a HiC interaction matrix of sample B (see \code{\link{f.load.one.sample}}).
#'@param binSize size of genomic bins in bp, if set to 0, restriction fragments are used instead.
#'@param rDir a directory where the figure is stored.
#'@param outfile name of the figure and the table (without file extension).
#'@param filterZero whether or not to filter the x percent of bins with the highest number of 0 entries.
#'@param filterThreshold if filterZero is given, this specifies the fraction of bins which shall be kept.
#'@param randomizeDiff randomize the matrices to test if the result is an artifact.
#'@return nothing, figures are stored on the HD.
#'@examples \dontrun{
#'f.plot.cor.difference(
#'  dataMatrixA = dataMatrixSampleA,
#'  dataMatrixB = dataMatrixSampleB,
#'  binSize = 1e5,
#'  rDir = "/path/to/where/the/figure/is/stored",
#'  outfile = "aNameForTheFigureWithoutExtension",
#'  filterZero = TRUE,
#'  filterThreshold = 0.95
#') }
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.compare.samples.relative.difference}}
#'@export
f.plot.cor.difference <- function(dataMatrixA, dataMatrixB, binSize, rDir, outfile, filterZero = TRUE, filterThreshold = 0.95, randomizeDiff = FALSE) {
	# create the matrix with the raw differences
	if (randomizeDiff)	{ rawDiff <- f.internal.randomize.matrix(dataMatrixA-dataMatrixB, binSize) }
	else 				{ rawDiff <- dataMatrixA-dataMatrixB }
	# get the bins which are ok
	if (filterZero) 	{ toKeep <- intersect(f.internal.find.non.zero.indices(dataMatrixA, filterThreshold), f.internal.find.non.zero.indices(dataMatrixB, filterThreshold)) }
	else 				{ toKeep <- 1:nrow(rawDiff) }
	# set the zero-filtered to zero (this is only for drawing):
	if (filterZero)		{ toRemove <- setdiff(1:nrow(rawDiff), toKeep); rawDiff[toRemove, ] <- 0; rawDiff[, toRemove] <- 0;}
	# normalize the difference
	normDiff <- rawDiff/((dataMatrixA+dataMatrixB)/2)
	normDiff[is.na(normDiff)] <- 0
	# correlate
	corDiff <- cor(normDiff)
	if (filterZero)		{ corDiff[toRemove, ] <- 0; corDiff[, toRemove] <- 0;}
	# plotting parameters
	seList <- f.get.se.list(binSize)
	axAt <- as.vector(unlist(f.internal.axis.maker.on.index(binSize, 100*binSize, seList)))
	axLab <- as.vector(unlist(f.internal.axis.maker(binSize, 100*binSize, seList)))/1e6
	grid <- do.call("rbind", seList)[,1]
	if (GLOBAL_VARIABLE_USE_SVG_AND_RSVG_CONVERT) {
		svg(file.path(rDir, paste(outfile, ".svg", sep = '')), height = 19, width = 19)
	} else {
		tiff(file.path(rDir, paste(outfile, ".tiff", sep = '')), width = 2400, height = 2400)
	}
	par(oma=c(5,5,0,0), mar = c(5,5,0,0))
	image(1:nrow(corDiff),1:ncol(corDiff), corDiff, col = colorRampPalette(c("#d8b365", "white", "#5ab4ac"))(9), zlim=c(-0.5,0.5), useRaster = TRUE, yaxt = "n", xaxt = "n", xlab = "", ylab = "")
	axis(1, at = axAt, labels = axLab, outer = FALSE, line = 2, lwd = 2, cex.axis = 1.5, las = 2)
	axis(2, at = axAt, labels = axLab, outer = FALSE, line = 2, lwd = 2, cex.axis = 1.5, las = 1)
	abline(h = grid, v = grid, lwd = 2, col = "black")
	dev.off()
	if (GLOBAL_VARIABLE_USE_SVG_AND_RSVG_CONVERT) { system(paste("rsvg-convert -a -d 100 -p 100 ", file.path(rDir, paste(outfile, ".svg", sep = ''))," > ", file.path(rDir, paste(outfile, ".png", sep = '')), sep = '')) }
}

#'@title Perform multiple pairwise comparisons between HiC samples using correlated differences.
#'@param dataMatrixList a list of HiC matrices with the samplename as key to the entry (see \code{\link{f.load.samples}}).
#'@param binSize size of genomic bins in bp, if set to 0, restriction fragments are used instead.
#'@param rDir a directory where the figure is stored.
#'@param outfilePrefix a prefix for the figure and table files.
#'@param filterZero whether or not to filter the x percent of bins with the highest number of 0 entries.
#'@param filterThreshold if filterZero is given, this specifies the fraction of bins which shall be kept.
#'@param randomizeDiff randomize the matrices to test if the result is an artifact.
#'@return nothing, figures are stored on the HD.
#'@examples \dontrun{
#'f.compare.samples.cor.difference(
#'  dataMatrixList = dataMatrices,
#'  binSize = 1e5,
#'  rDir = "/path/to/where/the/figures/are/stored",
#'  outfilePrefix = "aPrefixForTheFileNames",
#'  filterZero = TRUE,
#'  filterThreshold = 0.95
#') }
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.plot.relative.difference}}
#'@export
f.compare.samples.cor.difference <- function(dataMatrixList, binSize, rDir, outfilePrefix = "corDiff_", filterZero = TRUE, filterThreshold = 0.95, randomizeDiff = FALSE) {
	samplesForDiff <-  names(dataMatrixList)
	out <- list()
	for (i in 1:(length(samplesForDiff)-1)) {
		sampleA <- samplesForDiff[i]
		for (j in (i+1):length(samplesForDiff)) {
			sampleB <- samplesForDiff[j]
			compStr <- paste(sampleA, sampleB, sep = "_vs_")
			cat(paste(compStr, '\n', sep = ''))
			out[[compStr]] <- f.plot.cor.difference(dataMatrixList[[sampleA]], dataMatrixList[[sampleB]], binSize, rDir, paste(outfilePrefix, compStr, sep = ''), filterZero, filterThreshold, randomizeDiff)
		}
	}
	#return(out)
}

#'@title Calculate and plot the interaction decay exponent (see reference below for details).
#'@param dataMatrix a HiC interaction matrix (see \code{\link{f.load.one.sample}}).
#'@param binSize size of genomic bins in bp. 0 (restriction fragments) does not work.
#'@param rDir a directory where the figure is stored.
#'@param outfile name of the figure and the table (without file extension).
#'@param distance specifies the distance up to which the interaction frequency decay shall be considered.
#'Must be smaller than the smallest region (e.g. smaller than the smalles chromosome if \code{regionTable = data.frame()})
#'@param regionTable defines specific regions in the genome which shall be assessed. Per default,
#'each chromosome is first tested separately and a common IDE is calculated as the average between
#'all individual IDEs. The table must have three columns (chrom, start, and end). User defined names
#'can be given using an optional fourth column (name). An example is given in the tutorial \url{https://github.com/MWSchmid/HiCdat}.
#'@param filterZero whether or not to filter the x percent of bins with the highest number of 0 entries.
#'@param filterThreshold if filterZero is given, this specifies the fraction of bins which shall be kept.
#'@return a list with the results of the fit and the table with the individual values.
#'@examples \dontrun{
#'f.distance.decay(
#'  dataMatrix = dataMatrixSampleX,
#'  binSize = 1e5,
#'  rDir = "/path/to/where/the/figure/is/stored",
#'  outfile = "aNameForTheFigureWithoutExtension",
#'  distance = 10e6,
#'  regionTable = data.frame(),
#'  filterZero = TRUE,
#'  filterThreshold = 0.95
#') }
#'@references Liebermann-Aiden, E. and van Berkum, N. L. and Williams, L. and Imakaev, M. and Ragoczy, T. and Telling, A. and 
#'Amit, I. and Lajoie, B. R. and Sabo, P. J. and Dorschner, M. O. and Sandstrom, R. and Bernstein, B. and Bender, M. A. and 
#'Groudine, M. and Gnirke, A. and Stamatoyannopoulos, J. and Mirny, L. A. and Lander, E. S. and Dekker, J. (2009)
#'Comprehensive mapping of long-range interactions reveals folding principles of the human genome. \emph{Science} \bold{326}, 289--293.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@export
f.distance.decay <- function(dataMatrix, binSize, rDir, outfile, distance, regionTable = data.frame(), filterZero = TRUE, filterThreshold = 0.95) {
	if (nrow(regionTable) == 0) { chromSizes <- f.get.relevant.chrom.sizes(); regionTable <- data.frame(chrom = names(chromSizes), start = rep(0, length(chromSizes)), end = unlist(chromSizes), stringsAsFactors = FALSE) }
	if (ncol(regionTable) == 3) { regionTable$name <- paste(regionTable$chrom, regionTable$start, regionTable$end, sep = '_') }
	numRegs <- nrow(regionTable)
	if (filterZero) { zeroes <- f.internal.find.zero.indices(dataMatrix, filterThreshold); dataMatrix[zeroes,] <- NA; dataMatrix[,zeroes] <- NA;}
	seList <- f.get.se.list(binSize)
	chromSizes <- f.internal.get.chrom.sizes()
	binDistance <- 2:floor(distance/binSize)
	distVecList <- list()
	lineCols <- c(colorRampPalette(c("#d8b365", "#5ab4ac"))(numRegs), "black")
	names(lineCols) <- c(regionTable$name, "average")
	resultMatrix <- matrix(0, nrow = numRegs, ncol = 3, dimnames = list(regionTable$name, c("intercept", "slope", "pValue")))
	for (i in 1:numRegs) {
		chrom <- regionTable$chrom[i]
		start <- regionTable$start[i]
		end <- regionTable$end[i]
		regName <- regionTable$name[i]
		if (end > chromSizes[[chrom]]) {end <- chromSizes[[chrom]]}
		startBin <- f.translate.chrom.pos.to.index(chrom, start, seList, binSize)
		endBin <- f.translate.chrom.pos.to.index(chrom, end, seList, binSize)
		if (max(binDistance) > (endBin-startBin)) {cat(paste("ERROR - INTERRUPTED - region ", regName, " is too small - reduce the distance or remove this region\n", sep = '')); print(startBin);print(endBin); return(NULL);}
		z <- dataMatrix[startBin:endBin, startBin:endBin]
		m <- matrix(0, nrow = nrow(z), ncol = ncol(z))
		for(j in 1:(nrow(z)-1)) {
			band <- (row(m)==col(m)+j)
			m[band] <- mean(z[band], na.rm = TRUE)
		}
		diag(m) <- mean(diag(z), na.rm = TRUE)
		distVec <- m[binDistance,1]
		if(sum(is.na(distVec)) > 0) { cat(paste("removing additional NA values on chromosome", chrom, "- be careful with the results for this chromosome\n")); distVec[is.na(distVec)] <- 0 }
		distVecForTest <- distVec[distVec != 0]
		distanceForTest <- binDistance[distVec != 0]
		linMod <- lm(log10(distVecForTest) ~ log10(distanceForTest))
		distVecList[[regName]] <- distVec
		resultMatrix[regName, "intercept"] <- linMod$coefficients[1]
		resultMatrix[regName, "slope"] <- linMod$coefficients[2]
		resultMatrix[regName, "pValue"] <- anova(linMod)[["Pr(>F)"]][1]
	}
	vecTab <- do.call("cbind", distVecList)
	colnames(vecTab) <- regionTable$name
	# get the average slope
	meanVec <- apply(vecTab,1, function(x) mean(x, na.rm = TRUE))
	fullTab <- cbind()
	linMod <- lm(log10(meanVec) ~ log10(binDistance))
	pVal <- anova(linMod)[["Pr(>F)"]][1]
	icept <- linMod$coefficients[1]
	slope <- linMod$coefficients[2]
	# add the average to the other tables
	fullTab <- cbind(vecTab, meanVec)
	colnames(fullTab) <- c(regionTable$name, "average")
	resultMatrix <- rbind(resultMatrix, c(icept, slope, pVal))
	rownames(resultMatrix)[nrow(resultMatrix)] <- "average"
	# plotting
	if (GLOBAL_VARIABLE_USE_SVG) {
		svg(file.path(rDir, paste(outfile, ".svg",sep="")), width = 30, height = 15, pointsize = 25)
	} else {
		tiff(file.path(rDir, paste(outfile, ".tiff", sep = '')), width = 1600, height = 1600)
	}
	emptyPlot <- fullTab[,1]
	pointsToPlot <- (emptyPlot != 0)
	layout(matrix(c(1,2), nrow = 1))
	plot(log10(binDistance[pointsToPlot]), log10(emptyPlot[pointsToPlot]), bty = "n", main = "", ylab = "log10(contact frequency)", xlab = "log10(distance)", type = "n", ylim = c(-0.5, 2), xaxt = "n")
	for (regName in c(regionTable$name, "average")){
		toPlot <- fullTab[,regName]
		pointsToPlot <- (toPlot != 0)
		points(log10(binDistance[pointsToPlot]), log10(toPlot[pointsToPlot]), col = lineCols[regName], cex = 0.75, pch = 16)
		abline(resultMatrix[regName, "intercept"], resultMatrix[regName, "slope"], col = lineCols[regName], lwd = 1)
	}
	axis(1, at = seq(log10(2), log10(floor(distance/binSize)), length.out = 10), labels = round(seq(2, log10(distance), length.out = 10),1))
	plot(log10(binDistance[pointsToPlot]), log10(emptyPlot[pointsToPlot]), bty = "n", main = "", ylab = "", xlab = "", type = "n", ylim = c(-0.5, 2), xaxt = "n", yaxt = "n")
	legend("left", legend = paste(rownames(resultMatrix), "slope:", round(resultMatrix[,"slope"], 2)), col = lineCols, pch = 16, bty = "n")
	dev.off()
	out <- list(modelValues = resultMatrix, plotEntries = fullTab)
	return(out)
}

#'@title Perform principle component analysis on HiC data.
#'@param dataMatrix a HiC interaction matrix (see \code{\link{f.load.one.sample}}).
#'@param binSize size of genomic bins in bp. 0 (restriction fragments) does not work.
#'@param rDir a directory where the figure is stored.
#'@param outfile name of the figure (without file extension).
#'@param regionTable defines specific regions in the genome which shall be assessed. Per default,
#'each chromosome is tested separately. The table must have three columns (chrom, start, and end). User defined names
#'can be given using an optional fourth column (name). An example is given in the tutorial \url{https://github.com/MWSchmid/HiCdat}.
#'@param filterZero whether or not to filter the x percent of bins with the highest number of 0 entries.
#'@param filterThreshold if filterZero is given, this specifies the fraction of bins which shall be kept.
#'@param userLimits define your own min/max for the correlation values (related to the plot and its color mapping).
#'@note internal function called by \code{\link{f.principle.component.analysis.and.features}}.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.principle.component.analysis <- function(dataMatrix, binSize, rDir, outfile, regionTable = data.frame(), filterZero = TRUE, filterThreshold = 0.95, userLimits = c(-1, 1)) {
	if (nrow(regionTable) == 0) { chromSizes <- f.get.relevant.chrom.sizes(); regionTable <- data.frame(chrom = names(chromSizes), start = rep(0, length(chromSizes)), end = unlist(chromSizes), stringsAsFactors = FALSE) }
	if (ncol(regionTable) == 3) { regionTable$name <- paste(regionTable$chrom, regionTable$start, regionTable$end, sep = '_') }
	numRegs <- nrow(regionTable)
	if (filterZero) { binsToUse <- f.internal.find.non.zero.indices(dataMatrix, filterThreshold) }
	else { binsToUse <- 1:nrow(dataMatrix) }
	seList <- f.get.se.list(binSize)
	chromSizes <- f.internal.get.chrom.sizes()
	out <- list()
	if (GLOBAL_VARIABLE_USE_SVG) {
		svg(file.path(rDir, paste(outfile, ".svg", sep = '')), height = 9, width = 3*numRegs)
	} else {
		tiff(file.path(rDir, paste(outfile, ".tiff", sep = '')), width = 300*numRegs, height = 900)
	}
	layout(mat = matrix(c(1:(3*numRegs)), nrow = 3))
	par(oma = c(0,0,0,0), mar = c(0,0,0,0))
	for (i in 1:numRegs) {
		chrom <- regionTable$chrom[i]
		start <- regionTable$start[i]
		end <- regionTable$end[i]
		regName <- regionTable$name[i]
		if (end > chromSizes[[chrom]]) {end <- chromSizes[[chrom]]}
		startBin <- f.translate.chrom.pos.to.index(chrom, start, seList, binSize)
		endBin <- f.translate.chrom.pos.to.index(chrom, end, seList, binSize)
		validBins <- intersect(startBin:endBin, binsToUse)
		# normalize, correlate and then calculate the PCA
		z <- dataMatrix[validBins, validBins]
		m <- matrix(0, nrow=nrow(z), ncol=ncol(z))
		for(i in 1:(nrow(z)-1)) {
			band <- (row(m) == col(m)+i)
			if (mean(z[band]) == 0) { m[band] <- 1 }
			else { m[band] <- mean(z[band]) }
		}
		m <- m+t(m)
		if (mean(diag(z)) == 0) { diag(m) <- 1;} 
		else { diag(m) <- mean(diag(z)) }
		zNorm <- z/m
		CORz <- cor(zNorm)
		CORzForPlot <- matrix(0, nrow = length(startBin:endBin), ncol = length(startBin:endBin))
		CORzForPlot[validBins-startBin+1,validBins-startBin+1] <- CORz
		princp <- princomp(CORz)
		## TODO how should one turn the princomp. For AT, the first approach worked better, for mouse the second (where it is most likely due to the mappability)
# 		# turn the PCA in a way that negative means closed (on average more interaction with neighboring regions)
		orientedPrincp <- princp$loadings[,1]
		originalPlus <- which(orientedPrincp > 0)
		originalMinus <- which(orientedPrincp < 0)
		toCheck <- f.summary.along.diagonal(zNorm, 3, sum, FALSE)
		plusMean <- mean(toCheck[originalPlus])
		minusMean <- mean(toCheck[originalMinus])
		if (plusMean > minusMean) {
			cat(paste(outfile, " - switching sign of ", regName, " (plus/minus): ", plusMean, " / ", minusMean, '\n', sep =''))
			orientedPrincp <- -orientedPrincp
		}
# 		# well - try once to turn it according to the total number of interactions... (row/colsums)
# 		orientedPrincp <- princp$loadings[,1]
# 		originalPlus <- which(orientedPrincp > 0)
# 		originalMinus <- which(orientedPrincp < 0)
# 		toCheck <- apply(zNorm, 2, sum)
# 		plusMean <- mean(toCheck[originalPlus])
# 		minusMean <- mean(toCheck[originalMinus])
# 		if (plusMean < minusMean) {
# 			cat(paste(outfile, " - switching sign of ", regName, " (plus/minus): ", plusMean, " / ", minusMean, '\n', sep =''))
# 			orientedPrincp <- -orientedPrincp
# 		}
		# plot correlation matrix, first principle component, and a histrogramm with the correlation distribution to make sure the orientation is correct
		# cor matrix
		image(CORzForPlot, col = colorRampPalette(c("#ffeda0", "#feb24c","#f03b20"), interpolate="spline")(64), yaxt = "n", xaxt = "n", xlab = "", ylab = "", zlim = userLimits)
		# first principle component
		plot(validBins, orientedPrincp, type = 'l', xaxt = "n", yaxt = "n", main = "", xlab = "", ylab = "", xaxs = "i")
		abline(h=0)
		# histogram
		par(mar = c(5,4,1,1))
		plot(density(as.vector(CORz)), main = "", bty = "n", xaxs = "r", yaxs = "r", xlab = regName, ylab = "", las = 1, cex = 0.4, tck = 0.01)
		text(par("usr")[1] + par("usr")[2]/15, par("usr")[4] - par("usr")[4]/15, adj = c(0,1), labels = paste(c("median:", "mean:", "sd:"), c(round(median(as.vector(CORz)), digits = 3), round(mean(as.vector(CORz)), digits = 3), round(sd(as.vector(CORz)), digits = 3)), collapse = "\n", sep=" "))
		abline(v=0)
		par(mar = c(0,0,0,0))
		# store for return
		out[[regName]] <- cbind(validBins, orientedPrincp)
	}
	dev.off()
	return(out)
}

#'@title Correlate the first principle component to genomic/epigenetic features.
#'@param pca a list holding tables with valid bins an their value of the first principle component (see \code{\link{f.internal.principle.component.analysis}}).
#'@param annotation a dataframe holding the annotation of the bins (see \code{\link{f.read.annotation}}).
#'@param featureNames names of the features which are drawn (used to simplify the annotation names).
#'@param rDir a directory where the figure is stored.
#'@param outfilePrefix a prefix for the figure and table files.
#'@param pValueThreshold specifies the significance-threshold (only significant values are drawn in the heatmaps).
#'@note internal function called by \code{\link{f.principle.component.analysis.and.features}}.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.pca.and.feature.correlation.test <- function(pca, annotation, featureNames, rDir, outfilePrefix, pValueThreshold) {
	# pca: result from f.principle.component.analysis
	# annotation: the table holding fragment/bin annotations
	# featureNames: simplified feature names
	regNames <- names(pca)
	featuresToTest <- grep("^ann_|^sum_|^den_", colnames(annotation), value = TRUE)
	estimates <- matrix(0, nrow = length(regNames), ncol = length(featuresToTest), dimnames = list(regNames, featuresToTest))
	pValues <- matrix(0, nrow = length(regNames), ncol = length(featuresToTest), dimnames = list(regNames, featuresToTest))
	for (regName in regNames) {
		validBins <- pca[[regName]][,"validBins"]
		orientedPrincp <- pca[[regName]][,"orientedPrincp"]
		validAnnotation <- annotation[validBins,]
		for (feature in featuresToTest){
			temp <- cor.test(orientedPrincp,validAnnotation[,feature])
			estimates[regName, feature] <- temp$estimate
			pValues[regName, feature] <- temp$p.value
		}
	}
	# replace NAs and names
	estimates[is.na(estimates)] <- 0
	pValues[is.na(pValues)] <- 1
	colnames(estimates) <- featureNames
	colnames(pValues) <- featureNames
	# plot a heatmap with the significant fields
	forPlot <- estimates
	forPlot[pValues>pValueThreshold] <- 0
	xAxAt <- 1:nrow(forPlot)
	yAxAt <- 1:ncol(forPlot)
	xAxLab <- rownames(forPlot)
	yAxLab <- colnames(forPlot)
	if (GLOBAL_VARIABLE_USE_SVG) {
		svg(file.path(rDir, paste(outfilePrefix, "_PCA_significant_correlation.svg", sep = '')), height = 20, width = 20)
	} else {
		tiff(file.path(rDir, paste(outfilePrefix, "_PCA_significant_correlation.tiff", sep = '')), width = 1600, height = 1600)
	}
	par(oma = c(6,6,1,1), mar = c(6,6,1,1))
	image(x = xAxAt, y = yAxAt, z = forPlot, col = colorRampPalette(c("red", "white", "blue"))(11), xlab = "", ylab = "", xaxt = "n", yaxt = "n", oldstyle = TRUE, zlim = c(-1,1))
	axis(1, at = xAxAt, labels = xAxLab, outer=FALSE, line=2, lwd=2, las=2)
	axis(2, at = yAxAt, labels = yAxLab, outer=FALSE, line=2, lwd=2, las=1)
	dev.off()
	# plot the color legend for it
	if (GLOBAL_VARIABLE_USE_SVG) {
		svg(file.path(rDir, paste(outfilePrefix, "_PCA_significant_correlation_legend.svg", sep = '')), height = 20, width = 20)
	} else {
		tiff(file.path(rDir, paste(outfilePrefix, "_PCA_significant_correlation_legend.tiff", sep = '')), width = 800, height = 800)
	}
	par(oma = c(0,0,0,0), mar = c(0,0,0,0))
	legend_image <- as.raster(matrix(colorRampPalette(c("red", "white", "blue"))(11), ncol=1))
	legend_image <- rev(legend_image)
	plot(c(0,2),c(-1,1),type = 'n', axes = FALSE, xlab = '', ylab = '', main = '')
	text(x=1.5, y = seq(-1,1,l=9), labels = seq(-1,1,l=9))
	rasterImage(legend_image, 0, -1, 1,1)
	dev.off()
	# write tables
	write.table(t(estimates), file.path(rDir, paste(outfilePrefix, "_PCA_significant_correlation_estimates.txt", sep = '')), sep = '\t', quote = FALSE)
	write.table(t(pValues), file.path(rDir, paste(outfilePrefix, "_PCA_significant_correlation_pValues.txt", sep = '')), sep = '\t', quote = FALSE)
	# collect all info for later
	out <- list(E = estimates, P = pValues)
	return(out)
}

#'@title Test for enrichment/depletion of genomic/epigenetic features in two different compartments.
#'@param pca a list holding tables with valid bins an their value of the first principle component (see \code{\link{f.internal.principle.component.analysis}}).
#'@param annotation a dataframe holding the annotation of the bins (see \code{\link{f.read.annotation}}).
#'@param featureNames names of the features which are drawn (used to simplify the annotation names).
#'@param rDir a directory where the figure is stored.
#'@param outfilePrefix a prefix for the figure and table files.
#'@param pValueThreshold specifies the significance-threshold (only significant values are drawn in the heatmaps).
#'@note internal function called by \code{\link{f.principle.component.analysis.and.features}}.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.pca.and.feature.wilcox.test <- function(pca, annotation, featureNames, rDir, outfilePrefix, pValueThreshold) {
	# pca: result from f.principle.component.analysis
	# annotation: the table holding fragment/bin annotations
	# featureNames: simplified feature names
	regNames <- names(pca)
	featuresToTest <- grep("^ann_|^sum_|^den_", colnames(annotation), value = TRUE)
	estimates <- matrix(0, nrow = length(regNames), ncol = length(featuresToTest), dimnames = list(regNames, featuresToTest))
	pValues <- matrix(0, nrow = length(regNames), ncol = length(featuresToTest), dimnames = list(regNames, featuresToTest))
	AtypeTables <- list()
	BtypeTables <- list()
	for (regName in regNames) {
		validBins <- pca[[regName]][,"validBins"]
		orientedPrincp <- pca[[regName]][,"orientedPrincp"]
		validAnnotation <- annotation[validBins,]
		Atype <- which(orientedPrincp > 0) # called it once active/open
		Btype <- which(orientedPrincp < 0) # called it once passive/close
		# balanced set sizes ?
		#if (length(Atype) > length(Btype)) {Atype <- sample(Atype, length(Btype))}
		#if (length(Atype) < length(Btype)) {Btype <- sample(Btype, length(Atype))}
		for (feature in featuresToTest){
			if ((substr(feature,1,4) == "sum_") || (substr(feature,1,4) == "ann_")) { estimates[regName, feature] <- mean(validAnnotation[Atype,feature])-mean(validAnnotation[Btype,feature]) }
			else { estimates[regName, feature] <- log2(mean(validAnnotation[Atype,feature])/mean(validAnnotation[Btype,feature])) }
			pValues[regName, feature] <- wilcox.test(validAnnotation[Atype,feature], validAnnotation[Btype,feature], alternative = "two.sided", exact = FALSE)$p.value
		}
		AtypeTables[[regName]] <- annotation[Atype,]
		BtypeTables[[regName]] <- annotation[Btype,]
	}
	# replace NAs and names
	estimates[is.na(estimates)] <- 0
	pValues[is.na(pValues)] <- 1
	colnames(estimates) <- featureNames
	colnames(pValues) <- featureNames
	# plot a heatmap with the significant fields
	enrichmentLimit <- ifelse(abs(max(estimates))>abs(min(estimates)), abs(max(estimates)), abs(min(estimates)))
	forPlot <- estimates
	forPlot[pValues>pValueThreshold] <- 0
	xAxAt <- 1:nrow(forPlot)
	yAxAt <- 1:ncol(forPlot)
	xAxLab <- rownames(forPlot)
	yAxLab <- colnames(forPlot)
	if (GLOBAL_VARIABLE_USE_SVG) {
		svg(file.path(rDir, paste(outfilePrefix, "_PCA_significant_enrichment.svg", sep = '')), height = 20, width = 20)
	} else {
		tiff(file.path(rDir, paste(outfilePrefix, "_PCA_significant_enrichment.tiff", sep = '')), width = 1600, height = 1600)
	}
	par(oma = c(6,6,1,1), mar = c(6,6,1,1))
	image(x = xAxAt, y = yAxAt, z = forPlot, col = colorRampPalette(c("red", "white", "blue"))(11), xlab = "", ylab = "", xaxt = "n", yaxt = "n", oldstyle = TRUE, zlim = c(-enrichmentLimit,enrichmentLimit))
	axis(1, at = xAxAt, labels = xAxLab, outer=FALSE, line=2, lwd=2, las=2)
	axis(2, at = yAxAt, labels = yAxLab, outer=FALSE, line=2, lwd=2, las=1)
	dev.off()
	# plot the color legend for it
	if (GLOBAL_VARIABLE_USE_SVG) {
		svg(file.path(rDir, paste(outfilePrefix, "_PCA_significant_enrichment_legend.svg", sep = '')), height = 20, width = 20)
	} else {
		tiff(file.path(rDir, paste(outfilePrefix, "_PCA_significant_enrichment_legend.tiff", sep = '')), width = 800, height = 800)
	}
	par(oma = c(0,0,0,0), mar = c(0,0,0,0))
	legend_image <- as.raster(matrix(colorRampPalette(c("red", "white", "blue"))(11), ncol=1))
	legend_image <- rev(legend_image)
	plot(c(0,2),c(-1,1),type = 'n', axes = FALSE, xlab = '', ylab = '', main = '')
	text(x=1.5, y = seq(-1,1,l=9), labels = seq(-enrichmentLimit,enrichmentLimit,l=9))
	rasterImage(legend_image, 0, -1, 1,1)
	dev.off()
	# write tables
	write.table(t(estimates), file.path(rDir, paste(outfilePrefix, "_PCA_significant_enrichment_estimates.txt", sep = '')), sep = '\t', quote = FALSE)
	write.table(t(pValues), file.path(rDir, paste(outfilePrefix, "_PCA_significant_enrichment_pValues.txt", sep = '')), sep = '\t', quote = FALSE)
	# collect all info for later
	out <- list(E = estimates, P = pValues, AtypeTables = AtypeTables, BtypeTables = BtypeTables)
	return(out)
}

#'@title Plot enrichment/depletion of genomic/epigenetic features in two different compartments (boxplots).
#'@param AtypeTables observed values of annotation features in the first type of the compartments.
#'@param BtypeTables observed values of annotation features in the second type of the compartments.
#'@param featureNames names of the features which are drawn (used to simplify the annotation names).
#'@param rDir a directory where the figure is stored.
#'@param outfilePrefix a prefix for the figure file.
#'@note internal function called by \code{\link{f.principle.component.analysis.and.features}}.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.A.and.B.type.boxplot <- function(AtypeTables, BtypeTables, featureNames, rDir, outfilePrefix) {
	# AtypeTables: result from f.internal.pca.and.feature.wilcox.test >> out$AtypeTables
	# BtypeTables: result from f.internal.pca.and.feature.wilcox.test >>
	# featureNames: simplified feature names
	Atype <- do.call("rbind", AtypeTables)
	Btype <- do.call("rbind", BtypeTables)
	genericNames <- grep("^ann_|^sum_|^den_", colnames(Atype), value = TRUE)
	Atype <- Atype[,genericNames] # remove num, chrom, start, end
	Btype <- Btype[,genericNames] # remove num, chrom, start, end
	colnames(Atype) <- featureNames
	colnames(Btype) <- featureNames
	numA <- nrow(Atype)
	numB <- nrow(Btype)
	annCols <- featureNames[grep("^ann_", genericNames)]
	sumCols <- featureNames[grep("^sum_", genericNames)]
	denCols <- featureNames[grep("^den_", genericNames)]
	numAnnCols <- length(annCols)
	numSumCols <- length(sumCols)
	numDenCols <- length(denCols)
	numTotCols <- numSumCols + numAnnCols + numDenCols
	differentCols <- list()
	if (numAnnCols > 0) { differentCols$annCols <- annCols }
	if (numSumCols > 0) { differentCols$sumCols <- sumCols }
	if (numDenCols > 0) { differentCols$denCols <- denCols }
	numDifferentCols <- length(differentCols)
	# plotting
	if (GLOBAL_VARIABLE_USE_SVG) {
		svg(file.path(rDir, paste(outfilePrefix, "_PCA_A_and_B_type_boxplots.svg", sep = '')), height = numTotCols+3, width = 10)
	} else {
		tiff(file.path(rDir, paste(outfilePrefix, "_PCA_A_and_B_type_boxplots.tiff", sep = '')), height = 100*numTotCols+300, width = 1000)
	}
	layout(mat = matrix(c(1:numDifferentCols), ncol = 1), widths = rep(1, numDifferentCols), heights = c(numAnnCols+1, numSumCols+1, numDenCols+1))
	par(oma = c(2,10,0,0), mar = c(2,10,0,0))
	for (colType in names(differentCols)) {
		features <- differentCols[[colType]]
		numFeatures <- length(features)
		if (numFeatures > 1) {
			AtypeSummary <- apply(Atype[,features], 2, function(x) boxplot.stats(x))
			BtypeSummary <- apply(Btype[,features], 2, function(x) boxplot.stats(x))
		} else {
			AtypeSummary <- list(boxplot.stats(Atype[,features]))
			BtypeSummary <- list(boxplot.stats(Btype[,features]))
			names(AtypeSummary) <- features
			names(BtypeSummary) <- features
		}
		listForPlot <- list(stats = matrix(0, nrow = 5, ncol = 2*numFeatures),
							n = rep(c(numA, numB), numFeatures),
							conf = matrix(0, nrow = 2, ncol = 2*numFeatures),
							out = c(),# a vector with outlier y values. The x position is defined by the next vector "group"
							group = c(),
							names = paste(rep(features, each = 2), c("_A", "_B"), sep = '')
		)
		groupNumber <- 0
		for (fn in names(AtypeSummary)) {
			groupNumber <- groupNumber + 1
			listForPlot$stats[,groupNumber] <- AtypeSummary[[fn]]$stats
			listForPlot$conf[,groupNumber] <- AtypeSummary[[fn]]$conf
			listForPlot$out <- c(listForPlot$out, AtypeSummary[[fn]]$out)
			listForPlot$group <- c(listForPlot$group, rep(groupNumber, length(AtypeSummary[[fn]]$out)))
			groupNumber <- groupNumber + 1
			listForPlot$stats[,groupNumber] <- BtypeSummary[[fn]]$stats
			listForPlot$conf[,groupNumber] <- BtypeSummary[[fn]]$conf
			listForPlot$out <- c(listForPlot$out, BtypeSummary[[fn]]$out)
			listForPlot$group <- c(listForPlot$group, rep(groupNumber, length(BtypeSummary[[fn]]$out)))
		}
		bxp(listForPlot, xlab = "", ylab = "", main = "", bty = "n", pch = 19, cex = 0.5, horizontal = TRUE, las = 1, boxfill = rep(c("blue", "red"), numFeatures))
	}
	dev.off()
}

#'@title Analysis of the first principle component (PCA).
#'@param dataMatrix a HiC interaction matrix (see \code{\link{f.load.one.sample}}).
#'@param binSize size of genomic bins in bp. 0 (restriction fragments) does not work.
#'@param rDir a directory where the figure is stored.
#'@param outfilePrefix a prefix for the figure and table files.
#'@param annotation a dataframe holding the annotation of the bins (see \code{\link{f.read.annotation}}).
#'@param regionTable defines specific regions in the genome which shall be assessed. Per default,
#'each chromosome is tested separately. The table must have three columns (chrom, start, and end). User defined names
#'can be given using an optional fourth column (name). An example is given in the tutorial \url{https://github.com/MWSchmid/HiCdat}.
#'@param simplifiedNames is a list with the column names of the annotation as keys to simplified names as values
#'(e.g. ann_transposable_element_gene can be replaced by TE-gene).
#'@param filterZero whether or not to filter the x percent of bins with the highest number of 0 entries.
#'@param filterThreshold if filterZero is given, this specifies the fraction of bins which shall be kept.
#'@param pValueThreshold specifies the significance-threshold for the correlation and enrichment tests
#'(only significant values are drawn in the heatmaps).
#'@param userLimits define your own min/max for the correlation values (related to the plot and its color mapping).
#'@return a list with the correlation estimates and P-values (estimate, estimatePvalue), enrichment estimates and P-values (enrichment, enrichmentPvalue),
#'a list with all the first principle components for all the regions defined in \code{regionTable} (regToPC), and the values used for the enrichment test (AtypeTables, BtypeTables).
#'@examples \dontrun{
#'f.principle.component.analysis.and.features(
#'  dataMatrix = dataMatrixSampleX,
#'  binSize = 1e5,
#'  rDir = "/path/to/where/the/results/are/stored",
#'  outfilePrefix = "aPrefixForTheFileNames",
#'  annotation = annotationTable, # or data.frame() if only PCA is requested
#'  regionTable = data.frame(),
#'  simplifiedNames = list(),
#'  filterZero = TRUE,
#'  filterThreshold = 0.95,
#'  pValueThreshold = 0.05,
#'  userLimits = c(-1, 1)
#') }
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@export
f.principle.component.analysis.and.features <- function(dataMatrix, binSize, rDir, outfilePrefix, annotation = data.frame(), regionTable = data.frame(), simplifiedNames = list(), filterZero = TRUE, filterThreshold = 0.95, pValueThreshold = 0.05, userLimits = c(-1, 1)) {
	# check if there is an annotation, otherwise just do the PCA and return
	if (nrow(annotation) == 0) { return(f.internal.principle.component.analysis(dataMatrix, binSize, rDir, paste(outfilePrefix, "_first_component", sep = ''), regionTable, filterZero, filterThreshold, userLimits)) }
	# check the annotation table and set the names of the features
	featuresToTest <- grep("^ann_|^sum_|^den_", colnames(annotation), value = TRUE)
	featureNames <- featuresToTest
	if (length(featuresToTest) != (ncol(annotation)-4)) { cat("note that there are one or more invalid columns (without ann_, sum_, den_ in the name)\n"); return(FALSE); }
	if (length(simplifiedNames) > 0) { withinList <- featureNames %in% names(simplifiedNames); featureNames[withinList] <- unlist(simplifiedNames)[featureNames[withinList]]; } 
	# do the PCA
	cat("PCA...\n")
	pca <- f.internal.principle.component.analysis(dataMatrix, binSize, rDir, paste(outfilePrefix, "_first_component", sep = ''), regionTable, filterZero, filterThreshold, userLimits)
	# do the tests with the genomic and epigenetic features
	cat("correlation test...\n")
	corTest <- f.internal.pca.and.feature.correlation.test(pca, annotation, featureNames, rDir, outfilePrefix, pValueThreshold)
	cat("enrichment test...\n")
	coxTest <- f.internal.pca.and.feature.wilcox.test(pca, annotation, featureNames, rDir, outfilePrefix, pValueThreshold)
	# draw a boxplot for the A and B type compartments
	cat("boxplot drawing...\n")
	f.internal.A.and.B.type.boxplot(coxTest$AtypeTables, coxTest$BtypeTables, featureNames, rDir, outfilePrefix)
	# return the matrices
	out <- list( estimate = corTest$E, estimatePvalue = corTest$P, enrichment = coxTest$E, enrichmentPvalue = coxTest$P, regToPC = pca, AtypeTables = coxTest$AtypeTables, BtypeTables = coxTest$BtypeTables )
	return(out)
}

#'@title Calculate the sum of interaction frequencies of a certain set of bins.
#'@param dataMatrix a HiC interaction matrix (see \code{\link{f.load.one.sample}}).
#'@param seTab a table with the start and end bins of the regions of interest.
#'@note internal function called by \code{\link{f.test.interaction.frequencies}}.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.calculate.interaction.sums <- function(dataMatrix, seTab) {
	# seTab: a table with two columns: start and end bins of a certain region (bin numbers)
	out <- 0
	for (i in 1:(nrow(seTab)-1)) {
		curBins <- seTab[i, 1]:seTab[i, 2]
		for (j in (i+1):nrow(seTab)) {
			otherBins <- seTab[j, 1]:seTab[j, 2]
			out <- out + sum(dataMatrix[curBins, otherBins])
		}
	}
	return(out)
}

#'@title Test a set of regions for increased interaction frequency among each other compared to random sets of regions.
#'@param dataMatrix a HiC interaction matrix (see \code{\link{f.load.one.sample}}).
#'@param binSize size of genomic bins in bp. 0 (restriction fragments) does not work.
#'@param repetitions number of random sets to used for the empirical distribution.
#'@param testRegionsTable is a table with genomic regions of interest. Columns must be chrom, start, end. Rownames can be freely chosen.
#'@param regionDefinitionTable defines specific regions in the genome from which the random sets are sampled. If a table is supplied,
#'the regions of interest are assigned to the defined regions (unassigned will be removed) and sampling happens within the defined region.
#'If no regions are defined, the regions of interest are assigned to the whole chromosomes. Mapping of regions is done on the whole length. 
#'Hence, a test region must be entirely within a defined region to be assigned to it - unassigned test regions are removed.
#'Note that the defined regions should not overlap (test regions can do so). The table must have three columns (chrom, start, and end).
#'User defined names can be given using an optional fourth column (name). An example is given in the A. thaliana tutorial.
#'@return a vector holding the sum of interactions between the regions of interest, the mean and standard deviation of the sampled interaction sums, and the corresponding P-value
#'@examples \dontrun{
#'f.test.interaction.frequencies(
#'  dataMatrix = dataMatrixSampleX,
#'  binSize = 1e5,
#'  repetitions = 1e4,
#'  testRegionsTable = tableWithRegionsOfInterest,
#'  regionDefinitionTable = data.frame()
#') }
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@export
f.test.interaction.frequencies <- function(dataMatrix, binSize, repetitions, testRegionsTable, regionDefinitionTable = data.frame()) {
	cat("preparing...\n")
	# set up the region definition (create a default if no user specified region exists, add names if not given)
	if (nrow(regionDefinitionTable) == 0) { chromSizes <- f.internal.get.chrom.sizes(); regionDefinitionTable <- data.frame(chrom = names(chromSizes), start = rep(0, length(chromSizes)), end = unlist(chromSizes), stringsAsFactors = FALSE) }
	if (ncol(regionDefinitionTable) == 3) { regionDefinitionTable$name <- paste(regionDefinitionTable$chrom, regionDefinitionTable$start, regionDefinitionTable$end, sep = '_') }
	rownames(regionDefinitionTable) <- regionDefinitionTable$name
	# map the test regions to the defined regions
	testRegionsTable$assignedTo <- "empty"
	for (i in 1:nrow(testRegionsTable)) {
		temp <- (regionDefinitionTable$chrom == testRegionsTable$chrom[i]) & (regionDefinitionTable$start <= testRegionsTable$start[i]) & (regionDefinitionTable$end >= testRegionsTable$end[i])
		testRegionsTable$assignedTo[i] <- ifelse(sum(temp) == 1, regionDefinitionTable$name[temp], "none")
	}
	numRemoved <- sum(testRegionsTable$assignedTo == "none")
	cat(paste("removed ", numRemoved, " unassigned regions: ", paste(rownames(testRegionsTable)[testRegionsTable$assignedTo == "none"], collapse = ','), '\n', sep = '')) #cat(paste("removed ", numRemoved, " unassigned regions\n", sep = ''))
	testRegionsTable <- testRegionsTable[testRegionsTable$assignedTo!="none",]#subset(testRegionsTable, assignedTo != "none")
	# get a seList
	seList <- f.get.se.list(binSize)
	# get the interaction sums of the test regions
	testBinStart <- f.translate.chrom.pos.vector.to.index(testRegionsTable$chrom, testRegionsTable$start, seList, binSize)
	testBinEnd <- f.translate.chrom.pos.vector.to.index(testRegionsTable$chrom, testRegionsTable$end, seList, binSize)
	testSum <- f.internal.calculate.interaction.sums(dataMatrix, cbind(testBinStart, testBinEnd))
	# sample for each test regions one table with random binStarts and binEnds
	if (nrow(testRegionsTable) > 50) { cat("sampling (# == 1 test region):\n") }
	else { cat("sampling (# == 1 test region):") }
	sampledTabs <- list()
	for (i in 1:nrow(testRegionsTable)) {
		if ((i %% 10) == 0) { cat(" ") }
		if ((i %% 50) == 0) { cat("\n") }
		cat("#")
		testReg <- testRegionsTable$assignedTo[i]
		testSize <- testRegionsTable$end[i] - testRegionsTable$start[i]
		defChrom <- regionDefinitionTable[testReg, "chrom"]
		defStart <- regionDefinitionTable[testReg, "start"]
		defEnd <- regionDefinitionTable[testReg, "end"]
		numPossiblePos <- (defEnd-defStart)-testSize
		sampledPositions <- sample.int(numPossiblePos, repetitions)+defStart+floor(testSize/2)
		sampledBinStart <- f.translate.chrom.pos.vector.to.index(defChrom, sampledPositions-floor(testSize/2), seList, binSize)
		sampledBinEnd <- f.translate.chrom.pos.vector.to.index(defChrom, sampledPositions+floor(testSize/2), seList, binSize)
		sampledTabs[[i]] <- cbind(sampledBinStart, sampledBinEnd)
	}
	sampledTabs <- do.call("rbind", sampledTabs)
	cat("\n")
	# do the tests
	cat("calculating P-values (# == 1000 tests): ")
	# does not work like this: sampledSums <- aggregate(sampledTabs, by = list(repNum = rep(1:repetitions, nrow(testRegionsTable))), function(x) f.internal.calculate.interaction.sums(x, dataMatrix))
	sampledSums <- rep(0, repetitions)
	startLines <- repetitions*(0:(nrow(testRegionsTable)-1))
	for (i in 1:repetitions) {
		curLines <- startLines+i
		sampledSums[i] <- f.internal.calculate.interaction.sums(dataMatrix, sampledTabs[curLines,])
		if ((i %% 1000) == 0) { cat("#") }
	}
	cat("\n")
	#out <- list(testSum = log2(testSum+1), sampledMean = mean(log2(sampledSums+1)), sampledSd = sd(log2(sampledSums+1)), pValue = sum(sampledSums>testSum)/repetitions, sampledSums = sampledSums)
	out <- c(testSum = log2(testSum+1), sampledMean = mean(log2(sampledSums+1)), sampledSd = sd(log2(sampledSums+1)), pValue = sum(sampledSums>testSum)/repetitions)
	return(out)
}

#'@title Summarize the annotation along bins using bin indices.
#'@param annotation a dataframe holding the annotation of the bins (see \code{\link{f.read.annotation}}).
#'@param startBin first bin of the interval to be summarized.
#'@param endBin last bin of the interval to be summarized.
#'@param dataColsForSum names of the columns which are summarized using sum().
#'@param dataColsForMean names of the columns which are summarized using mean().
#'@note internal function called by \code{\link{f.test.regions.for.feature.enrichment.fragment.based}} or \code{\link{f.test.regions.for.feature.enrichment.bin.based}}.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.summarize.binned.annotation <- function(annotation, startBin, endBin, dataColsForSum, dataColsForMean) {
	subDat <- annotation[startBin:endBin,]
	if (is.vector(subDat)) {
		if (length(dataColsForSum) >= 1) { sumResults <- subDat[dataColsForSum] }
		else { sumResults <- NA }
		if (length(dataColsForMean) >= 1) { meanResults <- subDat[dataColsForMean] }
		else { meanResults <- NA }
	} else {
		if (length(dataColsForSum) > 1) { sumResults <- apply(subDat[,dataColsForSum], 2, sum) }
		else if (length(dataColsForSum) == 1) { sumResults <- sum(subDat[,dataColsForSum]) }
		else { sumResults <- NA }
		if (length(dataColsForMean) > 1) { meanResults <- apply(subDat[,dataColsForMean], 2, mean) }
		else if (length(dataColsForMean) == 1) { meanResults <- mean(subDat[,dataColsForMean]) }
		else { meanResults <- NA }
	}
	out <- c(sumResults, meanResults)
	out <- out[!is.na(out)]
	return(out)
}

#'@title Summarize the annotation along fragments/bins using genomic coordinates.
#'@param annotation a dataframe holding the annotation of the bins (see \code{\link{f.read.annotation}}).
#'@param c chromosome of the region to be summarized.
#'@param s start position of the region to be summarized.
#'@param e end position of the region to be summarized.
#'@param dataColsForSum names of the columns which are summarized using sum().
#'@param dataColsForMean names of the columns which are summarized using mean().
#'@note internal function called by \code{\link{f.test.regions.for.feature.enrichment.fragment.based}} or \code{\link{f.test.regions.for.feature.enrichment.bin.based}}.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.summarize.annotation <- function(annotation, c, s, e, dataColsForSum, dataColsForMean) {
	subDat <- annotation[((annotation$chrom == c) & (annotation$start >= s) & (annotation$end <= e)),] #subset(annotation, (chrom == c) & (start >= s) & (end <= e))
	if (length(dataColsForSum) > 1) { sumResults <- apply(subDat[,dataColsForSum], 2, sum) }
	else if (length(dataColsForSum) == 1) { sumResults <- sum(subDat[,dataColsForSum]) }
	else { sumResults <- NA }
	if (length(dataColsForMean) > 1) { meanResults <- apply(subDat[,dataColsForMean], 2, mean) }
	else if (length(dataColsForMean) == 1) { meanResults <- mean(subDat[,dataColsForMean]) }
	else { meanResults <- NA }
	out <- c(sumResults, meanResults)
	out <- out[!is.na(out)]
	return(out)
}

#'@title Summarize the annotation along fragments/bins using genomic coordinates without checking for the chromosome.
#'@param annotation a dataframe holding the annotation of the bins (see \code{\link{f.read.annotation}}).
#'@param s start position of the region to be summarized.
#'@param e end position of the region to be summarized.
#'@param dataColsForSum names of the columns which are summarized using sum().
#'@param dataColsForMean names of the columns which are summarized using mean().
#'@note internal function called by \code{\link{f.test.regions.for.feature.enrichment.fragment.based}} or \code{\link{f.test.regions.for.feature.enrichment.bin.based}}.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.internal.summarize.annotation}}
f.internal.summarize.annotation.no.chrom.check <- function(annotation, s, e, dataColsForSum, dataColsForMean) {
	subDat <- annotation[((annotation$start >= s) & (annotation$end <= e)),] #subset(annotation, (start >= s) & (end <= e))
	subDat <- annotation[((annotation$chrom == c) & (annotation$start >= s) & (annotation$end <= e)),] #subset(annotation, (chrom == c) & (start >= s) & (end <= e))
	if (length(dataColsForSum) > 1) { sumResults <- apply(subDat[,dataColsForSum], 2, sum) }
	else if (length(dataColsForSum) == 1) { sumResults <- sum(subDat[,dataColsForSum]) }
	else { sumResults <- NA }
	if (length(dataColsForMean) > 1) { meanResults <- apply(subDat[,dataColsForMean], 2, mean) }
	else if (length(dataColsForMean) == 1) { meanResults <- mean(subDat[,dataColsForMean]) }
	else { meanResults <- NA }
	out <- c(sumResults, meanResults)
	out <- out[!is.na(out)]
	return(out)
}

#'@title Calculate the enrichment given observed and sampled values.
#'@param testValues observed values.
#'@param sampledValues values obtained from random sampling.
#'@param countDataWasLogged specifies if the count data was already log-transformed. Generally TRUE if \code{useLog = TRUE} while loading the annotation with \code{\link{f.read.annotation}}.
#'@note internal function called by \code{\link{f.test.regions.for.feature.enrichment.fragment.based}} or \code{\link{f.test.regions.for.feature.enrichment.bin.based}}.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.calculate.region.enrichment <- function(testValues, sampledValues, countDataWasLogged = TRUE) {
	# NOTE that sampledValues must have ncol == number of ann/sum/den features and nrow == number of repetitions
	# the naming is taken from f.process.one.kee, the enrichment is calculated differently for the two groups
	# summed up (ann and sum) features will be log2(observed/expected)
	# meaned (den) features will be observed/expected
	# note - log for all is better for the heatmap (+ and -)
	obsToExp <- log2(testValues/apply(sampledValues, 2, mean))
	if (countDataWasLogged) {
		dataColsForSum <- grep("^sum_|^ann_", colnames(sampledValues), value = TRUE)
		if (length(dataColsForSum) > 1) { obsToExp[dataColsForSum] <- testValues[dataColsForSum] - apply(sampledValues[,dataColsForSum], 2, mean) }
		if (length(dataColsForSum) == 1) { obsToExp[dataColsForSum] <- testValues[dataColsForSum] - mean(sampledValues[,dataColsForSum]) }
	}
	#obsToExp <- log2(testValues/apply(sampledValues, 2, mean))
	#toLog <- grep("^sum_|^ann_", names(obsToExp), value = TRUE)
	#obsToExp[toLog] <- log2(obsToExp[toLog])
	return(obsToExp)
}

#'@title Plot significant enrichment of genomic/epigenetic features in a heatmap.
#'@param observed observed values.
#'@param enrichment enrichment compared to random sampling (see also \code{\link{f.internal.calculate.region.enrichment}}).
#'@param pValues the P-values specifying the significance of the enrichment.
#'@param pValueThreshold specifies the significance-threshold (only significant values are drawn in the heatmaps).
#'@param rDir a directory where the figure is stored.
#'@param outfilePrefix a prefix for the figure and table files.
#'@param featureNames names of the features which are drawn (used to simplify the annotation names).
#'@note internal function called by \code{\link{f.test.regions.for.feature.enrichment.fragment.based}} or \code{\link{f.test.regions.for.feature.enrichment.bin.based}}.
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
f.internal.plot.and.save.test.regions.for.feature.enrichment.heatmap <- function(observed, enrichment, pValues, pValueThreshold, rDir, outfilePrefix, featureNames) {
	# observed, enrichment, pValues are results from f.test.regions.for.feature.enrichment.fragment.based
	# featureNames: simplified feature names
	# for shorter writing:
	colnames(pValues) <- featureNames
	colnames(observed) <- featureNames
	colnames(enrichment) <- featureNames
	write.table(pValues, file.path(rDir, paste(outfilePrefix, "_pValues.txt", sep ='')), sep = '\t', quote = FALSE)
	write.table(observed, file.path(rDir, paste(outfilePrefix, "_observed.txt", sep ='')), sep = '\t', quote = FALSE)
	write.table(enrichment, file.path(rDir, paste(outfilePrefix, "_enrichment.txt", sep ='')), sep = '\t', quote = FALSE)
	## plot a heatmap with the significant fields
	forPlot <- enrichment
	forPlot[(pValues>pValueThreshold) & (pValues<(1-pValueThreshold))] <- 0
	xAxAt <- 1:nrow(forPlot)
	yAxAt <- 1:ncol(forPlot)
	xAxLab <- rownames(forPlot)
	yAxLab <- colnames(forPlot)
	maxValToPlot <- ceiling(max(abs(forPlot)))
	if (maxValToPlot > 10) {
		zLow <- -10
		zHigh <- 10
		forPlot[forPlot < zLow] <- zLow
		forPlot[forPlot > zHigh] <- zHigh
	} else {
		zLow <- -maxValToPlot
		zHigh <- maxValToPlot
	}
	if (GLOBAL_VARIABLE_USE_SVG) {
		svg(file.path(rDir, paste(outfilePrefix, "_enrichment.svg", sep = '')), height = 20, width = 20)
	} else {
		tiff(file.path(rDir, paste(outfilePrefix, "_enrichment.tiff", sep = '')), height = 1600, width = 1600)
	}
	par(oma = c(6,6,1,1), mar = c(6,6,1,1))
	image(x = xAxAt, y = yAxAt, z = forPlot, col = colorRampPalette(c("red", "white", "blue"))(11), xlab = "", ylab = "", xaxt = "n", yaxt = "n", oldstyle = TRUE, zlim = c(zLow,zHigh))
	axis(1, at = xAxAt, labels = xAxLab, outer=FALSE, line=2, lwd=2, las=2)
	axis(2, at = yAxAt, labels = yAxLab, outer=FALSE, line=2, lwd=2, las=1)
	dev.off()
	## and a color legend for it
	if (GLOBAL_VARIABLE_USE_SVG) {
		svg(file.path(rDir, paste(outfilePrefix, "_enrichment_legend.svg", sep = '')), height = 20, width = 20)
	} else {
		tiff(file.path(rDir, paste(outfilePrefix, "_enrichment_legend.tiff", sep = '')), height = 800, width = 800)
	}
	par(oma = c(0,0,0,0), mar = c(0,0,0,0))
	legend_image <- as.raster(matrix(colorRampPalette(c("red", "white", "blue"))(11), ncol=1))
	legend_image <- rev(legend_image)
	plot(c(0,2),c(-1,1),type = 'n', axes = FALSE, xlab = '', ylab = '', main = '')
	text(x=1.5, y = seq(-1,1,l=2*zHigh+1), labels = seq(zLow,zHigh,l=2*zHigh+1))
	rasterImage(legend_image, 0, -1, 1, 1)
	dev.off()
	out <- list(observed = observed, enrichment = enrichment, pValues = pValues)
	return(out)
}

#'@title Test if a set of regions has a specific genomic/epigenetic makeup compared to random sets of regions.
#'@param annotation a dataframe holding the annotation of the bins (see \code{\link{f.read.annotation}}).
#'@param rDir a directory where the figure is stored.
#'@param outfilePrefix a prefix for the figure and table files.
#'@param repetitions number of random sets to used for the empirical distribution.
#'@param testRegionsTable is a table with genomic regions of interest. Columns must be chrom, start, end. Rownames can be freely chosen.
#'@param regionDefinitionTable defines specific regions in the genome from which the random sets are sampled. If a table is supplied,
#'the regions of interest are assigned to the defined regions (unassigned will be removed) and sampling happens within the defined region.
#'If no regions are defined, the regions of interest are assigned to the whole chromosomes. Mapping of regions is done on the whole length. 
#'Hence, a test region must be entirely within a defined region to be assigned to it - unassigned test regions are removed.
#'Note that the defined regions should not overlap (test regions can do so). The table must have three columns (chrom, start, and end).
#'User defined names can be given using an optional fourth column (name). An example is given in the A. thaliana tutorial.
#'@param simplifiedNames is a list with the column names of the annotation as keys to simplified names as values
#'(e.g. ann_transposable_element_gene can be replaced by TE-gene).
#'@param pValueThreshold specifies the significance-threshold (only significant values are drawn in the heatmaps).
#'@param countDataWasLogged specifies if the count data was already log-transformed. Generally TRUE if \code{useLog = TRUE} while loading the annotation with \code{\link{f.read.annotation}}.
#'@return a list holding the observed values, the enrichment compared to the random sets, and the corresponding P-values.
#'@examples \dontrun{
#'enrichmentResult <- f.test.regions.for.feature.enrichment.fragment.based(
#'  annotation = annotationTableOnRestrictionFragments,
#'  rDir = "/path/to/where/the/results/are/stored",
#'  outfilePrefix = "aPrefixForTheFileNames",
#'  repetitions = 1e4,
#'  testRegionsTable = tableWithRegionsOfInterest,
#'  regionDefinitionTable = data.frame(),
#'  simplifiedNames = list(),
#'  pValueThreshold = 0.05,
#'  countDataWasLogged = TRUE
#') }
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.test.regions.for.feature.enrichment.bin.based}}
#'@export
f.test.regions.for.feature.enrichment.fragment.based <- function(annotation, rDir, outfilePrefix, repetitions, testRegionsTable, regionDefinitionTable = data.frame(), simplifiedNames = list(), pValueThreshold = 0.05, countDataWasLogged = TRUE) {
	cat("preparing...\n")
	# check the annotation table and set the names of the features
	featuresToTest <- grep("^ann_|^sum_|^den_", colnames(annotation), value = TRUE)
	featureNames <- featuresToTest
	if (length(featuresToTest) != (ncol(annotation)-4)) { cat("note that there are one or more invalid columns (without ann_, sum_, den_ in the name)\n"); return(FALSE); }
	if (length(featuresToTest) == 1) { 
		cat("testing only one feature - duplicating it for simplicity.\n")
		dummyName <- paste(featuresToTest, "duplicate", sep = '_')
		annotation[[dummyName]] <- annotation[[featuresToTest]]
		featuresToTest <- c(featuresToTest, dummyName)
		featureNames <- featuresToTest
	}
	if (length(simplifiedNames) > 0) { withinList <- featureNames %in% names(simplifiedNames); featureNames[withinList] <- unlist(simplifiedNames)[featureNames[withinList]]; } 
	# set up the region definition (create a default if no user specified region exists, add names if not given)
	if (nrow(regionDefinitionTable) == 0) { chromSizes <- f.internal.get.chrom.sizes(); regionDefinitionTable <- data.frame(chrom = names(chromSizes), start = rep(0, length(chromSizes)), end = unlist(chromSizes), stringsAsFactors = FALSE) }
	if (ncol(regionDefinitionTable) == 3) { regionDefinitionTable$name <- paste(regionDefinitionTable$chrom, regionDefinitionTable$start, regionDefinitionTable$end, sep = '_') }
	rownames(regionDefinitionTable) <- regionDefinitionTable$name
	# map the test regions to the defined regions
	testRegionsTable$assignedTo <- "empty"
	for (i in 1:nrow(testRegionsTable)) {
		temp <- (regionDefinitionTable$chrom == testRegionsTable$chrom[i]) & (regionDefinitionTable$start <= testRegionsTable$start[i]) & (regionDefinitionTable$end >= testRegionsTable$end[i])
		testRegionsTable$assignedTo[i] <- ifelse(sum(temp) == 1, regionDefinitionTable$name[temp], "none")
	}
	numRemoved <- sum(testRegionsTable$assignedTo == "none")
	cat(paste("removed ", numRemoved, " unassigned regions: ", paste(rownames(testRegionsTable)[testRegionsTable$assignedTo == "none"], collapse = ','), '\n', sep = '')) #cat(paste("removed ", numRemoved, " unassigned regions\n", sep = ''))
	testRegionsTable <- testRegionsTable[testRegionsTable$assignedTo!="none",]#subset(testRegionsTable, assignedTo != "none")
	# for each region, sample <repetitions> random regions, calculate individual pValues as well
	cat("sampling and testing (# == 1 test region): ")
	sampledTabs <- list()
	observedTabs <- list()
	individualPvalueTabs <- list()
	enrichmentTabs <- list()
	for (i in 1:nrow(testRegionsTable)) {
		cat("#")
		testReg <- testRegionsTable$assignedTo[i]
		testSize <- testRegionsTable$end[i] - testRegionsTable$start[i]
		defChrom <- regionDefinitionTable[testReg, "chrom"]
		defStart <- regionDefinitionTable[testReg, "start"]
		defEnd <- regionDefinitionTable[testReg, "end"]
		numPossiblePos <- (defEnd-defStart)-testSize
		subAnno <- annotation[((annotation$chrom == defChrom) & (annotation$start >= defStart) & (annotation$end <= defEnd)),] #subset(annotation, (chrom == defChrom) & (start >= defStart) & (end <= defEnd))
		# depending on the annotation data one would like to have sums or means
		dataColsForSum <- grep("^sum_|^ann_", colnames(subAnno), value = TRUE)
		dataColsForMean <- grep("^den_", colnames(subAnno), value = TRUE)
		resultsHeader <- c(dataColsForSum, dataColsForMean)
		# calculate the real values
		real <- f.internal.summarize.annotation.no.chrom.check(subAnno, testRegionsTable$start[i], testRegionsTable$end[i], dataColsForSum, dataColsForMean)
		names(real) <- resultsHeader
		observedTabs[[i]] <- real
		# sample random positions and calculate their values this will give a matrix with repetitions rows and length(resultsHeader) columns
		# TO AVOID NA - OVERSAMPLE BY 20 %, remove NA and then sample within the working set
		toAdd <- floor(testSize/2)
		sampledPositions <- sample.int(numPossiblePos, repetitions+ceiling(1.2*repetitions))+defStart+toAdd
		sampled <- t(sapply(sampledPositions, function(x) f.internal.summarize.annotation.no.chrom.check(subAnno, x-toAdd, x+toAdd, dataColsForSum, dataColsForMean)))
		sampled <- apply(sampled, 2, function(x) sample(x[!is.na(x)], repetitions))
		sampledTabs[[i]] <- sampled
		# calculate individual P value
		individualPvalueTabs[[i]] <- apply(t(sampled) > real, 1, sum)/repetitions # this here is working - just be aware with the t() and this stuff - below is the safer version
		# enrichment
		enrichmentTabs[[i]] <- f.internal.calculate.region.enrichment(real, sampled, countDataWasLogged)
	}
	cat("\n")
	cat("combining results...\n")
	sampledTabs <- do.call("rbind", sampledTabs)
	observedTabs <- do.call("rbind", observedTabs)
	individualPvalueTabs <- do.call("rbind", individualPvalueTabs)
	enrichmentTabs <- do.call("rbind", enrichmentTabs)
	rownames(observedTabs) <- rownames(testRegionsTable)
	rownames(individualPvalueTabs) <- rownames(testRegionsTable)
	rownames(enrichmentTabs) <- rownames(testRegionsTable)
	colsToKeep <- colnames(sampledTabs)
	sampledMean <- aggregate(sampledTabs, by = list(sn = rep(c(1:repetitions), nrow(testRegionsTable))), mean)
	sampledMean <- sampledMean[,colsToKeep]
	combinedMean <- apply(observedTabs, 2, mean)
	combinedPvalue <- apply(t(sampledMean) > combinedMean, 1, sum)/repetitions
	combinedEnrichment <- f.internal.calculate.region.enrichment(combinedMean, sampledMean, countDataWasLogged)
	observedTabs <- rbind(observedTabs, combinedMean)
	pValueTabs <- rbind(individualPvalueTabs, combinedPvalue)
	enrichmentTabs <- rbind(enrichmentTabs, combinedEnrichment)
	cat("saving results...\n")
	# note that out is like list(observed = observedTabs, enrichment = enrichmentTabs, pValues = pValueTabs) but with the generic feature names replaced by simplifiedNames
	out <- f.internal.plot.and.save.test.regions.for.feature.enrichment.heatmap(observedTabs, enrichmentTabs, pValueTabs, pValueThreshold, rDir, outfilePrefix, featureNames)
	return(out)
}

#'@title Test if a set of regions has a specific genomic/epigenetic makeup compared to random sets of regions.
#'@param annotation a dataframe holding the annotation of the bins (see \code{\link{f.read.annotation}}).
#'@param binSize size of genomic bins in bp. 0 (restriction fragments) does not work, use \code{\link{f.test.regions.for.feature.enrichment.fragment.based}} instead.
#'@param rDir a directory where the figure is stored.
#'@param outfilePrefix a prefix for the figure and table files.
#'@param repetitions number of random sets to used for the empirical distribution.
#'@param testRegionsTable is a table with genomic regions of interest. Columns must be chrom, start, end. Rownames can be freely chosen.
#'@param regionDefinitionTable defines specific regions in the genome from which the random sets are sampled. If a table is supplied,
#'the regions of interest are assigned to the defined regions (unassigned will be removed) and sampling happens within the defined region.
#'If no regions are defined, the regions of interest are assigned to the whole chromosomes. Mapping of regions is done on the whole length. 
#'Hence, a test region must be entirely within a defined region to be assigned to it - unassigned test regions are removed.
#'Note that the defined regions should not overlap (test regions can do so). The table must have three columns (chrom, start, and end).
#'User defined names can be given using an optional fourth column (name). An example is given in the A. thaliana tutorial.
#'@param simplifiedNames is a list with the column names of the annotation as keys to simplified names as values
#'(e.g. ann_transposable_element_gene can be replaced by TE-gene).
#'@param pValueThreshold specifies the significance-threshold (only significant values are drawn in the heatmaps).
#'@param countDataWasLogged specifies if the count data was already log-transformed. Generally TRUE if \code{useLog = TRUE} while loading the annotation with \code{\link{f.read.annotation}}.
#'@return a list holding the observed values, the enrichment compared to the random sets, and the corresponding P-values.
#'@examples \dontrun{
#'enrichmentResult <- f.test.regions.for.feature.enrichment.bin.based(
#'  annotation = annotationTableOnRestrictionFragments,
#'  rDir = "/path/to/where/the/results/are/stored",
#'  outfilePrefix = "aPrefixForTheFileNames",
#'  repetitions = 1e4,
#'  testRegionsTable = tableWithRegionsOfInterest,
#'  regionDefinitionTable = data.frame(),
#'  simplifiedNames = list(),
#'  pValueThreshold = 0.05,
#'  countDataWasLogged = TRUE
#') }
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@seealso \code{\link{f.test.regions.for.feature.enrichment.fragment.based}}
#'@export
f.test.regions.for.feature.enrichment.bin.based <- function(annotation, binSize, rDir, outfilePrefix, repetitions, testRegionsTable, regionDefinitionTable = data.frame(), simplifiedNames = list(), pValueThreshold = 0.05, countDataWasLogged = TRUE) {
	cat("preparing...\n")
	# check the annotation table and set the names of the features
	featuresToTest <- grep("^ann_|^sum_|^den_", colnames(annotation), value = TRUE)
	featureNames <- featuresToTest
	if (length(featuresToTest) != (ncol(annotation)-4)) { cat("note that there are one or more invalid columns (without ann_, sum_, den_ in the name)\n"); return(FALSE); }
	if (length(featuresToTest) == 1) { 
		cat("testing only one feature - duplicating it for simplicity.\n")
		dummyName <- paste(featuresToTest, "duplicate", sep = '_')
		annotation[[dummyName]] <- annotation[[featuresToTest]]
		featuresToTest <- c(featuresToTest, dummyName)
		featureNames <- featuresToTest
	}
	if (length(simplifiedNames) > 0) { withinList <- featureNames %in% names(simplifiedNames); featureNames[withinList] <- unlist(simplifiedNames)[featureNames[withinList]]; } 
	# set up the region definition (create a default if no user specified region exists, add names if not given)
	if (nrow(regionDefinitionTable) == 0) { chromSizes <- f.internal.get.chrom.sizes(); regionDefinitionTable <- data.frame(chrom = names(chromSizes), start = rep(0, length(chromSizes)), end = unlist(chromSizes), stringsAsFactors = FALSE) }
	if (ncol(regionDefinitionTable) == 3) { regionDefinitionTable$name <- paste(regionDefinitionTable$chrom, regionDefinitionTable$start, regionDefinitionTable$end, sep = '_') }
	rownames(regionDefinitionTable) <- regionDefinitionTable$name
	# map the test regions to the defined regions
	testRegionsTable$assignedTo <- "empty"
	for (i in 1:nrow(testRegionsTable)) {
		temp <- (regionDefinitionTable$chrom == testRegionsTable$chrom[i]) & (regionDefinitionTable$start <= testRegionsTable$start[i]) & (regionDefinitionTable$end >= testRegionsTable$end[i])
		testRegionsTable$assignedTo[i] <- ifelse(sum(temp) == 1, regionDefinitionTable$name[temp], "none")
	}
	numRemoved <- sum(testRegionsTable$assignedTo == "none")
	cat(paste("removed ", numRemoved, " unassigned regions: ", paste(rownames(testRegionsTable)[testRegionsTable$assignedTo == "none"], collapse = ','), '\n', sep = '')) #cat(paste("removed ", numRemoved, " unassigned regions\n", sep = ''))
	testRegionsTable <- testRegionsTable[testRegionsTable$assignedTo!="none",] #subset(testRegionsTable, assignedTo != "none")
	# get a seList
	seList <- f.get.se.list(binSize)
	# get the start and end bins of test regions
	testBinStart <- f.translate.chrom.pos.vector.to.index(testRegionsTable$chrom, testRegionsTable$start, seList, binSize)
	testBinEnd <- f.translate.chrom.pos.vector.to.index(testRegionsTable$chrom, testRegionsTable$end, seList, binSize)
	# for each region, sample <repetitions> random regions, calculate individual pValues as well
	cat("sampling and testing (# == 1 test region): ")
	seList <- f.get.se.list(binSize)
	sampledTabs <- list()
	observedTabs <- list()
	individualPvalueTabs <- list()
	enrichmentTabs <- list()
	for (i in 1:nrow(testRegionsTable)) {
		cat("#")
		testReg <- testRegionsTable$assignedTo[i]
		testSize <- testRegionsTable$end[i] - testRegionsTable$start[i]
		defChrom <- regionDefinitionTable[testReg, "chrom"]
		defStart <- regionDefinitionTable[testReg, "start"]
		defEnd <- regionDefinitionTable[testReg, "end"]
		numPossiblePos <- (defEnd-defStart)-testSize
		dataColsForSum <- grep("^sum_|^ann_", colnames(annotation), value = TRUE)
		dataColsForMean <- grep("^den_", colnames(annotation), value = TRUE)
		resultsHeader <- c(dataColsForSum, dataColsForMean)
		# calculate the real values
		real <- f.internal.summarize.binned.annotation(annotation, testBinStart[i], testBinEnd[i], dataColsForSum, dataColsForMean)
		names(real) <- resultsHeader
		observedTabs[[i]] <- real
		# sample and get feature summary
		toAdd <- floor(testSize/2)
		sampledPositions <- sample.int(numPossiblePos, repetitions)+defStart+toAdd
		sampledBinStart <- f.translate.chrom.pos.vector.to.index(defChrom, sampledPositions-toAdd, seList, binSize)
		sampledBinEnd <- f.translate.chrom.pos.vector.to.index(defChrom, sampledPositions+toAdd, seList, binSize)
		sampled <- t(apply(cbind(sampledBinStart, sampledBinEnd), 1, function(x) f.internal.summarize.binned.annotation(annotation, x[1], x[2], dataColsForSum, dataColsForMean)))
		sampledTabs[[i]] <- sampled
		# calculate individual P value
		individualPvalueTabs[[i]] <- apply(t(sampled) > real, 1, sum)/repetitions # this here is working - just be aware with the t() and this stuff - below is the safer version
		# enrichment
		enrichmentTabs[[i]] <- f.internal.calculate.region.enrichment(real, sampled, countDataWasLogged)
	}
	cat("\n")
	cat("combining results...\n")
	sampledTabs <- do.call("rbind", sampledTabs)
	observedTabs <- do.call("rbind", observedTabs)
	individualPvalueTabs <- do.call("rbind", individualPvalueTabs)
	enrichmentTabs <- do.call("rbind", enrichmentTabs)
	rownames(observedTabs) <- rownames(testRegionsTable)
	rownames(individualPvalueTabs) <- rownames(testRegionsTable)
	rownames(enrichmentTabs) <- rownames(testRegionsTable)
	colsToKeep <- colnames(sampledTabs)
	sampledMean <- aggregate(sampledTabs, by = list(sn = rep(c(1:repetitions), nrow(testRegionsTable))), mean)
	sampledMean <- sampledMean[,colsToKeep]
	combinedMean <- apply(observedTabs, 2, mean)
	combinedPvalue <- apply(t(sampledMean) > combinedMean, 1, sum)/repetitions
	combinedEnrichment <- f.internal.calculate.region.enrichment(combinedMean, sampledMean, countDataWasLogged)
	observedTabs <- rbind(observedTabs, combinedMean)
	pValueTabs <- rbind(individualPvalueTabs, combinedPvalue)
	enrichmentTabs <- rbind(enrichmentTabs, combinedEnrichment)
	cat("saving results...\n")
	# note that out is like list(observed = observedTabs, enrichment = enrichmentTabs, pValues = pValueTabs) but with the generic feature names replaced by simplifiedNames
	out <- f.internal.plot.and.save.test.regions.for.feature.enrichment.heatmap(observedTabs, enrichmentTabs, pValueTabs, pValueThreshold, rDir, outfilePrefix, featureNames)
	return(out)
}

#'@title Identify domains within chromosomes with HiCseg.
#'@param dataMatrix a HiC interaction matrix (see \code{\link{f.load.one.sample}}).
#'@param binSize size of genomic bins in bp, if set to 0, restriction fragments are used instead.
#'@param rDir a directory where the figure is stored.
#'@param outfilePrefix a prefix for the figure and table files.
#'@param minAverageDomainSize the minimal average domain size in base pairs. 
#'This argument is used to calculate the maximal number of domains within a chromosome given its size.
#'@param distributionType describes the distribution of the data: "B" is for Negative Binomial distribution, "P" is for the Poisson distribution,
#'and "G" is for the Gaussian distribution. In general, take Gaussian for normalized data and Poisson/Negative Binomial for raw data.
#'@param modelType "D" for block-diagonal and "Dplus" for the extended block-diagonal model.
#'@param useLog tells if the data shall be transformed using log2(data + 1).
#'@param regionDefinitionTable defines specific regions in the genome where domains shall be searched. 
#'If no regions are defined, domains are searched on the entire chromosomes. 
#'The table must have three columns (chrom, start, and end). User defined names can be given using an optional fourth column (name).
#'@return list with the region names as key to the output from the HiCseg function.
#'@note This is a simple wrapper function for the HiCseg package. To use it, you need to install this package first.
#'@references L\'evy-Leduc, C. and Delattre, M. and Mary-Huard, T. and Robin, S. (2014)
#'Two-dimensional segmentation for analyzing Hi-C data. \emph{Bioinformatics} \bold{30}, i386--i392.
#'@examples \dontrun{
#'domainResults <- f.identify.domains.with.HiCseg(
#'  dataMatrix = dataMatrixSampleX,
#'  binSize = 1e5, 
#'  rDir = "/path/to/where/the/results/are/stored",
#'  outfilePrefix = "aPrefixForTheFileNames",
#'  minAverageDomainSize = 1e6,
#'  distributionType = "G",
#'  modelType = "D",
#'  useLog = TRUE,
#'  regionDefinitionTable = data.frame()
#') }
#'@author Marc W. Schmid \email{marcschmid@@gmx.ch} and Stefan Grob \email{stefan@@grob.org}.
#'@export
f.identify.domains.with.HiCseg <- function(dataMatrix, binSize, rDir, outfilePrefix = "domains", minAverageDomainSize = 1e6, distributionType = "G", modelType = "D", useLog = TRUE, regionDefinitionTable = data.frame()) {
	#require("HiCseg")
	# set up the region definition (create a default if no user specified region exists, add names if not given)
	if (nrow(regionDefinitionTable) == 0) { chromSizes <- f.get.relevant.chrom.sizes(); regionDefinitionTable <- data.frame(chrom = names(chromSizes), start = rep(0, length(chromSizes)), end = unlist(chromSizes), stringsAsFactors = FALSE) }
	if (ncol(regionDefinitionTable) == 3) { regionDefinitionTable$name <- paste(regionDefinitionTable$chrom, regionDefinitionTable$start, regionDefinitionTable$end, sep = '_') }
	rownames(regionDefinitionTable) <- regionDefinitionTable$name
	# output shall be a list with the region names as key to the output from the HiCseg function
	out <- list()
	for (rn in rownames(regionDefinitionTable)) {
		rc <- regionDefinitionTable[rn, "chrom"]
		rs <- regionDefinitionTable[rn, "start"]
		re <- regionDefinitionTable[rn, "end"]
		temp <- f.extract.subset(dataMatrix, binSize, rc, rc, xStart = rs, yStart = rs, xEnd = re, yEnd = re)
		if (useLog) {temp <- log2(temp+1)}
		res <- HiCseg_linkC_R(ncol(temp), ceiling((re-rs)/minAverageDomainSize), distributionType, temp, modelType)
		if (sum(res[["t_hat"]]==0) == 0) { cat(paste("There may be more domains in the region ", rn, "\n", sep = '')) }
		# save the results of t_hat and J in a table
		out[[rn]] <- data.frame(tHat = res$t_hat, logLikelihood = res$J, chrom = rep(rc,length(res$t_hat)), boundary = rs+(res$t_ha-1)*binSize)
		out[[rn]] <- out[[rn]][out[[rn]]$tHat != 0, ]
		write.table(out[[rn]], file.path(rDir, paste(outfilePrefix, "_", rn, ".txt", sep = '')), sep = '\t', quote = FALSE, row.names = FALSE)
		# draw a figure with the domains (use correlated data plot)
		cm <- cor(f.internal.normalize.distance(temp))
		diag(cm) <- mean(cm[(row(cm)==col(cm)+1)])
		t_hat <- c(1, res$t_hat[res$t_hat!=0])
		linesMat <- matrix(0, nrow(cm), ncol(cm))
		for (i in 1:(length(t_hat)-1)) {
			cs <- t_hat[i]:t_hat[i+1]
			linesMat[t_hat[i],cs] <- 1
			linesMat[cs,t_hat[i]] <- 1
			linesMat[t_hat[i+1],cs] <- 1
			linesMat[cs,t_hat[i+1]] <- 1
		}
		# some plot parameters
		axLabelPos <- seq(1, ncol(temp), length.out = 10)
		axLabels <- (rs + (axLabelPos-1)*binSize)
		if ((re-rs) > 10e6) {axLabels <- axLabels/1e6}
		axLabels <- round(axLabels,2)
		colorSet <- colorRampPalette(c("#ffeda0", "#feb24c","#f03b20"))(64)
		if (GLOBAL_VARIABLE_USE_SVG_AND_RSVG_CONVERT) {
			svg(file.path(rDir, paste(outfilePrefix, "_", rn, ".svg", sep = '')), height = 19, width = 19)
		} else {
			tiff(file.path(rDir, paste(outfilePrefix, "_", rn, ".tiff", sep = '')), width = 2400, height = 2400)
		}
		par(oma=c(5,5,0,0), mar = c(5,5,0,0))
		image(1:ncol(cm), 1:ncol(cm), cm, col = colorSet, useRaster=TRUE, yaxt = "n", xaxt = "n", xlab = "", ylab = "")
		axis(1, at = axLabelPos, labels = axLabels, outer=FALSE, line=2, lwd=2, cex.axis=1.5, las=2)
		axis(2, at = axLabelPos, labels = axLabels, outer=FALSE, line=2, lwd=2, cex.axis=1.5, las=1)
		image(1:ncol(cm), 1:ncol(cm), linesMat, col = c(rgb(1,1,1,0), rgb(0,0,0,1)), add = TRUE)
		dev.off()
	}
	return(out)
}


