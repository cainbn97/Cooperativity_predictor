#!/usr/bin/Rscript

## Author: cainu5
## Generates logos and analysis of NRLB results
## 09/08/20


## Load libraries
library("ggseqlogo", lib.loc="/users/cainu5/Rpackages")
library("NRLBtools", lib.loc="/users/cainu5/Rpackages")


options(error=traceback)

## Bring in arguments
args = commandArgs(trailingOnly=TRUE)
TF = args[1]

## Read output file
Filename = paste(TF,"-","NRLBConfig.csv",sep="")

## Load models
Models = load.models(fileName = Filename)


Indices = length(Models$Values)
thres = 2.5

## Create logo files for each motif
for ( i in 1 ) #:Indices)
	{

	print("Model: ", i)

	# create logo of trimmed motif
	LogoName = paste(TF,"_original",i,".png",sep="")
	png(filename = LogoName)
	p = logo(model = Models, index = i)
	print(p)
	dev.off()

	## Edge/trim the motifs
	fit.output = Models$Values[[i]]

	## CHECK WHETHER THE EXP AND LOG ARE NECESSARY
	# this follows motif analysis in logo.R from package
	motif = exp(as.numeric(fit.output$NB))
	error = exp(as.numeric(fit.output$NE))
	k = length(motif)/4
	dim(motif) = c(4, k)
	dim(error) = c(4, k)
	motif = log(motif)
	rownames(motif) = c("A","C","G","T")
	error = log(error)
	rownames(error) = c("A","C","G","T")
	
	# centralize energy motif
	motif = apply(motif, 2, function(column) column-mean(column))

	# determine the sum of the motif (strength of binding or antibinding)
	motif_abs = abs(motif)
	totals = apply(motif_abs, 2, function(column) sum(column))

	# # edge beginning
	count = 0
	while ( totals[1] < thres )
		{
		count = count + 1
		# calculate new kmer length
		k = length(motif)/4

		# remove first kmer
		motif = motif[,2:k]

		# recalculate averages vector
		motif_abs = abs(motif)
		totals = apply(motif_abs, 2, function(column) sum(column))

		}
	L.del = count

	# edge end
	k = length(motif)/4
	count = 0
	while ( totals[k] < thres )
		{
		count = count + 1
		# calculate new kmer length
		k = length(motif)/4

		# remove last kmer
		motif = motif[,1:k-1]

		# recalculate averages vector
		motif_abs = abs(motif)
		totals = apply(motif_abs, 2, function(column) sum(column))
		
		# reset k value
		k = length(motif)/4

		}
	
	R.del = count
	
	
	
	## spacers? - start at center and test outward for a spacer
	
	# 0 = not spacer; 1 = spacer
	Spacer = rep(0,k)
	count = 0
	for ( l in length(totals) ) 
		{
		if ( totals[l] < thres )
			{
			while ( FOR == TRUE || REV == TRUE )
				{
				count = count + 1
				# forward direction
				if ( totals[ l + count ] < thres && FOR == TRUE )
					{ 
					Spacer[l + count] = 1
					# only turn this to 1 if an adjacent is also < thres
					Spacer[l] = 1
					}
				else
					{ FOR = FALSE }
				if ( totals[l - count] < thres  && REV = TRUE )
					{ 
					Spacer[l - count] = 1 
					# only turn this to 1 if an adjacent is also < thres
					Spacer[l] = 1
					}
				else
					{ REV = FALSE }
				}			
			}
		}	

	

	# create logo of trimmed motif
	LogoName = paste(TF,"_trimmed",i,".png",sep="")
	png(filename = LogoName)
	p = logo(model = Models, index = i, l.del = L.del, r.del = R.del )
	print(p)
	dev.off()
	
	
	# score the motif 
	
	
	
	

	}

