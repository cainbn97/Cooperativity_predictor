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


Indices = length(Models$Values)-1
print(paste("Indices: ",Indices))
prop = 30 # threshold for spacers/edging
thres_whole = 1 # threshold for entire motif
summary_file = paste(TF,"_summary_results.txt",sep = "")
file.create(summary_file)
expo = c("Motif", "Score", "nt_left", "nt_right")
write.table(t(expo), sep = "\t", summary_file, append = TRUE, row.names = FALSE, col.names = FALSE)

## Create logo files for each motif
for ( i in 1:Indices )
	{
	# ## create logo of original motif
	# LogoName = paste(TF,"_original",i,".png",sep="")
	# print(paste("Index number:  ",i))
	# png(filename = LogoName)
	# p = logo(model = Models, index = i)
	# print(p)
	# dev.off()

	## Extract motif information
	fit.output = Models$Values[[i]]

	## this follows motif analysis in logo.R from NRLB package
	
	
	print(paste("Index number:  ",i))
	
	
	## format motif
	motif = exp(as.numeric(fit.output$NB))
	k = length(motif)/4
	dim(motif) = c(4, k)
	motif = log(motif)
	rownames(motif) = c("A","C","G","T")
	
	
	## format error
	error = exp(as.numeric(fit.output$NE))
	dim(error) = c(4, k)
	error = log(error)
	rownames(error) = c("A","C","G","T")
	
	
	# print(paste("kmer length: ",k))
	next_req = FALSE
	
	## centralize energy motif
	motif = apply(motif, 2, function(column) column-mean(column))

	## determine the sum of the motif (strength of binding or antibinding)
	motif_abs = abs(motif)
	
	totals = apply(motif_abs, 2, function(column) sum(column))
	thres = max(totals)*prop/100
	

	## filter out bad motifs
	print(paste("Mean totals: ",mean(totals)))
	if ( mean(totals) < thres_whole || k < 3 )
		{ 
		print("Motif did not pass the total score threshold or length")
		next 
		}


	## edge beginning
	count = 0
	while ( totals[1] < thres )
		{
		count = count + 1

		# remove first kmer
		motif = motif[,2:k]
		error = error[,2:k]
		# recalculate averages vector
		motif_abs = abs(motif)
		
		totals = apply(motif_abs, 2, function(column) sum(column))
	
		# calculate new kmer length	
		k = dim(motif)[2]
		if ( k <= 3 )
			{
			print("Motif too small after trimming")
			next_req = TRUE
			break
			}

		}
	L.del = count

	## edge end
	k = length(motif)/4
	count = 0
	while ( totals[k] < thres )
		{
		count = count + 1
		## calculate new kmer length
		k = length(motif)/4

		## remove last kmer
		motif = motif[,1:k-1]
		error = error[,1:k-1]

		# recalculate averages vector
		motif_abs = abs(motif)
		totals = apply(motif_abs, 2, function(column) sum(column))
		
		# reset k value
		k = dim(motif)[2]
		if ( k <= 3 )
			{
			print("Motif too small after trimming")
			next_req = TRUE
			break
			}

		}
	
	## need to break out of two loops
	if ( next_req == TRUE )
		{next}
		
	R.del = count
	
	
	
	## spacers? - start at center and test outward for a spacer
	
	Spacer = rep("n",k)
	
	## set variables for start of the loop
	count = 0
	FOR = TRUE
	REV = TRUE
	
	for ( l in length(totals) ) 
		{
		if ( totals[l] < thres )
			{
			print(here)
			while ( FOR == TRUE || REV == TRUE )
				{
				count = count + 1
				# forward direction
				if ( totals[ l + count ] < thres && FOR == TRUE )
					{ 
					Spacer[l + count] = "y"
					# only turn this to 1 if an adjacent is also < thres
					Spacer[l] = "y"
					print("THERE")
					print(Spacer)
					}
				else
					{ 
					FOR = FALSE
					}
				if ( totals[l - count] < thres && REV == TRUE )
					{ 
					Spacer[l - count] = "y"
					# only turn this to 1 if an adjacent is also < thres
					Spacer[l] = "y"
					print("HERE")
					print(Spacer)
					}
				else
					{ 
					REV = FALSE 
					}
				}			
			}
		}	

	
	# create logo of trimmed motif
	LogoName = paste(TF,"_trimmed",i,".png",sep="")
	png(filename = LogoName)
	p = logo(model = Models, index = i, l.del = L.del, r.del = R.del )
	print(p)
	dev.off()
	
	print(error)
	print(motif)
	
	error_percent = abs(error/motif)*100 
	
	colnames(error_percent) = Spacer
	
	# find score; this is like golf currently
	score = mean(error_percent[,"n"])
	print(score)
	
	# write to a file
	output_summ = file(summary_file)
	
	## EVENTUALLY ADD NUMBER AND LENGTH OF SPACER
	
	## create dataframe with export information
	expo = c(i, round(score,0), L.del, R.del)
	write.table(t(expo), sep = "\t", summary_file, append = TRUE, row.names = FALSE, col.names = FALSE)
		
	}

	