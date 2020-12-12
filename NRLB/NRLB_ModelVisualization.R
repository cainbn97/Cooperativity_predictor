#!/usr/bin/Rscript

## Author: cainu5
##
## Filters out weak motifs
## Edges motifs so core motif remains
## Defines spacers
## Creates confidence score for motif
## Exports scores, logos of original and trimmed motifs, PWMS and list of 
## motif exports for DiffLogo clustering
## TODO - cluster motifs with DiffLogo so only unique motifs are kept
## TODO - define number and length of spacers and add to output
##
## 09/18/20


## Load libraries
library("ggseqlogo", lib.loc="/users/cainu5/Rpackages")
library("NRLBtools", lib.loc="/users/cainu5/Rpackages")
# library("proxy", lib.loc="/users/cainu5/Rpackages")
# library("cba", lib.loc="/users/cainu5/Rpackages")
# library("DiffLogo", lib.loc="/users/cainu5/Rpackages")

## Create an error traceback
options(error=traceback)

## Bring in arguments
args = commandArgs(trailingOnly=TRUE)
TF = args[1]

## Determine csv file name
Filename = paste(TF,"-","NRLBConfig.csv",sep="")

## Load models
Models = load.models(fileName = Filename)
Indices = length(Models$Values)-1

## Define thresholds
prop = 30 # threshold for spacers/edging
thres_whole = 1 # threshold for entire motif

## Create summary file - append is not true, so fresh file is created
summary_file = paste(TF,"_summary_results.txt",sep = "")
expo = c("Motif", "Score", "nt_left", "nt_right","Spacer?")
write.table(t(expo), sep = "\t", summary_file, row.names = FALSE, 
	col.names = FALSE, quote = FALSE )
	
## Create error log for filtered motifs
error_file = paste(TF,"_error.txt", sep = "")
expo_err = c("Motif", "Reason")
write.table(t(expo_err), sep = "\t", error_file, row.names = FALSE,
	col.names = FALSE, quote = FALSE )

for ( i in 1:Indices )
	{
	
	## Extract motif information
	fit.output = Models$Values[[i]]

	print(paste("Index number:  ",i))
	
	## format motif - like NRLB logo script
	motif = exp(as.numeric(fit.output$NB))
	k = length(motif)/4
	dim(motif) = c(4, k)
	motif = log(motif)
	rownames(motif) = c("A","C","G","T")
	
	## format error - like NRLB logo script
	error = exp(as.numeric(fit.output$NE))
	dim(error) = c(4, k)
	error = log(error)
	rownames(error) = c("A","C","G","T")
	
	## set variable for if need to break out of two loops
	next_req = FALSE
	
	## centralize energy motif around x-axis
	motif = apply(motif, 2, function(column) column-mean(column))

	## find sum of each column of absolute value of matrix - strength of binding
	## and antibinding
	motif_abs = abs(motif)
	totals = apply(motif_abs, 2, function(column) sum(column))
	thres = max(totals)*prop/100

	## filter out bad motifs
	if ( mean(totals) < thres_whole )
		{ 
		err = paste("Motif has an average column weight of", round(mean(totals),3))
		expo_err = c(i, err)
		write.table(t(expo_err), sep = "\t", error_file, append = TRUE, 
			row.names = FALSE, col.names = FALSE, quote = FALSE)
		next 
		}
		
	if ( k < 3 )
		{
		expo_err = c(i,"Motif is less than 3 nucleotides")
		write.table(t(expo_err), sep = "\t", error_file, append = TRUE, 
			row.names = FALSE, col.names = FALSE, quote = FALSE)
		next
		}

	## edge left
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
			expo_err = c(i,"Motif too small after trimming")
			write.table(t(expo_err), sep = "\t", error_file, append = TRUE, 
				row.names = FALSE, col.names = FALSE, quote = FALSE )
			next_req = TRUE
			break
			}
		}

	## need to break out of two loops
	if ( next_req == TRUE )
		{next}

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
			expo_err = c(i,"Motif too small after trimming")
			write.table(t(expo_err), sep = "\t", error_file, append = TRUE, 
				row.names = FALSE, col.names = FALSE,quote = FALSE)
			next_req = TRUE
			break
			}
		}
	
	## need to break out of two loops
	if ( next_req == TRUE )
		{next}
	R.del = count
	
	## set variables for start of the loop
	Spacer = rep("n",k)
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

	## create logo of trimmed motif
	LogoName = paste(TF,"_trimmed",i,".png",sep="")
	png(filename = LogoName)
	p = logo(model = Models, index = i, l.del = L.del, r.del = R.del )
	print(p)
	dev.off()

	
	## export motif to text file for diffLogo clustering
	motif_name = paste(TF, "_motif_",i,".txt",sep = "")
	write.table(motif, sep = "\t", motif_name, col.names = FALSE, 
		row.names = FALSE)
	
	## export motif name to motif list for diffLogo readout
	motif_list = paste(TF,"_motif_list.txt")
	write.table(motif_name, sep = "\n", motif_list, 
		row.names = FALSE, col.names = FALSE, append = TRUE, quote = FALSE)
	
	error_percent = abs(error/motif)*100 
	
	## label columns as spacer or not spacer
	colnames(error_percent) = Spacer
	
	## find score of non-spacers; this is like golf currently
	score = mean(error_percent[,"n"])
	
	## create dataframe with export information
	expo = c(i, round(score,0), L.del, R.del, any(Spacer == "y"))
	write.table(t(expo), sep = "\t", summary_file, append = TRUE, 
		row.names = FALSE, col.names = FALSE)

	}
	
	
# ## Prep motifs for clustering
# motifs = readLines(motif_list)
# print(motifs)
# pwm = list()
# for ( m in length(motifs) )
	# {
	# print(m)
	# print(motifs[m])
	# pwm[[m]] = getPwmFromPwmFile(motifs[m])
	# print(pwm[[m]])
	# }

# # diffLogoTable(motif_list)

	