#!/usr/bin/Rscript

## Author: cainu5
##
## Runs chi-square and post-hoc test based on COSMO counts
##
## 12/23/20

library("chisq.posthoc.test", lib.loc = "/users/cainu5/Rpackages")
library("readr")

## Create an error traceback
options(error=traceback)

## Read in table
file = list.files(pattern = 'COSMO_counts_motif1_motif2')

## NEED SOME KIND OF FILE CHECK HERE
# if ( length(file) == 0 )
	# {
	# exit()
	# }
	
counts = read.table(file, sep = "")
Cycle0_4 = counts[c("0","4"),]

## Read in long_consensus sequence from SELEX datasets
TF = basename(getwd())
consensus_path = paste("/users/cainu5/SELEX_analysis/testing/",TF,"/long_motif_consensus.txt", sep = "")

## Grab spacer value from string
long_consensus = read_file(consensus_path)
print(long_consensus, quote = FALSE)
spacer = parse_number(long_consensus)
Column = paste("X",spacer, sep = "")

Cycle0_4_norm = data.frame(matrix(0, nrow = 2, ncol = dim(counts)[2]), row.names = c("0","4"))

## Make counts percentage of total dimers at each cycle - normalization
for ( j in 1:dim(counts)[2] ) 
{
	for ( i in c("0","4") )
	{
		Cycle0_4_norm[i,j] = Cycle0_4[i,j]/sum(Cycle0_4[i,])*100
	}
	if ( Cycle0_4_norm["0", j] < Cycle0_4_norm["4",j] )
	{
		print(paste(paste("Dimers with a spacer of ", j), "is overrepresented in Cycle 4."), quote = FALSE)
	}
}

## Perform chi-square
chi_tot = chisq.test(Cycle0_4_norm, y = NULL, rescale.p = TRUE)
print(chi_tot)

## Perform post-hoc - no multiple corrections currently being used
chi_pair = chisq.posthoc.test(Cycle0_4_norm, method = "none")
print(chi_pair)

## Print relevant p-value based on found spacer from Homer SELEX analysis
print(subset(chi_pair, ( chi_pair$Value == "p values" & chi_pair$Dimension == 0 ))[Column])
