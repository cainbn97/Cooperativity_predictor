#!/usr/bin/Rscript

## Author: cainu5
##
## Runs chi-square on all motif arrangments to find significant
## Set for monomer runs
##
## 01/09/20

library("chisq.posthoc.test", lib.loc = "/users/cainu5/Rpackages")
library("readr")

## Create an error traceback
options(error=traceback)

## Read in table
file_list = sort(list.files(pattern = 'COSMO_counts'))

## Read in long_consensus sequence from SELEX datasets
TF = basename(getwd())
consensus_path = paste("/users/cainu5/SELEX_analysis/testing/",TF,"/monomer_site.txt", sep = "")

## Grab spacer value and motifs from string
long_consensus = read_file(consensus_path)
motif1 = substr(long_consensus,1,4); print(motif1, quote = FALSE)
# motif2 = substr(long_consensus,7,10); print(motif2, quote = FALSE)
# print(long_consensus, quote = FALSE)
# spacer = parse_number(long_consensus)
# Column = paste("X",spacer, sep = "")

## print TF information
COSMO_output_file = "/users/cainu5/SELEX_analysis/COSMO_output/COSMO_run_summary_all_motif_arrangments_0.8.txt"
sink(COSMO_output_file, append = TRUE)
cat(paste(TF, motif1, sep = "\t"))
sink()

for ( f in 1:length(file_list))
	{
	print(file_list[f], quote = FALSE)
	
	## parse motif arrangment from file name
	split_string = unlist(strsplit(file_list[f],"_"))
	slot_1 = split_string[4]
	slot_2 = unlist(strsplit(split_string[5],".txt"))[1]
		
	## Parse matrix from heatmap script
	counts = read.table(file_list[f], sep = "")
	Cycle0_4 = counts[c("0","4"),]

	Cycle0_4_norm = data.frame(matrix(0, nrow = 2, ncol = dim(counts)[2]), row.names = c("0","4"))

	## Make counts percentage of total dimers at each cycle - look at proportions of dimers not enrichment
	for ( j in 1:dim(counts)[2] ) 
	{
		for ( i in c("0","4") )
		{
			Cycle0_4_norm[i,j] = Cycle0_4[i,j]/sum(Cycle0_4[i,])*100
		}
	}
	
	## Perform chi-square
	chi_tot = chisq.test(Cycle0_4_norm, y = NULL, rescale.p = TRUE)
	print(chi_tot)

	## Perform post-hoc with bonferroni corrections
	chi_pair_bon = chisq.posthoc.test(Cycle0_4_norm, method = "bonferroni")
	print(chi_pair_bon)

	## Print relevant significant bonferroni p-values
	chi_spac_bon = subset(chi_pair_bon, ( chi_pair_bon$Value == "p values" & chi_pair_bon$Dimension == 0 ))[3:12]
	
	## Print out '-' if cycle 0 was more prevalent
	for ( j in 1:dim(counts)[2] )
	{
		if ( Cycle0_4_norm["0",j] > Cycle0_4_norm["4",j] )
		{ 
			print('here')
			chi_spac_bon[j] = '-'
		}
	}
	print(Cycle0_4_norm)
	chi_min_bon = chi_spac_bon[which(chi_spac_bon < 0.05)]
	print(chi_min_bon)

	## Write data to file
	sink(COSMO_output_file, append = TRUE)
	cat("\t")
	cat(paste(slot_1, slot_2, sep = "\t"))
	for ( i in 1:length(chi_spac_bon) )
		{
		cat("\t")
		# cat(colnames(chi_min_bon)[i], sep = "\t")
		cat(unlist(chi_spac_bon[i]), sep="\t")
		}
	cat("\t")
	sink()

	}
sink(COSMO_output_file, append = TRUE)
cat("\n")
sink()