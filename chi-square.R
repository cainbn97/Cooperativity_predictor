#!/usr/bin/Rscript

## Author: cainu5
##
## Runs chi-square and post-hoc test based on COSMO counts
##
## 06/10/21

library("chisq.posthoc.test")
library("readr")
library("outliers")

## Create an error traceback
options(error=traceback)

TF = basename(getwd())
c = 1
setwd('PWMs_of_dimers')
dimers = list.files(pattern = '^dimer*')
# print(dimers)

## Prep output log
setwd("..")
COSMO_output_file = paste(getwd(),'/',TF,'_COSMO_run_summary_wgrubbs.txt', sep = "")
sink(COSMO_output_file, append = FALSE)
cat(paste('TF','PWM','Dimer site','Chi-square', 'Chi-square post hoc no correction', 'Chi-square post hoc bonferroni correction', 'Grubbs test', 'Top spacer found', sep = "\t"))
cat("\n")


if ( readLines('long_motif_consensus.txt')[1] == 'N/A')
{
	cat(paste(TF,'No dimer site found', sep = "\t"))
	cat('\n')
} else if ( length(dimers) == 0 ) {
	sink(COSMO_output_file, append = TRUE)
	cat(paste(TF, 'Dimer prevalence < 5%', sep = "\t"))
	cat('\n')
} 

sink()

for ( dimer in sort(dimers) )
	{
	print(dimer)
	dimer_name = sub('\\.motif$', '', dimer)
	## Read in table
	pattern2search = paste(dimer_name, '_', 'COSMO_counts_motif1_motif2_FF.txt', sep = "")
	File = list.files(pattern = pattern2search)
	counts = read.table(File, sep = "")
	Cycle0_4 = rbind( head(counts, n = 1), tail(counts, n = 1) )

	## Read in long_consensus sequence from SELEX datasets
	consensus_path = paste(getwd(),"/long_motif_consensus.txt", sep = "")

	## Grab spacer value from string
	long_consensus = readLines(consensus_path)[c]
	print(long_consensus, quote = FALSE)
	spacer = parse_number(long_consensus)
	Column = paste("X",spacer, sep = "")
	
	Cycle0_4 = data.frame(Cycle0_4, row.names = c("start","end"))
	print(Cycle0_4)
	Cycle0_4_norm = data.frame(matrix(0, nrow = 2, ncol = dim(counts)[2]), row.names = c("start","end"))

	## Make counts percentage of total dimers at each cycle - look at proportions of dimers not enrichment
	for ( j in 1:dim(counts)[2] ) 
	{
		for ( i in c("start","end") )
		{
			Cycle0_4_norm[i,j] = Cycle0_4[i,j]/sum(Cycle0_4[i,])*100
		}
		if ( Cycle0_4_norm["start", j] < Cycle0_4_norm["end",j] )
		{
			print(paste("Dimers with a spacer of ", j, "is overrepresented in Cycle 4."), quote = FALSE)
		}
	}
	print(Cycle0_4_norm)
	## Perform chi-square
	chi_tot = chisq.test(Cycle0_4_norm, y = NULL, rescale.p = TRUE)
	print(chi_tot)
	chi_tot_pvalue = chi_tot$p.value

	## Perform post-hoc
	chi_pair_bon = chisq.posthoc.test(Cycle0_4_norm, method = "bonferroni", round = 30)
	print(chi_pair_bon)
	chi_pair = chisq.posthoc.test(Cycle0_4_norm, method = "none")
	print(chi_pair)

	## Print relevant p-value based on found spacer from Homer SELEX analysis
	chi_spac_bon = subset(chi_pair_bon, ( chi_pair$Value == "p values" & chi_pair$Dimension == "start" ))[Column]
	chi_spac = subset(chi_pair, ( chi_pair$Value == "p values" & chi_pair$Dimension == "start" ))[Column]

	## run grubbs test
	grubbs_result = grubbs.test(as.numeric(Cycle0_4["end",]), type = 10, opposite = FALSE)
	G = grubbs_result$statistic["G"]
	U = grubbs_result$statistic["U"]
	p_grubbs = grubbs_result$p.value
	alt = grubbs_result$alternative
	top_spac = parse_number(colnames(Cycle0_4["end",][which(Cycle0_4["end",] == parse_number(alt))]))

	## Read in statistics file and print out relevant p-value
	# stat_file = list.files(pattern ='_stats.txt')
	# stats = read.table(paste('Cycle4/',TF,'_4_stats.txt',sep = ""), sep = "\t", skip = 1, header = TRUE)
	# spac_orient_p = paste('motif1.jpwm|motif2.jpwm|FF|',spacer,sep = "")
	# COSMO_zstats = subset(stats, ( stats$TF1.TF2..F.R..D == spac_orient_p))[2:7]
	# p_spac = 'N/A' #subset(stats, ( stats$TF1.TF2..F.R..D == spac_orient_p))[8]

	## Write data to file
	sink(COSMO_output_file, append = TRUE)
	cat(paste(TF, dimer_name, long_consensus, round(chi_tot_pvalue,4) , chi_spac, chi_spac_bon, p_grubbs, top_spac, sep = "\t"))
	cat("\t")
	# cat(unlist(COSMO_zstats),sep = "\t")
	cat("\n")
	sink()
	
	c = c + 1
	}
