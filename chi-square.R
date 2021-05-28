#!/usr/bin/Rscript

## Author: cainu5
##
## Runs chi-square and post-hoc test based on COSMO counts
##
## 12/23/20

library("chisq.posthoc.test", lib.loc = "/users/cainu5/Rpackages")
library("readr")
library("outliers")

## Create an error traceback
options(error=traceback)

## Read in table
file = list.files(pattern = 'COSMO_counts_motif1_motif2')
	
counts = read.table(file, sep = "")
Cycle0_4 = counts[c("0","4"),]

## Read in long_consensus sequence from SELEX datasets
TF = basename(getwd())
consensus_path = paste(getwd(),"/long_motif_consensus.txt", sep = "")

## Grab spacer value from string
long_consensus = read_file(consensus_path)
print(long_consensus, quote = FALSE)
spacer = parse_number(long_consensus)
Column = paste("X",spacer, sep = "")

Cycle0_4_norm = data.frame(matrix(0, nrow = 2, ncol = dim(counts)[2]), row.names = c("0","4"))

## Make counts percentage of total dimers at each cycle - look at proportions of dimers not enrichment
for ( j in 1:dim(counts)[2] ) 
{
	for ( i in c("0","4") )
	{
		Cycle0_4_norm[i,j] = Cycle0_4[i,j]/sum(Cycle0_4[i,])*100
	}
	if ( Cycle0_4_norm["0", j] < Cycle0_4_norm["4",j] )
	{
		print(paste("Dimers with a spacer of ", j, "is overrepresented in Cycle 4."), quote = FALSE)
	}
}

## Perform chi-square
chi_tot = chisq.test(Cycle0_4_norm, y = NULL, rescale.p = TRUE)
print(chi_tot)
chi_tot_pvalue = chi_tot$p.value

## Perform post-hoc
chi_pair_bon = chisq.posthoc.test(Cycle0_4_norm, method = "bonferroni")
print(chi_pair_bon)
chi_pair = chisq.posthoc.test(Cycle0_4_norm, method = "none")
print(chi_pair)

## Print relevant p-value based on found spacer from Homer SELEX analysis
chi_spac_bon = subset(chi_pair_bon, ( chi_pair$Value == "p values" & chi_pair$Dimension == 0 ))[Column]
chi_spac = subset(chi_pair, ( chi_pair$Value == "p values" & chi_pair$Dimension == 0 ))[Column]

## run grubbs test
grubbs_result = grubbs.test(as.numeric(Cycle0_4["4",]), type = 10, opposite = FALSE)
G = grubbs_result$statistic["G"]
U = grubbs_result$statistic["U"]
p_grubbs = grubbs_result$p.value
alt = grubbs_result$alternative
top_spac = parse_number(colnames(Cycle0_4["4",][which(Cycle0_4["4",] == parse_number(alt))]))

## Read in statistics file and print out relevant p-value
stat_file = list.files(pattern ='_stats.txt')
stats = read.table(paste('Cycle4/',TF,'_4_stats.txt',sep = ""), sep = "\t", skip = 1, header = TRUE)
spac_orient_p = paste('motif1.jpwm|motif2.jpwm|FF|',spacer,sep = "")
COSMO_zstats = subset(stats, ( stats$TF1.TF2..F.R..D == spac_orient_p))[2:7]
p_spac = subset(stats, ( stats$TF1.TF2..F.R..D == spac_orient_p))[8]

## Write data to file
COSMO_output_file = paste(getwd(),'/',TF,'_COSMO_run_summary_wgrubbs.txt', sep = "")
sink(COSMO_output_file, append = TRUE)
cat(paste('TF','Dimer site','Chi-square', 'Chi-square post hoc no correction', 'Chi-square post hoc bonferroni correction', 'Grubbs test', 'Top spacer found', 'Z-test result (p)', sep = "\t"))
cat("\n")
cat(paste(TF, long_consensus, round(chi_tot_pvalue,4) , chi_spac, chi_spac_bon, round(p_grubbs,4), top_spac, p_spac, sep = "\t"))
cat("\t")
cat(unlist(COSMO_zstats),sep = "\t")
for ( i in c("0","1","2","3","4") )
	{
	cat("\t")
	cat(unlist(counts[i,]), sep="\t")
	}
cat("\n")
sink()

