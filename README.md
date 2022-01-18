# Predict homeodomain transcription factor cooperativity with HT-SELEX

This program predicts cooperativity of transcription factors using HT-SELEX by both assessing dimer versus individual half-site enrichment and spacer preference as done in *Link paper here!*.

We suggest running both the [`enrichment analysis`](#Enrichment analysis) and [`COSMO analysis`](#COSMO analysis) as they test different aspects of cooperativity. 

## Prerequisites

*Prerequisite versions that program was tested on are listed. Other versions may be usable.*
1. [HOMER 4.9](http://homer.ucsd.edu/homer/introduction/install.html)
	- Requires below line changes in findKnownMotifs.pl to eliminate divide by zero error. 


| Line number | Original                                                                     | Required new line
| 152 --------| $totalNumTargets = floor($numTargets/($percentTargets\*0.01)+0.5);           | $totalNumTargets = floor($numTargets/($percentTargets\*0.01+0.0005));
| 153         | $totalNumBackground = floor($numBackground/($percentBackground\*0.01)+0.5);  | $totalNumBackground = floor($numBackground/($percentBackground\*0.01+0.0005));



2. COSMO (will need to link)
	- Python 2.7
		- [MOODS 1.0.2.1](https://www.cs.helsinki.fi/group/pssmfind/)
		- numpy 1.16.6
		- scipy 1.2.3
	- We suggest creating a virtual environment with these dependencies. 
3. Python 3.7.1
	- os
	- numpy
	- glob
	- pandas 
	- re
	- argparse
	- matplotlib
	- scipy
	- seaborn
	- bs4
4. R 3.2.0
	- [readr](https://cran.r-project.org/web/packages/readr/index.html)
	- [outliers](https://cran.r-project.org/web/packages/outliers/index.html)
	- [chisq.posthoc.test](https://cran.r-project.org/web/packages/chisq.posthoc.test/index.html)

## Usage

User can select to download HT-SELEX data associated with [DNA-Binding Specificities of Human Transcription Factors](http://dx.doi.org/10.1016/j.cell.2012.12.009) for analysis or use their own data. 

### Downloading HT-SELEX data from Jolma 2013

```bash
Cooperativity_predictor.sh OPTIONS [ANALYSIS NAME] [Download link of Cycle 1]
```

For example,
```bash
Cooperativity_predictor.sh GSX2 ftp.sra.ebi.ac.uk/vol1/run/ERR195/ERR195221/GSX2_TCCAAC20NCG_Y_1.fastq.gz
```
The program will automatically determine the download links for the remaining cycles. The initial library download link must be inputted separately.

### Using downloaded or own HT-SELEX data

```bash
Cooperativity_predictor.sh OPTIONS [ANALYSIS NAME]
```

Program expects the files to be in the following layout:

```
	Analysis name	
	├── Cycle1	
	│   └── \[Cycle1\].fastq.gz		
	├── Cycle2	
	│   └── \[Cycle2\].fastq.gz	
	├── Cycle3	
	│   └── \[Cycle3\].fastq.gz	
	├── Cycle4	
	   └── \[Cycle4\].fastq.gz

```

### Program run options

Program run options:
Default: all analyses run
\[-h\] HOMER de novo motif analysis and dimer site selection run.
\[-e\] Enrichment analysis only run, requires Homer de novo motif files.
\[-c\] COSMO only run, requires Homer de novo motif files.

Additional options:
\[-z\] Initial library cycle download link for normalization. If no download link
		 or file \[-b\] provided, program uses cycle 1 library for normalization.
\[-b\] Enter a fastq.gz file from the initial library of HT-SELEX. If no download link \[-z\]
		 or file provided, program uses cycle 1 library for normalization
\[-t\] Motif threshold for MOODS. Must be a value between 0 and 1.
		Default is 0.8.
\[-p\]: Enter a motif file if you would like to use a motif generator
		other than HOMER. This will overwrite -h option.

## Program outputs

### Dimer search

1. dimer_description_check.txt
	- Gives Site:Non-site and site specific nucleotide impacts for each dimer site
2. *If COSMO analysis specificed* top_dimer_kmer_motifs_\[motif\]
	- Contains Jaspar formatted PWMs for COSMO analysis
3. long_motif_consensus.txt
	- This is more of an intermediate file. It contains the found dimer sites. 

**Output organization**

&emsp; Analysis
&emsp;├── top_dimer_kmer_motifs_\[motif\]
&emsp;└── motif1.jpwm
&emsp;│	└── motif2.jpwm
&emsp;├── dimer_description_check.txt
&emsp;├── long_motif_consensus.txt
&emsp;...
	


*Program will call a motif a dimer site if nucleotide impact for each site is greater than 0.6 and the site:non-site ratio is greater than 1.6.*

### Enrichment analysis

1. Folders containing known motif results for each half site and dimer site at each cycle
2. Analysis_Enrichment_analysis_run_summary.txt
3. \[Analysis\]\_NatLog_2_Enrichment_plot_\[motif\].txt

*Any motif with an enrichment factor of greater than 2 is indicative of cooperative behavior.*

**Output organization**
	
	Analysis

	├── Cycle1

	│   └── GSX2_1_dimer1_homer

	│	└── GSX2_1_site1_dimer1_homer

	│	└── GSX2_1_site2_dimer1_homer

	├── Analysis_Enrichment_analysis_run_summary.txt

	├── \[Analysis\]\_NatLog_2_Enrichment_plot_\[motif\].txt

	...

	
**Example of enrichment plot**

![Enrichment_sample](Enrichment_plot_sample.png)

### COSMO analysis

1. *Analysis*\_COSMO_run_summary_wgrubbs.txt
2. Raw COSMO output and tabbed files separated by motif combination for each cycle and dimer site motif
3. Dimer site counts found by COSMO for each motif combination for each dimer site motif
	- x-axis: spacer length
	- y-axis: cycle number
4. Heatmap illustrating dimer counts across cycles and spacer lengths

*Any motif with a p-value of < 0.05 for the chi-square test for independence and Grubb's test is indicative of cooperative behavior.*

**Output organization**

	Analysis

	├── top_dimer_kmer_motifs_dimer1

	│   └── *Analysis*\_dimer1_motif1_motif1_FF_cosmo_output.png

	│	└── *Analysis*\_dimer1_motif1_motif2_FF_cosmo_output.png

	│	└── *Analysis*\_dimer1_motif2_motif1_FF_cosmo_output.png

	│	└── *Analysis*\_dimer1_motif2_motif2_FF_cosmo_output.png

	│	└── motif1.jpwm

	│	└── motif2.jpwm

	│	└── *Analysis*\_dimer1_motif1_motif1_FF.txt

	│	└── *Analysis*\_dimer1_motif1_motif2_FF.txt

	│	└── *Analysis*\_dimer1_motif2_motif1_FF.txt

	│	└── *Analysis*\_dimer1_motif2_motif2_FF.txt

	├── Cycle1

	│   └── Analysis_1_dimer1_homer

	│		└── cosmo.counts.tab

	│		└── Cycle1_motif1_motif1_FF.tab

	│		└── Cycle1_motif1_motif2_FF.tab

	│		└── Cycle1_motif2_motif1_FF.tab

	│		└── Cycle2_motif2_motif1_FF.tab

	│	└── GSX2_1_site1_dimer1_homer

	│	└── GSX2_1_site2_dimer1_homer

	...

	
**Example of COSMO heatmap**

![COSMO_heatmap](COSMO_heatmap_sample.png)

## Authors

| Contributor                       | Institution                 | Remarks
|-----------------------------------|-----------------------------|-------------------------
| [Brittany Cain](mailto:Brittany.Cain@cchmc.org)       | Cincinnati Children's Hosp.   | Pipeline author



## License

MIT. See [`LICENSE.txt`](LICENSE.txt).
