#!/bin/bash 

# Author: cainu5
# SELEX analysis for cooperativity
# 05/28/21

## Still need to add monomer COSMO run if no dimer is present ##

## Grab important file paths
BASEDIR="$(pwd)"
CODEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

## Implement a switch to allow COSMO and enrichment script to be run separately
cflag='';eflag='';hflag=''
	
while getopts "ceht:m:" opt; do
	case $opt in
	c) cflag=true 
		echo "	COSMO analysis selected"
		;;
	e) eflag=true 
		echo "	Enrichment analysis selected"
		;;
	h) hflag=true 
		echo "	HOMER de novo motif analysis selected"
		;;
	t)
		THRES=${OPTARG}
		(( t >= 0 && t <= 1 )) || 
		echo "	[-t] Enter a value between 0 and 1 for moods thresholding. Default is 0.8" 
		;;
	m)
		MODE=${OPTARG}
		(( m == 1 || m == 2)) || 
		echo "	[-m]: Enter 1 for a monomer to monomer COSMO run. Enter 2 to 
				run COSMO on the top two 4 mers from dimer site. Default is 2"
		;;
	?)
		echo 	"Script usage:
			Program run options
			[-c] COSMO run
			[-h] HOMER de novo motif analysis only run
			[-e] Enrichment analysis run
			Default: all analyses run
			[-t] Motif threshold for MOODS. Must be a value between 0 and 1. 
					Default is 0.8.
			[-m] COSMO motif run option. Enter 1 for a monomer to monomer COSMO run. Enter 2 to 
					run COSMO on the top two 4 mers from dimer site. Default is 2."
		;;
	esac
done

shift $((OPTIND-1))

## Bring in arguments
TF=${1:?Expected TF as argument #1}
Target_link=${2:?Expected ENA ftp download link for cycle 1 as argument #2}
Zero_link=${3:?Expected ENA ftp download link for cycle 0 as argument #3}

cd "$BASEDIR"

## Set flag defaults if no program options were selected
if [ "$cflag" != "true" ] && [ "$eflag" != "true" ] && [ "$hflag" != "true" ]
then
	echo "	De novo motif, enrichment, and COSMO analysis will be performed."
	cflag=true
	eflag=true
	hflag=true
fi

#####################################
## ACQUIRE FASTQ FROM ENA DATABASE ##
#####################################

## Check if folder already exists
if [ -e "$TF"/Cycle4/*.fastq.gz ]
then
	echo "	FASTQ files exist. Starting analysis..."
	## Could be cool to not require links if fastq files are present
	
	
else
	echo "	FASTQ files do not exist. Downloading files..."
	
	## Make necessary directories and download files
	mkdir "$TF"; cd "$TF"
	mkdir Cycle1; mkdir Cycle2; mkdir Cycle3; mkdir Cycle4

	## Read sample accession of cycle 1 and extrapolate other cycle accessions
	SampAcc_Cycle1=$( echo "$Target_link" | cut -d / -f 5 | cut -d R -f 3 )
	SampAcc_Cycle2=$(( SampAcc_Cycle1 + 1 ))
	SampAcc_Cycle3=$(( SampAcc_Cycle1 + 2 ))
	SampAcc_Cycle4=$(( SampAcc_Cycle1 + 3 ))

	## Piece together target links
	LINK1="$Target_link"
	LINK2=$( echo $Target_link |\
		cut -d / -f 1,2,3,4 )"/ERR${SampAcc_Cycle2}/"$( echo $Target_link |\
		cut -d / -f 6 | cut -d _ -f 1,2,3 )"_2.fastq.gz"
	LINK3=$( echo $Target_link |\
		cut -d / -f 1,2,3,4 )"/ERR${SampAcc_Cycle3}/"$( echo $Target_link |\
		cut -d / -f 6 | cut -d _ -f 1,2,3 )"_3.fastq.gz"
	LINK4=$( echo $Target_link |\
		cut -d / -f 1,2,3,4 )"/ERR${SampAcc_Cycle4}/"$( echo $Target_link |\
		cut -d / -f 6 | cut -d _ -f 1,2,3 )"_4.fastq.gz"


	cd "$BASEDIR"/"$TF"
	cd Cycle1; curl -flOJ "$Target_link"
	cd ../Cycle2; curl -flOJ "$LINK2"
	cd ../Cycle3; curl -flOJ "$LINK3"
	cd ../Cycle4; curl -flOJ "$LINK4"

fi

## Check file integrities and grab downloaded fastq file names in a vector
Rounds=(1 2 3 4)
for Cycle in ${Rounds[@]}
do
	cd "$BASEDIR"/"$TF"/"Cycle${Cycle}"
	Fastq_target["$Cycle"]=$( ls *"${Cycle}.fastq.gz" )
done

# check file integrities
cd "$BASEDIR"
for Cycle in ${Rounds[@]}
do
	if gzip -t "$TF"/"Cycle${Cycle}"/"${Fastq_target[$Cycle]}"
	then
		echo "	${Fastq_target[$Cycle]} downloaded successfully."
	else
		echo '	Download of ${Fastq_target["$Cycle"]} was unsuccessful'
		rm -f "${Fastq_target[$Cycle]}"
		exit 1
	fi
done

## Download cycle 0 file if necessary
cd "$BASEDIR"
ZeroTag=$( echo $Zero_link | cut -d '_' -f 2 )
cd "$BASEDIR"

if [ -e */"ZeroCycle_${ZeroTag}_0_0.fastq.gz" ]
then
	echo ""
else
	mkdir "$ZeroTag"; cd "$ZeroTag"
	curl -flOJ "$Zero_link"
fi

cd "$BASEDIR"/"$ZeroTag"
Fastq_bg=$( ls *.fastq.gz )
if gzip -t "${Fastq_bg}"
then
	echo "	${Fastq_bg} downloaded successfully."
else 
	echo "	Download of Zero tag file was unsuccessful"
	cd ..
	rm -rf "$ZeroTag"
	exit 1
fi

#####################################
###### HOMER de novo analysis  ######
#####################################

if [ "$hflag" = "true" ]
then
	## Load important modules
	module load homer/4.9-wrl
	module load python3

	## Convert zipped fastq files to fasta for zero cycle
	cd "$BASEDIR"/"$ZeroTag"
	zcat "$Fastq_bg" > "${ZeroTag}.fastq"
	paste - - - - < "${ZeroTag}.fastq" | cut -f 1,2 | sed 's/^@/>/' |\
	tr "\t" "\n" > "${ZeroTag}.fa"
	rm "${ZeroTag}.fastq"

	## Convert zipped fastq files to fasta for cycle 4, extract 50k sequences (sampling)
	cd "$BASEDIR"/"$TF"/"Cycle4"
	zcat "${Fastq_target[4]}" > "${TF}_${ZeroTag}_4.fastq" 
	paste - - - - < "${TF}_${ZeroTag}_4.fastq" | shuf | cut -f 1,2 |\
	sed 's/^@/>/' | tr "\t" "\n" | head -100000 > "${TF}_${ZeroTag}_4.fa"
	rm "${TF}_${ZeroTag}_4.fastq"

	## De novo motif analysis of cycle 4 short
	echo "	Starting de novo motif analysis for Cycle3 for 6-8 bp sequences"
	findMotifs.pl "${TF}_${ZeroTag}_4.fa" fasta "${TF}_4_homer_denovo_short" \
	-fasta "$BASEDIR"/"${ZeroTag}"/"${ZeroTag}.fa" \
	-noredun -len 6,8 -noknown -p 4

	python "$CODEDIR"/htmltotext.py \
		"${TF}_4_homer_denovo_short"/homerResults.html \
		"${TF}_4_homer_denovo_short"/homerResults.txt 

	## De novo motif analysis of cycle 4 long
	echo "	Starting de novo motif analysis for Cycle4 for 16-18 bp sequences"
	echo "	This could take a while..."
	findMotifs.pl "${TF}_${ZeroTag}_4.fa" fasta "${TF}_4_homer_denovo_long" \
	-fasta "$BASEDIR"/"${ZeroTag}"/"${ZeroTag}.fa" \
	-noredun -len 16,18 -noknown -p 4
	rm "${TF}_${ZeroTag}_4.fa"

	python "$CODEDIR"/htmltotext.py \
		"${TF}_4_homer_denovo_long"/homerResults.html \
		"${TF}_4_homer_denovo_long"/homerResults.txt 
		
	## Define spacer based on consensus sequence of top motif, generate COSMO motif files
	cd "$BASEDIR"/"$TF"
	python "$CODEDIR"/Consensus_sequence_search.py -c


	## Copy relevant motifs to TF directory
	cd "$BASEDIR"/"$TF"/Cycle4
	l_motif_number=$( head "${TF}_4_homer_denovo_long"/D_site_motif.txt )	
	cp "${TF}_4_homer_denovo_short/homerResults/motif1.motif" \
		"$BASEDIR"/"$TF"/"monomer.motif"
	cp "${TF}_4_homer_denovo_long/homerResults/${l_motif_number}" \
		"$BASEDIR"/"$TF"/"dimer.motif"
	rm "${TF}_4_homer_denovo_long"/D_site_motif.txt	

fi

#####################################
####### Enrichment analysis  ########
#####################################

if [ "$eflag" = "true" ]
then
	 
	## Load important modules
	mpurge
	module load homer/4.9-wrl
	module load python3
		
	## Convert zipped fastq files to fasta for zero cycle
	cd "$BASEDIR"/"$ZeroTag"
	zcat "$Fastq_bg" > "${ZeroTag}.fastq"
	paste - - - - < "${ZeroTag}.fastq" | cut -f 1,2 | sed 's/^@/>/' |\
	tr "\t" "\n" > "${ZeroTag}.fa"
	rm "${ZeroTag}.fastq"
	
	for Cycle in ${Rounds[@]}
	do

		## Convert fastq file to fasta file
		cd "$BASEDIR"/"$TF"/"Cycle${Cycle}"
		zcat "${Fastq_target[$Cycle]}" > "${TF}_${ZeroTag}_${Cycle}.fastq"
		paste - - - - < "${TF}_${ZeroTag}_${Cycle}.fastq" | cut -f 1,2 |\
		sed 's/^@/>/' |	 tr "\t" "\n" > "${TF}_${ZeroTag}_${Cycle}.fa"
		rm "${TF}_${ZeroTag}_${Cycle}.fastq"

		## Find prevalence of Cycle 4 short motifs
		echo "	Searching monomers in cycle $Cycle..."
		findMotifs.pl "${TF}_${ZeroTag}_${Cycle}.fa" fasta  \
		"${TF}_${Cycle}_monomer_homer" \
		-fasta "$BASEDIR"/"${ZeroTag}"/"${ZeroTag}.fa" -nomotif \
		-mknown "$BASEDIR"/"$TF"/"monomer.motif" -noweight -p 4 \
		-maskMotif "$BASEDIR"/"$TF"/"dimer.motif"
		
		## Convert knownResults.html to text file
		python "$CODEDIR"/htmltotext.py \
			"${TF}_${Cycle}_monomer_homer"/knownResults.html \
			"${TF}_${Cycle}_monomer_homer"/knownResults.txt 
		
		echo "	Searching dimers in cycle $Cycle..."
		findMotifs.pl "${TF}_${ZeroTag}_${Cycle}.fa" fasta  \
		"${TF}_${Cycle}_dimer_homer" \
		-fasta "$BASEDIR"/"${ZeroTag}"/"${ZeroTag}.fa" -nomotif \
		-mknown "$BASEDIR"/"$TF"/"dimer.motif" -noweight -p 4 \

		## Convert knownResults.html to text file
		python "$CODEDIR"/htmltotext.py \
			"${TF}_${Cycle}_dimer_homer"/knownResults.html \
			"${TF}_${Cycle}_dimer_homer"/knownResults.txt 


		## Remove fasta and fastq files - easy to redownload
		rm "${TF}_${ZeroTag}_${Cycle}.fa"
		# rm "${Fastq_target[$Cycle]}"
	done 

	rm "$BASEDIR"/"${ZeroTag}"/"${ZeroTag}.fa"
	cd "$BASEDIR"/"$TF"
	python "$CODEDIR"/Dimer_enrichment_calculator.py

fi 


#####################################
###########     COSMO     ###########
#####################################

if [ "$cflag" = "true" ]
then

	module purge
	module load python3

	## Set variables
	DIST=10; echo "	A max spacer distance is 10."
	if [ -z ${THRES+x} ]
	then
		echo "	A default of 0.8 being used for motif threshold."
		THRES=0.8
	else
		echo '	Motif threshold set to "$THRES". '
	fi

	if [ -z ${MODE+x} ]
	then
		echo "	Default COSMO mode 2 being used. COSMO will use the top two kmers from dimer site."
		MODE=2
	else
		echo '	Mode set to "$MODE".'
	fi
	cd "$BASEDIR"/"$TF" 

	if [ "$MODE" -eq 1 ]
	then
		python "$CODEDIR"/monomer_motif_trimmer.py --mon_length 4 --top

	else 
		if [ -d top_dimer_kmer_motifs ]
		then
			echo "	Motifs for COSMO have already been generated."
		else
			python "$CODEDIR"/Consensus_sequence_search.py -c		
		fi
		
	fi

	## Activate virtual environment that contains MOODs
	module purge
	cd "$CODEDIR"/../COSMO/cosmo
	source "$CODEDIR"/../COSMO/cosmo/venv/bin/activate
	
	cd "$BASEDIR"/"$ZeroTag"
	zcat "$Fastq_bg" > "${ZeroTag}.fastq"
	paste - - - - < "${ZeroTag}.fastq" | cut -f 1,2 | sed 's/^@/>/' |\
	tr "\t" "\n" | awk 'NR %2 {print ">chr" (NR+1)/2 ":1-30"} NR %2-1 {print $0}' \
		> "${ZeroTag}.fa"
	rm "${ZeroTag}.fastq"
	
	## Run COSMO on zero cycle
	Cycle=0
	"$CODEDIR"/../COSMO/cosmo/cosmo_v1.py --fasta "${ZeroTag}.fa" \
	--threshold $THRES --distance $DIST -p "$BASEDIR"/"$TF"/top_dimer_kmer_motifs/
	mkdir "$BASEDIR"/"$TF"/Cycle0
	mv cosmo.counts.tab "$BASEDIR"/"$TF"/Cycle0; cd "$BASEDIR"/"$TF"/Cycle0
	
	grep "motif1.jpwm|motif1.jpwm|FF" cosmo.counts.tab > "Cycle${Cycle}_motif1_motif1_FF.tab"
	grep "motif2.jpwm|motif2.jpwm|FF" cosmo.counts.tab > "Cycle${Cycle}_motif2_motif2_FF.tab"
	grep "motif1.jpwm|motif2.jpwm|FF" cosmo.counts.tab > "Cycle${Cycle}_motif1_motif2_FF.tab"
	grep "motif2.jpwm|motif1.jpwm|FF" cosmo.counts.tab > "Cycle${Cycle}_motif2_motif1_FF.tab"

	for Cycle in ${Rounds[@]}
	do 
		echo "	Starting cycle ${Cycle}"

		## Convert fastq file to fasta file
		cd "$BASEDIR"/"$TF"/"Cycle${Cycle}"
		zcat "${Fastq_target[$Cycle]}" > "${TF}_${ZeroTag}_${Cycle}.fastq"
		paste - - - - < "${TF}_${ZeroTag}_${Cycle}.fastq" | cut -f 1,2 |\
			sed 's/^@/>/' |	 tr "\t" "\n" | \
			awk 'NR %2 {print ">chr" (NR+1)/2 ":1-30"} NR %2-1 {print $0}' \
			> "${TF}_${ZeroTag}_${Cycle}.fa"	
		rm "${TF}_${ZeroTag}_${Cycle}.fastq"
		
		## Run foreground scan
		# "$CODEDIR"/../COSMO/cosmo/cosmo_v1.py --fasta \
			# "${TF}_${ZeroTag}_${Cycle}.fa" \
			# --threshold $THRES --distance $DIST -p ../top_dimer_kmer_motifs/
		
		## Run 30 background scans
		if [ "$Cycle" -eq 4 ]
		then
			echo "	Running background scans on Cycle 4 for statistics."
			mkdir 'COSMO_bg_scans'; cd 'COSMO_bg_scans'
			shuffle_count=($(seq 1 1 30))
			for shuff in ${shuffle_count[@]}
			do
				echo "	Running bg scan "$shuff" "
				# "$CODEDIR"/../COSMO/cosmo/cosmo_v1.py --fasta \
				# ../"${TF}_${ZeroTag}_${Cycle}.fa" --number $shuff --threshold $THRES \
				# --distance $DIST --pwmdir ../../top_dimer_kmer_motifs/ --scramflag
			done

			## Run stats
			mv ../cosmo.counts.tab .
			"$CODEDIR"/../COSMO/cosmo/cosmostats_v1.py -N 30 \
				> "${TF}_${Cycle}_stats.txt"
			mv "${TF}_${Cycle}_stats.txt" ..
			mv cosmo.counts.tab ..
			cd ..
		fi

		## Organize output - not all palindromic - orientation matters
		grep "motif1.jpwm|motif1.jpwm|FF" cosmo.counts.tab > "Cycle${Cycle}_motif1_motif1_FF.tab"
		grep "motif2.jpwm|motif2.jpwm|FF" cosmo.counts.tab > "Cycle${Cycle}_motif2_motif2_FF.tab"
		grep "motif1.jpwm|motif2.jpwm|FF" cosmo.counts.tab > "Cycle${Cycle}_motif1_motif2_FF.tab"
		grep "motif2.jpwm|motif1.jpwm|FF" cosmo.counts.tab > "Cycle${Cycle}_motif2_motif1_FF.tab"

		rm *.fa
	done	


	## Make heat maps - shows enrichment - no relevant statistics	
	deactivate
	cd "$BASEDIR"/"$TF"
	module load python3
	
	echo "	Generating heatmaps..."
	python "$CODEDIR"/heatmap_plotter.py 
	

	## Perform statistical analysis
	echo "	Running stats..."

	module load R
	if [ "$MODE" -eq 1 ]
	then
		Rscript "$CODEDIR"/chi-square_all_motif_arrangments.R

	else	
		Rscript "$CODEDIR"/chi-square.R
	fi
	
	mkdir COSMO_count_matrices;	mv *COSMO_counts* COSMO_count_matrices
fi 