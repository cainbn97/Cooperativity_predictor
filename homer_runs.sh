#!/bin/bash 

# Author: cainu5
# SELEX analysis for cooperativity
# Main script that runs all other scripts
# 05/28/21

## Still need to add monomer COSMO run if no dimer is present ##

## Grab important file paths
BASEDIR="$(pwd)"
CODEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

## Implement a switch to allow COSMO and enrichment script to be run separately
cflag='';eflag='';hflag=''
	
while getopts "ceht:m:z:" opt; do
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
	z)
		Zero_link=${OPTARG}
		ZeroTag=$( echo ${OPTARG} | cut -d '_' -f 2 )
		;;
	?)
		echo 	
		"Script usage:
		
		Program run options: 
		Default: all analyses run
		[-h] HOMER de novo motif analysis only run.
		[-e] Enrichment analysis only run, requires Homer de novo motif files.
		[-c] COSMO only run, requires Homer de novo motif files.
		
		Additional options
		[-z] Initial library cycle download link for normalization. If not provided,
			uses cycle 1 library for normalization.
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

cd "$BASEDIR"

## Set flag defaults if no program options were selected
if [ "$cflag" != "true" ] && [ "$eflag" != "true" ] && [ "$hflag" != "true" ]
then
	echo "	De novo motif, enrichment, and COSMO analysis will be performed."
	cflag=true
	eflag=true
	hflag=true
fi

## Check for de novo motif analyses if -h not given
if [ "$hflag" != "true" ] && [ -e "$TF"/Cycle4/"$TF"_4_homer_denovo_long/homerResults.txt ]
then
	echo ""
else
	echo "	Cannot find a previous de novo motif analysis. Will run HOMER de novo"
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
	cd Cycle1; wget -nv "$Target_link"
	cd ../Cycle2; wget -nv "$LINK2"
	cd ../Cycle3; wget -nv "$LINK3"
	cd ../Cycle4; wget -nv "$LINK4"

fi

## Check file integrities and grab downloaded fastq file names in a vector
Rounds=(1 2 3 4)
for Cycle in ${Rounds[@]}
do
	cd "$BASEDIR"/"$TF"/"Cycle${Cycle}"
	Fastq_target["$Cycle"]=$( ls *.fastq.gz )
done

# check file integrities
cd "$BASEDIR"
for Cycle in ${Rounds[@]}
do
	if gzip --test "$TF"/"Cycle${Cycle}"/"${Fastq_target[$Cycle]}"
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

if [ ! -e */"ZeroCycle_${ZeroTag}_0_0.fastq.gz" ] && [ ! -z ${Zero_link+x} ]
then
	## Zero cycle not yet downloaded but link was given: download initial cycle	
	mkdir "$ZeroTag"; cd "$ZeroTag"
	wget -nv "$Zero_link"
	cd "$BASEDIR"/"$ZeroTag"
	Fastq_bg=$( ls *.fastq.gz )
	if gzip --test "${Fastq_bg}"
	then
		echo "	${Fastq_bg} downloaded successfully."
	else 
		echo "	Download of Zero tag file was unsuccessful"
		cd ..
		rm -rf "$ZeroTag"
		exit 1
	fi

elif [ -z ${Zero_link+x} ]
then
	## Zero link was not given. Copy first cycle library to use as control
	ZeroTag=$( echo ${Target_link} | tr '_', '\n' | grep .0N )
	Fastq_bg="${Fastq_target[1]}"
	mkdir "$ZeroTag"
	cp "$BASEDIR"/"$TF"/Cycle1/"${Fastq_target[1]}" "$ZeroTag"

elif [ -e */"ZeroCycle_${ZeroTag}_0_0.fastq.gz" ]
then
	## Zero cycle library already downloaded. Check file integrity
	cd "$BASEDIR"/"$ZeroTag"
	Fastq_bg=$( ls *.fastq.gz )
	if gzip --test "${Fastq_bg}"
	then
		echo "	${Fastq_bg} downloaded successfully."
	else 
		echo "	Download of Zero tag file was unsuccessful"
		cd ..
		rm -rf "$ZeroTag"
		exit 1
	fi
fi


#####################################
###### HOMER de novo analysis  ######
#####################################

if [ "$hflag" = "true" ]
then
	echo "
	
	
		"

	## Load important modules
	module load homer/4.9-wrl
	module load python3

	## Convert zipped fastq files to fasta for zero cycle
	cd "$BASEDIR"/"$ZeroTag"
	zcat "$Fastq_bg" > "${ZeroTag}.fastq"
	paste - - - - < "${ZeroTag}.fastq" | shuf | cut -f 1,2 | sed 's/^@/>/' |\
	tr "\t" "\n" | head -100000 > "${ZeroTag}.fa"
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

	cd Cycle4
	cp "${TF}_4_homer_denovo_short/homerResults/motif1.motif" \
		"$BASEDIR"/"$TF"/"monomer.motif"
	
	if [ -e "${TF}_4_homer_denovo_long"/D_site_motif.txt ]
	then
		cat "${TF}_4_homer_denovo_long"/D_site_motif.txt | while read line;
		do
			## Copy relevant motifs to TF directory
			l_motif_number=$( echo $line | cut -d . -f 1 | cut -d f -f 2 )
			cp "${TF}_4_homer_denovo_long/homerResults/$line" \
				"$BASEDIR"/"$TF"/dimer_"$l_motif_number".motif
			
			## Convert half sites to motif files		
			seq2profile.pl $( head -"$l_motif_number" "$BASEDIR"/"$TF"/long_motif_consensus.txt |\
				tail -1 | head -c 4 ) > "$BASEDIR"/"$TF"/site1_dimer_"$l_motif_number".motif
			
			seq2profile.pl $( head -"$l_motif_number" "$BASEDIR"/"$TF"/long_motif_consensus.txt |\
				tail -1 | tail -c 5 | head -c 4 ) >\
				"$BASEDIR"/"$TF"/site2_dimer_"$l_motif_number".motif
				
		done
		rm "${TF}_4_homer_denovo_long"/D_site_motif.txt
	fi
fi

#####################################
####### Enrichment analysis  ########
#####################################

if [ "$eflag" = "true" ]
then
	echo "	Starting enrichment analysis..."
	echo "


		"
		
	if [ -z ${Zero_link+x} ]
	then
		Rounds=( 2 3 4 )
	else
		Rounds=( 1 2 3 4 )
	fi

	## Stop run if no long consensus motifs were found
	cd "$BASEDIR"/"$TF"
	if [ $( head -1 long_motif_consensus.txt ) = 'N/A' ]
	then
		python "$CODEDIR"/Dimer_enrichment_calculator.py
		exit 1
	fi
	 
	## Load important modules
	module purge
	module load homer/4.9-wrl
	module load python3
		
	## Convert zipped fastq files to fasta for zero cycle
	cd "$BASEDIR"/"$ZeroTag"
	zcat "$Fastq_bg" > "${ZeroTag}.fastq"
	paste - - - - < "${ZeroTag}.fastq" | cut -f 1,2 | sed 's/^@/>/' |\
	tr "\t" "\n" > "${ZeroTag}.fa"
	rm "${ZeroTag}.fastq"
	
	
	## Search motifs at every cycle
	cd "$BASEDIR"/"$TF"	
	for Cycle in ${Rounds[@]}
	do
		echo "	Searching for monomers and dimers in cycle $Cycle ..."
		## Convert fastq file to fasta file
		cd "$BASEDIR"/"$TF"/"Cycle${Cycle}"
		zcat "${Fastq_target[$Cycle]}" > "${TF}_${ZeroTag}_${Cycle}.fastq"
		paste - - - - < "${TF}_${ZeroTag}_${Cycle}.fastq" | cut -f 1,2 |\
		sed 's/^@/>/' |	 tr "\t" "\n" > "${TF}_${ZeroTag}_${Cycle}.fa"
		rm "${TF}_${ZeroTag}_${Cycle}.fastq"

		## Search every dimer motif
		ls .. | grep dimer_..motif | grep -v site | while read dmotifs;
		do
			dmotif_number=$( echo $dmotifs | cut -d . -f 1 | cut -d f -f 2 )
			
			## Find prevalence of half site 1
			findMotifs.pl "${TF}_${ZeroTag}_${Cycle}.fa" fasta  \
			"${TF}_${Cycle}_site1_${dmotif_number}_mask_homer" \
			-fasta "$BASEDIR"/"${ZeroTag}"/"${ZeroTag}.fa" -nomotif \
			-mknown "$BASEDIR"/"$TF"/site1_"$dmotifs" -noweight -p 4 \
			-maskMotif "$BASEDIR"/"$TF"/"$dmotifs"
			
			## Convert knownResults.html to text file
			python "$CODEDIR"/htmltotext.py \
				"${TF}_${Cycle}_site1_${dmotif_number}_mask_homer"/knownResults.html \
				"${TF}_${Cycle}_site1_${dmotif_number}_mask_homer"/knownResults.txt 
						
			## Find prevalence of half site 2
			findMotifs.pl "${TF}_${ZeroTag}_${Cycle}.fa" fasta  \
			"${TF}_${Cycle}_site2_${dmotif_number}_mask_homer" \
			-fasta "$BASEDIR"/"${ZeroTag}"/"${ZeroTag}.fa" -nomotif \
			-mknown "$BASEDIR"/"$TF"/site2_"$dmotifs" -noweight -p 4 \
			-maskMotif "$BASEDIR"/"$TF"/"$dmotifs"
			
			## Convert knownResults.html to text file
			python "$CODEDIR"/htmltotext.py \
				"${TF}_${Cycle}_site2_${dmotif_number}_mask_homer"/knownResults.html \
				"${TF}_${Cycle}_site2_${dmotif_number}_mask_homer"/knownResults.txt 
			
			## Find prevalence of Cycle 4 long motifs
			findMotifs.pl "${TF}_${ZeroTag}_${Cycle}.fa" fasta  \
			"${TF}_${Cycle}_${dmotif_number}_homer" \
			-fasta "$BASEDIR"/"${ZeroTag}"/"${ZeroTag}.fa" -nomotif \
			-mknown "$BASEDIR"/"$TF"/"$dmotifs" -noweight -p 4 \

			## Convert knownResults.html to text file
			python "$CODEDIR"/htmltotext.py \
				"${TF}_${Cycle}_${dmotif_number}_homer"/knownResults.html \
				"${TF}_${Cycle}_${dmotif_number}_homer"/knownResults.txt 
		done
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
	echo "
	
	
		"
	echo "	Starting COSMO analysis..."
	module purge
	module load python3
	
	cd "$BASEDIR"/"$TF"
	
	if [ $( head -1 long_motif_consensus.txt ) = 'N/A' ]
	then
		echo "	No dimer site found. Exiting"
		exit 1
	fi

	## Set variables
	DIST=20; echo "	The max spacer distance is set to 10."
	if [ -z ${THRES+x} ]
	then
		echo "	A default of 0.8 being used for motif threshold."
		THRES=0.8
	else
		echo "	Motif threshold set to "$THRES". "
	fi

	if [ -z ${MODE+x} ]
	then
		echo "	Default COSMO mode 2 being used. COSMO will use the top two kmers from dimer site."
		MODE=2
	else
		echo "	Mode set to "$MODE"."
	fi
	
	cd "$BASEDIR"/"$TF" 
	if [ "$MODE" -eq 1 ]
	then
		python "$CODEDIR"/monomer_motif_trimmer.py --mon_length 4 --top

	else 
		## Check if the COSMO motifs have been generated
		if [ -d top_dimer_kmer_motifs_dimer_1 ]
		then
			echo "	Motifs for COSMO have already been generated."
		else
			python "$CODEDIR"/Consensus_sequence_search.py -c		
		fi
		
	fi
	
	echo ""
	
	## Activate virtual environment that contains MOODs
	module purge
	cd "$CODEDIR"/../COSMO/cosmo
	source "$CODEDIR"/../COSMO/cosmo/venv/bin/activate
	
	## Run COSMO on all cycles
	Rounds=(0 1 2 3 4)
	for Cycle in ${Rounds[@]}
	do
		echo "	Running COSMO on cycle $Cycle..."
		
		## Run COSMO on all dimer site motifs		
		ls "$BASEDIR"/"$TF" | grep top_dimer_kmer_motifs | while read dmotifs
		do 
			dimer=$( echo $dmotifs | cut -d _ -f 5,6 )				
			cd "$BASEDIR"
			
			## Run COSMO on zero cycle if applicable
			if [ -e "$ZeroTag"/"ZeroCycle_${ZeroTag}_0_0.fastq.gz" ] && [ "$Cycle" -eq 0 ]
			then		
				cd "$BASEDIR"/"$ZeroTag"
				zcat "$Fastq_bg" > "${ZeroTag}.fastq"
				paste - - - - < "${ZeroTag}.fastq" | cut -f 1,2 | sed 's/^@/>/' |\
					tr "\t" "\n" | awk 'NR %2 {print ">chr" (NR+1)/2 ":1-40"} NR %2-1 {print $0}' \
					> "${ZeroTag}.fa"
				rm "${ZeroTag}.fastq"
							
				"$CODEDIR"/../COSMO/cosmo/cosmo_v1.py --fasta "${ZeroTag}.fa" \
				--threshold "$THRES" --distance "$DIST" -p "$BASEDIR"/"$TF"/"$dmotifs"/
				
				if [ ! -d "$BASEDIR"/"$TF"/Cycle0 ]
				then
					mkdir "$BASEDIR"/"$TF"/Cycle0
				fi
				
				mkdir "$BASEDIR"/"$TF"/Cycle0/"$TF"_"$Cycle"_"$dimer"_homer
				mv cosmo.counts.tab "$BASEDIR"/"$TF"/Cycle0/"$TF"_"$Cycle"_"$dimer"_homer/
				cd "$BASEDIR"/"$TF"/Cycle0/"$TF"_"$Cycle"_"$dimer"_homer/
			fi
			
			if [ "$Cycle" -gt 0 ]
			then
				## Convert fastq file to fasta file for non-zero cycles
				cd "$BASEDIR"/"$TF"/Cycle"$Cycle"
				zcat "${Fastq_target[$Cycle]}" > "${TF}_${ZeroTag}_${Cycle}.fastq"
				paste - - - - < "${TF}_${ZeroTag}_${Cycle}.fastq" | cut -f 1,2 |\
					sed 's/^@/>/' |	 tr "\t" "\n" | \
					awk 'NR %2 {print ">chr" (NR+1)/2 ":1-40"} NR %2-1 {print $0}' \
					> "${TF}_${ZeroTag}_${Cycle}.fa"	
				rm "${TF}_${ZeroTag}_${Cycle}.fastq"
				
				## Run foreground scan
				"$CODEDIR"/../COSMO/cosmo/cosmo_v1.py --fasta \
					"$BASEDIR"/"$TF"/Cycle"$Cycle"/"${TF}_${ZeroTag}_${Cycle}.fa" \
					--threshold "$THRES" --distance "$DIST" -p ../"$dmotifs"/
					
				# Move results corresponding folder
				if [ -d "$TF"_"$Cycle"_"$dimer"_homer ]
				then
					mv cosmo.counts.tab "$TF"_"$Cycle"_"$dimer"_homer/ 
				else
					mkdir "$TF"_"$Cycle"_"$dimer"_homer
					mv cosmo.counts.tab "$TF"_"$Cycle"_"$dimer"_homer/
				fi
				
				## Run 30 background scans
				# if [ "$Cycle" -eq 4 ]
				# then
					# echo "	Running background scans on Cycle 4 for statistics."
					# mkdir 'COSMO_bg_scans'; cd 'COSMO_bg_scans'
					# shuffle_count=($(seq 1 1 30))
					# for shuff in ${shuffle_count[@]}
					# do
						# echo "	Running bg scan $shuff"
						# "$CODEDIR"/../COSMO/cosmo/cosmo_v1.py --fasta \
						# ../"${TF}_${ZeroTag}_${Cycle}.fa" --number "$shuff" --threshold "$THRES" \
						# --distance "$DIST" --pwmdir ../../"$dmotifs"/ --scramflag
					# done

					# ## Run stats
					# mv ../cosmo.counts.tab .
					# "$CODEDIR"/../COSMO/cosmo/cosmostats_v1.py -N 30 \
						# > "${TF}_${Cycle}_stats.txt"
					# mv "${TF}_${Cycle}_stats.txt" ..
					# mv cosmo.counts.tab ..
					# cd ..
				# fi

				## Organize output - not all palindromic - orientation matters
				
			fi
			
			cd "$BASEDIR"/"$TF"/Cycle"$Cycle"/"$TF"_"$Cycle"_"$dimer"_homer
			grep "motif1.jpwm|motif1.jpwm|FF" cosmo.counts.tab > "Cycle${Cycle}_motif1_motif1_FF.tab"
			grep "motif2.jpwm|motif2.jpwm|FF" cosmo.counts.tab > "Cycle${Cycle}_motif2_motif2_FF.tab"
			grep "motif1.jpwm|motif2.jpwm|FF" cosmo.counts.tab > "Cycle${Cycle}_motif1_motif2_FF.tab"
			grep "motif2.jpwm|motif1.jpwm|FF" cosmo.counts.tab > "Cycle${Cycle}_motif2_motif1_FF.tab"
			
		done
		
		## Remove intermediate files
		cd "$BASEDIR"/"$TF"/Cycle"$Cycle"
		rm -f *.fa
		
	done	

	## Make heat maps - shows enrichment - no relevant statistics - 
	## run all at once to avoid deactivating and inactivating COSMO a bunch
	deactivate
	module load python3

	echo ""
	echo "	Generating heatmaps..."
	cd "$BASEDIR"/"$TF"
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
	
	## Organize COSMO output
	ls | grep top_dimer_kmer_motifs | cut -d _ -f 5,6 | while read dmotifs
	do
		mv "$TF"_"$dmotifs"* top_dimer_kmer_motifs_"$dmotifs"/
	done
fi 