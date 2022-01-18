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
	
while getopts "ceht:z:p:b:" opt; do
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
		exit 1
		;;
	# m)
		# MODE=${OPTARG}
		# if [ "$MODE" -eq 1 ] 
		# then
			# echo "	Running COSMO monomer to monomer."
		# elif [ "$MODE" == 2 ]
		# then
			# echo "	Running COSMO on top two 4mers from dimer site."
		# else
			# echo "	[-m]: Enter 1 for a monomer to monomer COSMO run. Enter 2 to 
		# run COSMO on the top two 4 mers from dimer site. Default is 2" &&
			# exit 1
		# fi
		# ;;
	z)
		Zero_link=${OPTARG}
		ZeroTag=$( echo ${Zero_link} | cut -d '_' -f 2 )
		;;
	p)
		PWM_file=${OPTARG}
		if [ ! -d "$PWM_file" ]; then
			echo "	[-p]: Enter a directory of motifs if you would like to use a motif generator other than HOMER."
			exit 1
		else
			echo "	$PWM_file will be used for analysis."			
		fi
		;;
	b)
		Fastq_bg=${OPTARG}
		if [ ! -e "$Fastq_bg" ]; then
			echo " [-b]: Enter a fastq.gz file from the initial library of HT-SELEX."
			exit 1
		fi
		;;
	?)			
		echo "	Usage to download ENA files: homer_runs.sh OPTIONS [ANALYSIS NAME] [Download link of Cycle 1 of ENA]" 
		echo "	Usage without downloading ENA files: homer_runs.sh OPTIONS [ANALYSIS NAME]"
		echo "		Program expects: "
		echo "			TF"
		echo "			├── Cycle1"
		echo "			│   └── [Cycle1].fastq.gz"
		echo "			├── Cycle2"
		echo "			│   └── [Cycle2].fastq.gz"
		echo "			├── Cycle3"
		echo "			│   └── [Cycle3].fastq.gz"
		echo "			├── Cycle4"
		echo "			    └── [Cycle4].fastq.gz"
		echo ""
		echo ""
		echo "	Program run options: "
		echo "	Default: all analyses run"
		echo "	[-h] HOMER de novo motif analysis only run."
		echo "	[-e] Enrichment analysis only run, requires Homer de novo motif files."
		echo "	[-c] COSMO only run, requires Homer de novo motif files."
		echo "	"
		echo "	Additional options:"
		echo "	[-z] Initial library cycle download link for normalization. If no download link"
		echo "		 or file [-b] provided, program uses cycle 1 library for normalization."
		echo "	[-b] Enter a fastq.gz file from the initial library of HT-SELEX. If no download link [-z]"
		echo "		 or file provided, program uses cycle 1 library for normalization"
		echo "	[-t] Motif threshold for MOODS. Must be a value between 0 and 1. "
		echo "		Default is 0.8."
		# echo "	[-m] COSMO motif run option. Enter 1 for a monomer to monomer COSMO run. Enter 2 to "
		# echo "		run COSMO on the top two 4 mers from dimer site. Default is 2.:"
		echo "	[-p]: Enter a motif file if you would like to use a motif generator "
		echo "		other than HOMER. This will overwrite -h option."
		echo ""		
		echo "	Program requirements:"
		echo "	1. HOMER"
		echo "	2. COSMO"
		echo "	3. Python 2.7" 
		echo "	4. Python 3 "
		echo ""; echo ""
		exit 1
		;;
	esac
done

shift $((OPTIND-1))

## Bring in arguments
TF=${1:?Expected TF as argument #1}
if [ ! -z "$2" ]; then Target_link=${2:?Expected ENA ftp download link for cycle 1 as argument #2}; fi

cd "$BASEDIR"

## Set flag defaults if no program options were selected
if [ "$cflag" != "true" ] && [ "$eflag" != "true" ] && [ "$hflag" != "true" ]
then
	echo "	De novo motif, enrichment, and COSMO analysis will be performed."
	cflag=true
	eflag=true
	hflag=true
fi

if [ ! -z ${PWM_file+x} ]
then
	echo "	Motif file given. Program will not run HOMER de novo."
	hflag=""
fi

## Check for de novo motif analyses if -h not given
if [ -z ${PWM_file+x} ] && [ ! "$hflag" = true ] && [ ! -e "$TF"/Cycle4/"$TF"_4_homer_denovo_long/homerResults.txt ]
then
	echo "	Cannot find a previous de novo motif analysis. Exiting. 
	Please select -h to run Homer de novo motif analysis or -p to enter your own motif file."
	exit 1
fi

#####################################
## ACQUIRE FASTQ FROM ENA DATABASE ##
#####################################

## Piece together target links

if [ ! -z ${Target_link+x} ]
then
	## Read sample accession of cycle 1 and extrapolate other cycle accessions
	SampAcc_Cycle1=$( echo "$Target_link" | cut -d / -f 5 | cut -d R -f 3 )
	SampAcc_Cycle2=$(( SampAcc_Cycle1 + 1 ))
	SampAcc_Cycle3=$(( SampAcc_Cycle1 + 2 ))
	SampAcc_Cycle4=$(( SampAcc_Cycle1 + 3 ))

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
	LINKS=("$LINK1" "$LINK2" "$LINK3" "$LINK4")
fi


if [ ! -d "$TF" ]
then
	mkdir "$TF"
	
fi

cd "$BASEDIR"/"$TF"

Rounds=(1 2 3 4)

for Cycle in ${Rounds[@]}
do
	if [ -e Cycle"$Cycle"/*.fastq.gz ]
	then
		echo "	FASTQ file for cycle $Cycle exist..."
				
	else
		echo "	Downloading files for cycle $Cycle from Jolma 2013..."
		
		## Make necessary directories and download files		
		mkdir Cycle"$Cycle"; cd Cycle"$Cycle"; wget -nv ${LINKS[$(( $Cycle - 1 )) ]}
		cd "$BASEDIR"/"$TF"
	fi
done


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
cd "$BASEDIR"; echo ""; Rounds=( 1 2 3 4 )

if [ ! -e */"ZeroCycle_${ZeroTag}_0_0.fastq.gz" ] && [ ! -z ${Zero_link+x} ]
then
	## Zero cycle not yet downloaded but link was given: download initial cycle	
	echo "	Zero cycle file link given. Downloading file..."
	mkdir "$ZeroTag"; cd "$ZeroTag"
	wget -nv "$Zero_link"
	cd "$BASEDIR"/"$ZeroTag"
	Fastq_bg=$( ls *.fastq.gz )
	ZeroTag=$(find $(pwd) | grep ${Fastq_bg} )
	if gzip --test "${Fastq_bg}"
	then
		echo "	${Fastq_bg} downloaded successfully."
	else 
		echo "	Download of Zero tag file was unsuccessful"
		cd ..
		rm -rf "$ZeroTag"
		exit 1
	fi

elif [ -z ${Zero_link+x} ] && [ ! -e "$Fastq_bg" ]
then
	## No existing file. Zero link was not given. Copy first cycle library to use as control
	echo "No zero cycle given"
	ZeroTag=$( echo ${Target_link} | tr '_', '\n' | grep .0N )
	Fastq_bg="${Fastq_target[1]}"
	Rounds=( 2 3 4 )
	mkdir "$ZeroTag"
	cp "$BASEDIR"/"$TF"/Cycle1/"${Fastq_target[1]}" "$ZeroTag"
	ZeroTag=$(find $(pwd) | grep ${Fastq_bg} )

elif [ -e */"ZeroCycle_${ZeroTag}_0_0.fastq.gz" ]
then
	## Zero cycle library already downloaded. Check file integrity
	echo "	Zero cycle given. Checking file integrity..."
	cd "$BASEDIR"/"$ZeroTag"
	Fastq_bg=$( find "$ZeroTag"/*.fastq.gz )
	ZeroTag=$(find $(pwd) | grep ${Fastq_bg} )
	if gzip --test "${Fastq_bg}"
	then
		echo "	${Fastq_bg} file intact."
	else 
		echo "	Download of Zero tag file was unsuccessful"
		cd ..
		rm -rf "$ZeroTag"
		exit 1
	fi

elif [ -e "$Fastq_bg" ]
then
	## If initial library was given
	echo "	Zero cycle given. Checking file integrity..."
	ZeroTag=$(find $(pwd) | grep ${Fastq_bg} )
	if gzip --test "${Fastq_bg}"
	then
		echo "	${Fastq_bg} file intact."
	else 
		echo "	Download of Zero tag file was unsuccessful"
		cd ..
		rm -rf "$ZeroTag"
		exit 1
	fi	

fi

Zero_dir=$(dirname "$ZeroTag")


## Prep PWM directory
if [ -d "$BASEDIR"/"$PWM_file" ]
then
	echo ""
	cd "$BASEDIR"/"$TF"/Cycle4; mkdir -p "${TF}_4_homer_denovo_long/homerResults/PWMs"
	ls "$BASEDIR"/"$PWM_file" | grep .motif | grep -v site | grep -v RV | while read dmotifs;
	do 
		cp "$BASEDIR"/"$PWM_file"/"$dmotifs" \
			"$BASEDIR"/"$TF"/Cycle4/"${TF}_4_homer_denovo_long/homerResults/PWMs/"
	done
		
fi	

# #####################################
# ###### HOMER de novo analysis  ######
# #####################################

if [ "$hflag" = "true" ]
then
	echo "
	
	
		"

	## Load important modules
	module load homer/4.9
	module load python3

	## Convert zipped fastq files to fasta for zero cycle
	cd "$BASEDIR"/"$ZeroTag"
	zcat "$Fastq_bg" > "${ZeroTag}.fastq"
	paste - - - - < "${ZeroTag}.fastq" | shuf | cut -f 1,2 | sed 's/^@/>/' |\
	tr "\t" "\n" | head -100000 > "${ZeroTag}.fa"
	rm "${ZeroTag}.fastq"

	## Convert zipped fastq files to fasta for cycle 4, extract 50k sequences (sampling)
	cd "$BASEDIR"/"$TF"/"Cycle4"
	zcat "${Fastq_target[4]}" > "${TF}_4.fastq" 
	paste - - - - < "${TF}_4.fastq" | shuf | cut -f 1,2 |\
	sed 's/^@/>/' | tr "\t" "\n" | head -100000 > "${TF}_4.fa"
	rm "${TF}_4.fastq"

	## De novo motif analysis of cycle 4 long
	echo "	Starting de novo motif analysis for Cycle4 for 16-18 bp sequences"
	echo "	This could take a while..."
	findMotifs.pl "${TF}_4.fa" fasta "${TF}_4_homer_denovo_long" \
	-fasta "$BASEDIR"/"${ZeroTag}"/"${ZeroTag}.fa" \
	-noredun -len 16,18 -noknown -p 4
	rm "${TF}_4.fa"

	python "$CODEDIR"/htmltotext.py \
		"${TF}_4_homer_denovo_long"/homerResults.html \
		"${TF}_4_homer_denovo_long"/homerResults.txt 
	
	mkdir "${TF}_4_homer_denovo_long/PWMs"; cd "${TF}_4_homer_denovo_long/PWMs"
	cp ../motif[0-9]..motif .
		
	## Define spacer based on consensus sequence of top motif, generate COSMO motif files
	cd "$BASEDIR"/"$TF"
	python "$CODEDIR"/Consensus_sequence_search.py -c

	cd Cycle4
	
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

# #####################################
# ####### Enrichment analysis  ########
# #####################################

if [ "$eflag" = "true" ]
then
	echo "	Starting enrichment analysis..."
	echo "


		"

	## Load important modules
	module purge
	module load homer/4.9-wrl
	module load python3
	
	## Define spacer based on consensus sequence of top motif, generate COSMO motif files
	cd "$BASEDIR"/"$TF"
	if [ ! -s long_motif_consensus.txt ]
	then
		if [ -d "$BASEDIR"/"$PWM_file" ]
		then
			python "$CODEDIR"/Consensus_sequence_search.py -cp
		else
			python "$CODEDIR"/Consensus_sequence_search.py -c
		fi
		
		count=0
		mkdir "PWMs_of_dimers"; cd Cycle4
		if [ -e "${TF}_4_homer_denovo_long"/D_site_motif.txt ]
		then
			cat "${TF}_4_homer_denovo_long"/D_site_motif.txt | while read line;
			do
				count=$(( count + 1 ))
				## Copy relevant motifs to TF directory
				cp "${TF}_4_homer_denovo_long/homerResults/PWMs/$line" \
					"$BASEDIR"/"$TF"/"PWMs_of_dimers"/"dimer_${line}"
				
				## Convert half sites to motif files		
				seq2profile.pl $( head -"$count" "$BASEDIR"/"$TF"/long_motif_consensus.txt |\
					tail -1 | head -c 4 ) > "$BASEDIR"/"$TF"/"PWMs_of_dimers"/site1_dimer_"$line"
				
				seq2profile.pl $( head -"$count"  "$BASEDIR"/"$TF"/long_motif_consensus.txt |\
					tail -1 | tail -c 5 | head -c 4 ) >\
					"$BASEDIR"/"$TF"/"PWMs_of_dimers"/site2_dimer_"$line"
					
			done
		# rm "${TF}_4_homer_denovo_long"/D_site_motif.txt
		fi
	fi

	## Stop run if no long consensus motifs were found
	cd "$BASEDIR"/"$TF"
	if [ $( head -1 long_motif_consensus.txt ) = 'N/A' ]
	then
		python "$CODEDIR"/Dimer_enrichment_calculator.py
		module load R
		Rscript "$CODEDIR"/chi-square.R
		exit 1
	fi
	 
		
	## Convert zipped fastq files to fasta for zero cycle
	cd "$Zero_dir"
	zcat "$Fastq_bg" > "${TF}_0.fastq"
	paste - - - - < "${TF}_0.fastq" | cut -f 1,2 | sed 's/^@/>/' |\
	tr "\t" "\n" > "${TF}_0.fa"
	rm "${TF}_0.fastq"
	
	
	## Search motifs at every cycle
	cd "$BASEDIR"/"$TF"	
	for Cycle in ${Rounds[@]}
	do
		echo "	Searching for monomers and dimers in cycle $Cycle ..."
		## Convert fastq file to fasta file
		cd "$BASEDIR"/"$TF"/"Cycle${Cycle}"
		zcat "${Fastq_target[$Cycle]}" > "${TF}_${Cycle}.fastq"
		paste - - - - < "${TF}_${Cycle}.fastq" | cut -f 1,2 |\
		sed 's/^@/>/' |	 tr "\t" "\n" > "${TF}_${Cycle}.fa"
		rm "${TF}_${Cycle}.fastq"

		## Search every dimer motif
		ls ../"PWMs_of_dimers" | grep .motif | grep -v site | while read dmotifs;
		do
			dmotif_number=$( echo $dmotifs | cut -d . -f 1  )
			
			## Find prevalence of half site 1
			findMotifs.pl "${TF}_${Cycle}.fa" fasta  \
			"${TF}_${Cycle}_site1_${dmotif_number}_mask_homer" \
			-fasta "$Zero_dir"/"${TF}_0.fa" -nomotif \
			-mknown "$BASEDIR"/"$TF"/"PWMs_of_dimers"/site1_"$dmotifs" -noweight  \
			-maskMotif "$BASEDIR"/"$TF"/"PWMs_of_dimers"/"$dmotifs"
			
			## Convert knownResults.html to text file
			python "$CODEDIR"/htmltotext.py \
				"${TF}_${Cycle}_site1_${dmotif_number}_mask_homer"/knownResults.html \
				"${TF}_${Cycle}_site1_${dmotif_number}_mask_homer"/knownResults.txt 
			

			## Find prevalence of half site 2
			findMotifs.pl "${TF}_${Cycle}.fa" fasta  \
			"${TF}_${Cycle}_site2_${dmotif_number}_mask_homer" \
			-fasta "$Zero_dir"/"${TF}_0.fa" -nomotif \
			-mknown "$BASEDIR"/"$TF"/"PWMs_of_dimers"/site2_"$dmotifs" -noweight  \
			-maskMotif "$BASEDIR"/"$TF"/"PWMs_of_dimers"/"$dmotifs"
			
			## Convert knownResults.html to text file
			python "$CODEDIR"/htmltotext.py \
				"${TF}_${Cycle}_site2_${dmotif_number}_mask_homer"/knownResults.html \
				"${TF}_${Cycle}_site2_${dmotif_number}_mask_homer"/knownResults.txt 
			
			## Find prevalence of Cycle 4 long motifs
			findMotifs.pl "${TF}_${Cycle}.fa" fasta  \
			"${TF}_${Cycle}_${dmotif_number}_homer" \
			-fasta "$Zero_dir"/"${TF}_0.fa" -nomotif \
			-mknown "$BASEDIR"/"$TF"/"PWMs_of_dimers"/"$dmotifs" -noweight


			## Convert knownResults.html to text file
			python "$CODEDIR"/htmltotext.py \
				"${TF}_${Cycle}_${dmotif_number}_homer"/knownResults.html \
				"${TF}_${Cycle}_${dmotif_number}_homer"/knownResults.txt 
				
		done
	done 

	rm "$Zero_dir"/"${TF}_0.fa"
	rm "$BASEDIR"/"$TF"/Cycle"${Cycle}"/"${TF}_${Cycle}.fa"
	cd "$BASEDIR"/"$TF"
	python "$CODEDIR"/Dimer_enrichment_calculator.py

fi 


# #####################################
# ###########     COSMO     ###########
# #####################################

if [ "$cflag" = "true" ]
then
	echo "
	
	
		"
	echo "	Starting COSMO analysis..."
	module purge
	module load python3
	module load R
	
	cd "$BASEDIR"/"$TF"
	
	if [ -e long_motif_consensus.txt ]
	then
		echo "	Motifs for COSMO have been generated."
	else
		if [ -d "$BASEDIR"/"$PWM_file" ]
		then
			python "$CODEDIR"/Consensus_sequence_search.py -cp
		else
			python "$CODEDIR"/Consensus_sequence_search.py -c
		fi
	
		mkdir "PWMs_of_dimers"; cd Cycle4
		if [ -e "${TF}_4_homer_denovo_long"/D_site_motif.txt ]
		then
			cat "${TF}_4_homer_denovo_long"/D_site_motif.txt | while read line;
			do
				## Copy relevant motifs to TF directory
				cp "${TF}_4_homer_denovo_long/homerResults/PWMs/$line" \
					"$BASEDIR"/"$TF"/"PWMs_of_dimers"/"dimer_${line}"
			done
		fi
	fi
	
	
	if [ $( head -1 long_motif_consensus.txt ) = 'N/A' ]
	then
		echo "	No dimer site found. Exiting"
		Rscript "$CODEDIR"/chi-square.R
		exit 1
	fi

	## Set variables
	DIST=10; echo "	The max spacer distance is set to 10."
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
			dimer=$( echo $dmotifs | cut -d _ -f 5-10 | cut -d . -f 1 )				
			cd "$BASEDIR"
			
			## Run COSMO on zero cycle if applicable
			if [ -e "$Zero_dir"/"$Fastq_bg" ] && [ "$Cycle" -eq 0 ]
			then		
				cd "$Zero_dir"
				zcat "$Fastq_bg" > "${TF}_0.fastq"
				paste - - - - < "${TF}_0.fastq" | cut -f 1,2 | sed 's/^@/>/' |\
					tr "\t" "\n" | awk 'NR %2 {print ">chr" (NR+1)/2 ":1-40"} NR %2-1 {print $0}' \
					> "${TF}_0.fa"
				rm "${TF}_0.fastq"
							
				"$CODEDIR"/../COSMO/cosmo/cosmo_v1.py --fasta "${TF}_0.fa" \
				--threshold "$THRES" --distance "$DIST" -p "$BASEDIR"/"$TF"/"$dmotifs"/
				
				if [ ! -d "$BASEDIR"/"$TF"/Cycle0 ]
				then
					mkdir "$BASEDIR"/"$TF"/Cycle0
				fi
				
				mkdir "$BASEDIR"/"$TF"/Cycle0/"$TF"_"$Cycle"_"$dimer"
				mv cosmo.counts.tab "$BASEDIR"/"$TF"/Cycle0/"$TF"_"$Cycle"_"$dimer"/
				cd "$BASEDIR"/"$TF"/Cycle0/"$TF"_"$Cycle"_"$dimer"/
			fi
			
			if [ "$Cycle" -gt 0 ]
			then
				## Convert fastq file to fasta file for non-zero cycles
				cd "$BASEDIR"/"$TF"/Cycle"$Cycle"
				zcat "${Fastq_target[$Cycle]}" > "${TF}_${Cycle}.fastq"
				paste - - - - < "${TF}_${Cycle}.fastq" | cut -f 1,2 |\
					sed 's/^@/>/' |	 tr "\t" "\n" | \
					awk 'NR %2 {print ">chr" (NR+1)/2 ":1-40"} NR %2-1 {print $0}' \
					> "${TF}_${Cycle}.fa"	
				rm "${TF}_${Cycle}.fastq"
				
				## Run foreground scan
				"$CODEDIR"/../COSMO/cosmo/cosmo_v1.py --fasta \
					"$BASEDIR"/"$TF"/Cycle"$Cycle"/"${TF}_${Cycle}.fa" \
					--threshold "$THRES" --distance "$DIST" -p ../"$dmotifs"/
					
				# Move results corresponding folder
				if [ -d "$TF"_"$Cycle"_"$dimer" ]
				then
					mv cosmo.counts.tab "$TF"_"$Cycle"_"$dimer"/ 
				else
					mkdir "$TF"_"$Cycle"_"$dimer"
					mv cosmo.counts.tab "$TF"_"$Cycle"_"$dimer"/
				fi
				
				# Run 30 background scans
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

				# Organize output - not all palindromic - orientation matters
				
			fi
			
			cd "$BASEDIR"/"$TF"/Cycle"$Cycle"/"$TF"_"$Cycle"_"$dimer"
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
	module load R

	echo ""
	echo "	Generating heatmaps..."
	cd "$BASEDIR"/"$TF"
	python "$CODEDIR"/heatmap_plotter.py 
	

	## Perform statistical analysis
	echo "	Running stats..."

	if [ "$MODE" -eq 1 ]
	then
		Rscript "$CODEDIR"/chi-square_all_motif_arrangments.R

	else	
		Rscript "$CODEDIR"/chi-square.R
	fi
	
	## Organize COSMO output
	ls | grep top_dimer_kmer_motifs | while read dmotifs
	do
		dimer=$( echo $dmotifs | cut -d _ -f 5-10 | cut -d . -f 1 )
		mv "$TF"_"$dimer"* top_dimer_kmer_motifs_"$dimer"/
		# mv site?_"$dmotifs".motif top_dimer_kmer_motifs_"$dmotifs"/
	done
fi 