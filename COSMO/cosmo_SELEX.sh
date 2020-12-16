#!/bin/bash 

# Author: cainu5
# COSMO runs for SELEX data
# 12/10/20

(( TRACE )) && set -x
set -euo pipefail
trap 'Something went wrong. Script error on line #$LINENO."' ERR

## Bring in arguments for file download
TF=${1:?Expected target file as argument #1}
Target_link=${2:?Expected target file as argument #2}
Zero_link=${3:?Expected target file as argument #3}
Target=$TF

## Bring in arguments for COSMO
THRES=${4:?Expected motif threshold as argument #4}
DIST=${5:?Expected max spacer length as argument #5}
PWM_MODE=${6:?Expected PWM mode ( 1 == monomer motif, 2 == dimer motif) as argument #6}

cd ~/SELEX_analysis/COSMO_output

#####################################
## ACQUIRE FASTQ FROM ENA DATABASE ##
#####################################


## Check if folder already exists
# if [ -d "$TF" ]
# then
	# echo "Have you run COSMO for this dataset?"
	# exit 1
# fi

## Make necessary directories and download files
#mkdir $TF; 
cd $TF
#mkdir Cycle0; mkdir Cycle1; mkdir Cycle2; mkdir Cycle3; mkdir Cycle4

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

ZeroTag=$( echo $Zero_link | cut -d '_' -f 2 )

cd ~/SELEX_analysis/COSMO_output/"$Target"
# cd Cycle0; curl -flOJ "$Zero_link"
# cd ../Cycle1; curl -flOJ "$Target_link"
# cd ../Cycle2; curl -flOJ "$LINK2"
# cd ../Cycle3; curl -flOJ "$LINK3"
# cd ../Cycle4; curl -flOJ "$LINK4"

## Grab downloaded fastq file names in a vector
Rounds=(0 1 2 3 4)
for Cycle in ${Rounds[@]}
do
	cd ~/SELEX_analysis/COSMO_output/"$Target"/"Cycle${Cycle}"
	Fastq_target["$Cycle"]=$( ls *"${Cycle}.fastq.gz" )
	echo "$Cycle" "${Fastq_target[@]}"
done

## check file integrities
cd ~/SELEX_analysis/COSMO_output
for Cycle in ${Rounds[@]}
do
	if gzip -t "$TF"/"Cycle${Cycle}"/"${Fastq_target[$Cycle]}"
	then
		echo ""
	else
		echo 'Download of ${Fastq_target["$Cycle"]} was unsuccessful'
		exit 1
	fi
done



#####################################
##       Determine matrices        ##
#####################################

module load python3
#mkdir "$TF"/motifs
cd ~/SELEX_analysis/testing/"$Target"

if [ "$PWM_MODE" -eq 2 ]
then

	python ~/SELEX_analysis/code/Consensus_sequence_search.py 

elif [ "$PWM_MODE" -eq 1 ]
then
	
	python ~/SELEX_analysis/code/COSMO/monomer_motif_trimmer.py --mon_length 6 \
		--top

fi
module purge

#####################################
###########     COSMO     ###########
#####################################

## Activate virtual environment that contains MOODs
cd ~/SELEX_analysis/COSMO/cosmo
source ~/SELEX_analysis/COSMO/cosmo/venv/bin/activate

for Cycle in ${Rounds[@]}
do 
	echo "Starting cycle ${Cycle}"

	## Convert fastq file to fasta file
	cd ~/SELEX_analysis/COSMO_output/"$Target"/"Cycle${Cycle}"
	zcat "${Fastq_target[$Cycle]}" > "${Target}_${ZeroTag}_${Cycle}.fastq"
	paste - - - - < "${Target}_${ZeroTag}_${Cycle}.fastq" | cut -f 1,2 |\
		sed 's/^@/>/' |	 tr "\t" "\n" | \
		awk 'NR %2 {print ">chr" (NR+1)/2 ":1-20"} NR %2-1 {print $0}' \
		> "${Target}_${ZeroTag}_${Cycle}.fa"	
	rm "${Target}_${ZeroTag}_${Cycle}.fastq"

	## Run foreground scan
	~/SELEX_analysis/COSMO/cosmo/cosmo_v1.py --fasta \
		"${Target}_${ZeroTag}_${Cycle}.fa" \
		--threshold $THRES --distance $DIST -p ../motifs/
	
	## Run 30 background scans
	shuffle_count=($(seq 1 1 30))
	for shuff in ${shuffle_count[@]}
	do
		~/SELEX_analysis/COSMO/cosmo/cosmo_v1.py --fasta \
		"${Target}_${ZeroTag}_${Cycle}.fa" --number $shuff --threshold $THRES \
		--distance $DIST --pwmdir ../motifs/ --scramflag
	done

	## Run stats
	~/SELEX_analysis/COSMO/cosmo/cosmostats_v1.py -N 30 \
		> "${Target}_${Cycle}_stats.txt"

	## Organize output
	grep "motif1.jpwm" cosmo.counts.tab | grep -v "motif2.jpwm" | grep "FF" \
		> "Cycle${Cycle}_motif1_motif1_FF.tab"
	grep "motif1.jpwm" cosmo.counts.tab | grep "motif2.jpwm" | grep "FF" \
		> "Cycle${Cycle}_motif1_motif2_FF.tab"
	grep "motif2.jpwm" cosmo.counts.tab | grep -v "motif1.jpwm" | grep "FF" \
		> "Cycle${Cycle}_motif2_motif2_FF.tab"

done	


## Make heat maps - shows enrichment - no relevant statistics	
deactivate
cd ..
module load python3
python ~/SELEX_analysis/code/COSMO/heatmap_plotter.py 

#rm */*.fastq.gz
rm */*.fa



