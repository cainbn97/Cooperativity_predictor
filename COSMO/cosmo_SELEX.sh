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
# mkdir $TF; 
# cd $TF
# mkdir Cycle1; mkdir Cycle2; mkdir Cycle3; mkdir Cycle4

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
# curl -flOJ "$Zero_link"
# cd Cycle1; curl -flOJ "$Target_link"
# cd ../Cycle2; curl -flOJ "$LINK2"
# cd ../Cycle3; curl -flOJ "$LINK3"
# cd ../Cycle4; curl -flOJ "$LINK4"

## Grab downloaded fastq file names in a vector
Rounds=(1 2 3 4)
for Cycle in ${Rounds[@]}
do
	cd ~/SELEX_analysis/COSMO_output/"$Target"/"Cycle${Cycle}"
	Fastq_target["$Cycle"]=$( ls *"${Cycle}.fastq.gz" )
	echo "$Cycle" "${Fastq_target[@]}"
done

cd ~/SELEX_analysis/COSMO_output/"$TF"
Fastq_bg=$( ls *0_0.fastq.gz )

# check file integrities
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

if gzip -t "$TF/${Fastq_bg}"
then
	echo ""
else 
	echo "Download of Zero tag file was unsuccessful"
	exit 1
fi


#####################################
##       Determine matrices        ##
#####################################

module load python3

l_motif_number=$( cut -f 1 ~/SELEX_analysis/testing/"$Target"/"Cycle4"/"${Target}_4_homer_denovo_long"/top_long.txt )
#mkdir "$TF"/motifs

cd ~/SELEX_analysis/testing/"$Target"
python ~/SELEX_analysis/code/Consensus_sequence_search.py 

module purge

#####################################
###########     COSMO     ###########
#####################################

## Activate virtual environment that contains MOODs
cd ~/SELEX_analysis/COSMO/cosmo
source ~/SELEX_analysis/COSMO/cosmo/venv/bin/activate

## Convert fastq file to fasta file - Zero Cycle
# echo "$Fastq_bg"
# cd ~/SELEX_analysis/COSMO_output/"$Target"
# zcat "$Fastq_bg" > "ZeroCycle.fastq"
# paste - - - - < "ZeroCycle.fastq" | cut -f 1,2 |\
	# sed 's/^@/>/' |	 tr "\t" "\n" |\
	# awk 'NR %2 {print ">chr" (NR+1)/2 ":1-20"} NR %2-1 {print $0}'> "ZeroCycle.fa"	
# rm "ZeroCycle.fastq"
	
# ## Run "background scan" - do not use scrambling option in COSMO
# ~/SELEX_analysis/COSMO/cosmo/cosmo_v1.py -fa "ZeroCycle.fa" \
	# -t $THRES -d $DIST -p motifs/
# mv cosmo.counts.tab cosmo.counts.tab.1

for Cycle in ${Rounds[@]}
do 

	## Convert fastq file to fasta file
	cd ~/SELEX_analysis/COSMO_output/"$Target"/"Cycle${Cycle}"
	zcat "${Fastq_target[$Cycle]}" > "${Target}_${ZeroTag}_${Cycle}.fastq"
	paste - - - - < "${Target}_${ZeroTag}_${Cycle}.fastq" | cut -f 1,2 |\
		sed 's/^@/>/' |	 tr "\t" "\n" | \
		awk 'NR %2 {print ">chr" (NR+1)/2 ":1-20"} NR %2-1 {print $0}' \
		> "${Target}_${ZeroTag}_${Cycle}.fa"	
	rm "${Target}_${ZeroTag}_${Cycle}.fastq"

	## Run foreground scan
	~/SELEX_analysis/COSMO/cosmo/cosmo_v1.py -fa "${Target}_${ZeroTag}_${Cycle}.fa" \
		-t $THRES -d $DIST -p ../motifs/
	
	## Move zero cycle and perform stats
	cp ../cosmo.counts.tab.1 cosmo.counts.tab.1
	~/SELEX_analysis/COSMO/cosmo/cosmostats_v1.py -N 1 > "${Target}_${Cycle}_stats.txt"

done	
	
	

