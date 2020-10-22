#!/bin/bash 

# Author: cainu5
# Homer runs for SELEX data
# 10/09/20

(( TRACE )) && set -x
set -euo pipefail
trap 'Something went wrong. Script error on line #$LINENO."' ERR

module load homer/4.9-wrl
module load python3


## Bring in arguments
TF=${1:?Expected target file as argument #1}
Target_link=${2:?Expected target file as argument #2}
Zero_link=${3:?Expected target file as argument #4}
Target=$TF

cd ~/SELEX_analysis/testing

#####################################
## ACQUIRE FASTQ FROM ENA DATABASE ##
#####################################

## Check if folder already exists
if [ -d "$TF" ]
then
	echo "Have you run Homer for this dataset?"
	exit 1
fi

## Make necessary directories and download files
mkdir $TF; cd $TF
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

ZeroTag=$( echo $Zero_link | cut -d '_' -f 2 )
echo $ZeroTag; cd ~/SELEX_analysis/testing
if [ -e */"ZeroCycle_${ZeroTag}_0_0.fastq.gz" ]
then
	echo "Zero cycle has been downloaded"
else
	cd ~/SELEX_analysis/testing
	mkdir "$ZeroTag"; cd "$ZeroTag"
	curl -flOJ "$Zero_link"
fi

cd ~/SELEX_analysis/testing/"$Target"
cd Cycle1; curl -flOJ "$Target_link"
cd ../Cycle2; curl -flOJ "$LINK2"
cd ../Cycle3; curl -flOJ "$LINK3"
cd ../Cycle4; curl -flOJ "$LINK4"

# Grab downloaded fastq file names in a vector
Rounds=(1 2 3 4)
for Cycle in ${Rounds[@]}
do
	cd ~/SELEX_analysis/testing/"$Target"/"Cycle${Cycle}"
	Fastq_target["$Cycle"]=$( ls *"${Cycle}.fastq.gz" )
	echo "$Cycle" "${Fastq_target[@]}"
done

cd ~/SELEX_analysis/testing/"$ZeroTag"
Fastq_bg=$( ls *.fastq.gz )

## check file integrities
cd ~/SELEX_analysis/testing
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

if gzip -t "$ZeroTag/${Fastq_bg}"
then
	echo ""
else 
	echo "Download of Zero tag file was unsuccessful"
	exit 1
fi


# #####################################
# ########## Homer analysis ###########
# #####################################


## Convert zipped fastq files to fasta for zero cycle
cd ~/SELEX_analysis/testing/"$ZeroTag"
zcat "$Fastq_bg" > "${ZeroTag}.fastq"
paste - - - - < "${ZeroTag}.fastq" | cut -f 1,2 | sed 's/^@/>/' |\
tr "\t" "\n" > "${ZeroTag}.fa"
rm "${ZeroTag}.fastq"

## Convert zipped fastq files to fasta for cycle 2
cd ~/SELEX_analysis/testing/"$Target"/"Cycle2"
zcat "${Fastq_target[2]}" > "${Target}_${ZeroTag}_2.fastq"
paste - - - - < "${Target}_${ZeroTag}_2.fastq" | cut -f 1,2 |\
sed 's/^@/>/' |	 tr "\t" "\n" > "${Target}_${ZeroTag}_2.fa"
rm "${Target}_${ZeroTag}_2.fastq"

## De novo motif analysis of cycle 2 short
echo "Starting de novo motif analysis for Cycle2-short"
findMotifs.pl "${Target}_${ZeroTag}_2.fa" fasta "${Target}_2_homer_denovo_short" \
-fasta ~/SELEX_analysis/testing/"${ZeroTag}"/"${ZeroTag}.fa" \
-mcheck /data/weirauchlab/databank/appdata/HOMER/customizedTFs/v2.00/human_cisbp200.motif \
-noredun -len 6,8 -p 8 -noknown

python /data/weirauchlab/team/ches2d/MyTools/html2txt/code/html2txt.py  \
	"${Target}_2_homer_denovo_short"/homerResults.html \
	"${Target}_2_homer_denovo_short"/homerResults.txt 

## Pull out top relevant motif - searches for homeodomain output
grep -E -n \
	'POU2|ALX3|ALX4|ARX|BARHL2|BARX1|BSX|CART1|CDX1|CDX2|DLX1|DLX2|DLX3|DLX4|DLX5|DLX6|DMBX1|DPRX|DRGX|DUXA|EMX1|EMX2|EN1|EN2|ESX1|EVX1|EVX2|GBX1|GBX2|GSC|GSC2|GSX1|GSX2|HESX1|HMBOX1|HMX1|HMX2|HMX3|HNF1A|HNF1B|HOXA1|HOXA10|HOXA13|HOXA2|HOXB13HOXB2|HOXB3|HOXB5|HOXC10|HOXC11|HOXC12|HOXC13|HOXD11|HOXD12|HOXD13|HOXD8|IRX2|IRX5|ISL2|ISX|LBX2|LHX2|LHX6|LHX9|LMX1A|LMX1B|MEOX1|MEOX2|MIXL1|MNX1|MSX1|MSX2|NKX2-3|NKX2-8|NKX3-1|NKX3-2|NKX6-1|NKX6-2|NOTO|OTX1|OTX2|PDX1|PHOX2A|PHOX2B|PITX1|PITX3|PROP1|PRRX1|PRRX2|RAX|RAXL1|RHOXF1|SHOX|SHOX2|UNCX|VAX1|VAX2|VENTX|VSX1|VSX2' \
	"${Target}_2_homer_denovo_short"/homerResults.txt | head -1 \
	> "${Target}_2_homer_denovo_short"/top_short.txt
	
## Check if monomer was found
if [ -s "${Target}_2_homer_denovo_short"/top_short.txt ]
then
	echo 'TAAT motif found. Continuing to long de novo motif'
else
	echo 'No TAAT motif found. Exiting'
	rm "${Target}_${ZeroTag}_${Cycle}.fa"
	rm "${Fastq_target[$Cycle]}"
	rm ~/SELEX_analysis/testing/"${ZeroTag}"/"${ZeroTag}.fa"
	exit 1
fi

## De novo motif analysis of cycle 2 long
echo "Starting de novo motif analysis for Cycle2-long"
findMotifs.pl "${Target}_${ZeroTag}_2.fa" fasta "${Target}_2_homer_denovo_long" \
-fasta ~/SELEX_analysis/testing/"${ZeroTag}"/"${ZeroTag}.fa" \
-mcheck /data/weirauchlab/databank/appdata/HOMER/customizedTFs/v2.00/human_cisbp200.motif \
-noredun -len 16,18,20 -p 8 -noknown
rm "${Target}_${ZeroTag}_2.fa"

python /data/weirauchlab/team/ches2d/MyTools/html2txt/code/html2txt.py  \
	"${Target}_2_homer_denovo_long"/homerResults.html \
	"${Target}_2_homer_denovo_long"/homerResults.txt 

## Pull out top relevant long motif - searches for long homeodomain output
grep -E -n \
'M03076_2.00|M03104_2.00|M03107_2.00|M03109_2.00|M03116_2.00|M03127_2.00|M03144_2.00|M03150_2.00|M03155_2.00|M03159_2.00|M03162_2.00|M03170_2.00|M03184_2.00|M03188_2.00|M03197_2.00|M03246_2.00|M03249_2.00|M03255_2.00|M03260_2.00|M03288_2.00|M03300_2.00|M03302_2.00|M03304_2.00|M03310_2.00|M03311_2.00|M03790_2.00|M03792_2.00|M03807_2.00|M03819_2.00|M03821_2.00|M03828_2.00|M03864_2.00|M05120_2.00|M05121_2.00|M05233_2.00|M05234_2.00|M05235_2.00|M05236_2.00|M05315_2.00|M05316_2.00|M05349_2.00|M05356_2.00|M05444_2.00|M05477_2.00|M05482_2.00|M05483_2.00|M05491_2.00|M05501_2.00|M05503_2.00|M05504_2.00|M07976_2.00|M09139_2.00|M09144_2.00|M09151_2.00|M09168_2.00|M09176_2.00|M09179_2.00|M09191_2.00|M09200_2.00|M09219_2.00|M09220_2.00|M09222_2.00|M09227_2.00|M09228_2.00|M10652_2.00|M10665_2.00|M10690_2.00|M10691_2.00|M10693_2.00|M10704_2.00|M10717_2.00|M10741_2.00|M10759_2.00|M10767_2.00|M10771_2.00|M10784_2.00|M10788_2.00|M10823_2.00|M10826_2.00|M10827_2.00|M10829_2.00|M10835_2.00|M10843_2.00|M10848_2.00|M10849_2.00' \
"${Target}_2_homer_denovo_long"/homerResults.txt | head -1 \
> "${Target}_2_homer_denovo_long"/top_long.txt


## Check if dimer was found
if [ -s "${Target}_2_homer_denovo_long"/top_long.txt ]
then
	echo 'Dimer motif found. Continuing to searches.'
	## could also add feature to start de novo on cycle 3
else
	echo 'No dimer motif found. Exiting'
	rm "${Target}_${ZeroTag}_${Cycle}.fa"
	rm "${Fastq_target[$Cycle]}"
	rm ~/SELEX_analysis/testing/"${ZeroTag}"/"${ZeroTag}.fa"
	exit 1
fi


## Copy relevant motifs to TF directory
s_motif_number=$( cut -d : -f 1 "${Target}_2_homer_denovo_short"/top_short.txt )
l_motif_number=$( cut -d : -f 1 "${Target}_2_homer_denovo_long"/top_long.txt )
cp "${Target}_2_homer_denovo_short/homerResults/motif${s_motif_number}.motif" \
	~/SELEX_analysis/testing/"$Target"/"motif${s_motif_number}_short.motif"
cp "${Target}_2_homer_denovo_long/homerResults/motif${l_motif_number}.motif" \
	~/SELEX_analysis/testing/"$Target"/"motif${l_motif_number}_long.motif"
cat ~/SELEX_analysis/testing/"$Target"/"motif${s_motif_number}_short.motif" \
	~/SELEX_analysis/testing/"$Target"/"motif${l_motif_number}_long.motif" \
	> ~/SELEX_analysis/testing/"$Target"/"Cycle2.motif"
	

### Define spacer based on consensus sequence of top motif
# python Consensus_sequence_search.py \
# --target "${Target}"
# --l_motif "${l_motif_number}"

for Cycle in ${Rounds[@]}
do

	## Convert fastq file to fasta file
	cd ~/SELEX_analysis/testing/"$Target"/"Cycle${Cycle}"
	zcat "${Fastq_target[$Cycle]}" > "${Target}_${ZeroTag}_${Cycle}.fastq"
	paste - - - - < "${Target}_${ZeroTag}_${Cycle}.fastq" | cut -f 1,2 |\
	sed 's/^@/>/' |	 tr "\t" "\n" > "${Target}_${ZeroTag}_${Cycle}.fa"
	rm "${Target}_${ZeroTag}_${Cycle}.fastq"

	## Find prevalence of Cycle 2 short and long motifs
	echo "Beginning short motif homer run for cycle $Cycle"
	findMotifs.pl "${Target}_${ZeroTag}_${Cycle}.fa" fasta "${Target}_${Cycle}_homer" \
	-fasta ~/SELEX_analysis/testing/"${ZeroTag}"/"${ZeroTag}.fa" -nomotif \
	-mknown ~/SELEX_analysis/testing/"$Target"/"Cycle2.motif" -p 8

	## Convert knownResults.html to text file
	python /data/weirauchlab/team/ches2d/MyTools/html2txt/code/html2txt.py  \
		"${Target}_${Cycle}_homer"/knownResults.html \
		"${Target}_${Cycle}_homer"/knownResults.txt 


	## Remove fasta and fastq files - easy to redownload
	rm "${Target}_${ZeroTag}_${Cycle}.fa"
	rm "${Fastq_target[$Cycle]}"
done 

rm ~/SELEX_analysis/testing/"${ZeroTag}"/"${ZeroTag}.fa"

cd ~/SELEX_analysis/testing/"$Target"
python ~/SELEX_analysis/code/Dimer_enrichment_calculator.py

