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
Cycle=${2:?Expected target file as argument #2}
Target_link=${3:?Expected target file as argument #3}
Zero_link=${4:?Expected target file as argument #4}
Target=$TF

cd ~/NRLB/testing

## ACQUIRE FASTQ FROM ENA DATABASE ##



## Check if folder already exists
if [ -d "$TF" ]
then
	echo "Have you run Homer for this dataset?"
	exit 1
fi

mkdir $TF; cd $TF
mkdir Cycle1; mkdir Cycle2; mkdir Cycle3; mkdir Cycle4

# cd "$TF" || mkdir "$TF" && cd "$TF"
mkdir "Cycle${Cycle}"; cd "Cycle${Cycle}"
curl -flOJ "$Target_link"

cd ~/NRLB/testing/"$Target"/"Cycle${Cycle}"
Fastq_target=$( ls *"${Cycle}.fastq.gz" )
echo $Fastq_target

# # Check if zero cycle has been previously downloaded
ZeroTag=$( echo $Zero_link | cut -d '_' -f 2 )
echo $ZeroTag; cd ~/NRLB/testing
if [ -e */"ZeroCycle_${ZeroTag}_0_0.fastq.gz" ]
then
	echo "Zero cycle has been downloaded"
else
	cd ~/NRLB/testing
	mkdir "$ZeroTag"; cd "$ZeroTag"
	curl -flOJ "$Zero_link"
fi


cd ~/NRLB/testing/"$ZeroTag"
Fastq_bg=$( ls *.fastq.gz )
echo $Fastq_bg

cd ~/NRLB/testing
## check file integrities
if gzip -t "$TF"/"Cycle${Cycle}"/"${Fastq_target}"
then
	echo "Download of TF file was successful"
else
	echo "Download of TF file was unsuccessful"
	exit 1
fi

if gzip -t "$ZeroTag/${Fastq_bg}"
then
	echo "Download of Zero tag file was successful"
else 
	echo "Download of Zero tag file was unsuccessful"
	exit 1
fi


## Convert zipped fastq files to fasta
cd ~/NRLB/testing/"$ZeroTag"
zcat "$Fastq_bg" > "${ZeroTag}.fastq"
paste - - - - < "${ZeroTag}.fastq" | cut -f 1,2 | sed 's/^@/>/' |\
 tr "\t" "\n" > "${ZeroTag}.fa"
rm "${ZeroTag}.fastq"

cd ~/NRLB/testing/"$Target"/"Cycle${Cycle}"
zcat "$Fastq_target" > "${Target}_${ZeroTag}_${Cycle}.fastq"
paste - - - - < "${Target}_${ZeroTag}_${Cycle}.fastq" | cut -f 1,2 | sed 's/^@/>/' |\
 tr "\t" "\n" > "${Target}_${ZeroTag}_${Cycle}.fa"
rm "${Target}_${ZeroTag}_${Cycle}.fastq"

# Run homer for short
findMotifs.pl "${Target}_${ZeroTag}_${Cycle}.fa" fasta "${Target}_${Cycle}_homer_short" \
-fasta ~/NRLB/testing/"${ZeroTag}"/"${ZeroTag}.fa" \
-mcheck /data/weirauchlab/databank/appdata/HOMER/customizedTFs/v2.00/human_cisbp200.motif \
-noredun -len 6,8 -p 3 -noknown

	
# Calculate and output important info
python /data/weirauchlab/team/ches2d/MyTools/html2txt/code/html2txt.py  \
	"${Target}_${Cycle}_homer_short"/homerResults.html \
	"${Target}_${Cycle}_homer_short"/homerResults.txt 
	
grep -E \
	'ALX3|ALX4|ARX|BARHL2|BARX1|BSX|CART1|CDX1|CDX2|DLX1|DLX2|DLX3|DLX4|DLX5|DLX6|DMBX1|DPRX|DRGX|DUXA|EMX1|EMX2|EN1|EN2|ESX1|EVX1|EVX2|GBX1|GBX2|GSC|GSC2|GSX1|GSX2|HESX1|HMBOX1|HMX1|HMX2|HMX3|HNF1A|HNF1B|HOXA1|HOXA10|HOXA13|HOXA2|HOXB13HOXB2|HOXB3|HOXB5|HOXC10|HOXC11|HOXC12|HOXC13|HOXD11|HOXD12|HOXD13|HOXD8|IRX2|IRX5|ISL2|ISX|LBX2|LHX2|LHX6|LHX9|LMX1A|LMX1B|MEOX1|MEOX2|MIXL1|MNX1|MSX1|MSX2|NKX2-3|NKX2-8|NKX3-1|NKX3-2|NKX6-1|NKX6-2|NOTO|OTX1|OTX2|PDX1|PHOX2A|PHOX2B|PITX1|PITX3|PROP1|PRRX1|PRRX2|RAX|RAXL1|RHOXF1|SHOX|SHOX2|UNCX|VAX1|VAX2|VENTX|VSX1|VSX2' "${Target}_${Cycle}_homer_short"/homerResults.txt | head -1 |\
	> "${Target}_${Cycle}_homer_short"/top_short.txt ||\
	printf "%s\t%d\t%d\t%d" "No motif found" "0" "0" "1"  \
	> "${Target}_${Cycle}_homer_short"/top_short.txt

head -1 "${Target}_${Cycle}_homer_short"/homerResults.txt \
> "${Target}_${Cycle}_homer_short"/top_short.txt

Log10pvalue_short=$( cut -f -2 "${Target}_${Cycle}_homer_short"/top_short.txt |\
	cut -d % -f 1 )
Per_target_short=$( cut -f 3 "${Target}_${Cycle}_homer_short"/top_short.txt |\
	cut -d % -f 1)
Per_bg_short=$( cut -f 4 "${Target}_${Cycle}_homer_short"/top_short.txt |\
	cut -d % -f 1)
Fold_change_short=$( echo $Per_target_short $Per_bg_short | awk '{print $1/$2}' )
echo "%Target of Short Motif: $Per_target_short"
echo "%Bg of Short Motif: $Per_bg_short"
echo "Fold change of Short Motif: $Fold_change_short"

## Run homer for long
findMotifs.pl "${Target}_${ZeroTag}_${Cycle}.fa" fasta "${Target}_${Cycle}_homer_long" \
-fasta ~/NRLB/testing/"${ZeroTag}"/"${ZeroTag}.fa" \
-mcheck /data/weirauchlab/databank/appdata/HOMER/customizedTFs/v2.00/human_cisbp200.motif \
-noredun -len 16,18,20 -p 3 -noknown

python /data/weirauchlab/team/ches2d/MyTools/html2txt/code/html2txt.py  \
	"${Target}_${Cycle}_homer_long"/homerResults.html \
	"${Target}_${Cycle}_homer_long"/homerResults.txt 

grep -E 'M03184|M03255|M03127|M05483|M03828|M03155|M05482|M05491|M09200|M03790' \
	"${Target}_${Cycle}_homer_long"/homerResults.txt | head -1 \
	> "${Target}_${Cycle}_homer_long"/top_long.txt ||\
	head -1 "${Target}_${Cycle}_homer_long"/homerResults.txt \
	> "${Target}_${Cycle}_homer_long"/top_long.txt

Log10pvalue_long=$( cut -f -2 "${Target}_${Cycle}_homer_long"/top_long.txt |\
	cut -d % -f 1 )
Per_target_long=$( cut -f 3 "${Target}_${Cycle}_homer_long"/top_long.txt |\ 
	cut -d % -f 1 )
Per_bg_long=$( cut -f 4 "${Target}_${Cycle}_homer_long"/top_long.txt |\
	cut -d % -f 1 )
Fold_change_long=$( echo "scale=2;$Per_target_long/$Per_bg_long" | bc )
echo "%Target: $Per_target_long"
echo "%Bg: $Per_bg_long"
echo "Fold change: $Fold_change_long"

printf "%s\t%s\t%s\t%.2f%%\t%.2f%%\t%.2f\t%s\t%.2f%%\t%.2f%%\t%.2f\n" \
	"$Target" "$Cycle" \
	"$Log10pvalue_short" "$Per_target_short" "$Per_bg_short" "$Fold_change_short" \
	"$Log10pvalue_long" "$Per_target_long" "$Per_bg_long" "$Fold_change_long" \
	>> ~/Run_summary.txt

## Remove fasta and fastq files - easy to redownload
rm "${Target}_${ZeroTag}_${Cycle}.fa"
rm $Fastq_target
rm ~/NRLB/testing/"${ZeroTag}"/"${ZeroTag}.fa"