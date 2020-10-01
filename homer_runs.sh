#!/bin/bash 

# Author: cainu5
# Homer runs for SELEX data
# 09/29/20

(( TRACE )) && set -x
set -euo pipefail
trap 'Something went wrong. Script error on line #$LINENO."' ERR

module load homer/4.9-wrl

Target=${1:?Expected target file as argument #1}
Cycle=${2:?Expected target file as argument #1}

cd ~/NRLB/testing/"$Target"/"Cycle${Cycle}"
Fastq_target=$( ls *"${Cycle}.fastq.gz" )
echo $Fastq_target

## Find corresponding zero cycle
ZeroTag=$( echo "${Fastq_target}" | cut -d '_' -f 2 )
echo $ZeroTag

cd ~/NRLB/testing/"$ZeroTag"
Fastq_bg=$( ls *.fastq.gz )
echo $Fastq_bg

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

## Run homer for short
findMotifs.pl "${Target}_${ZeroTag}_${Cycle}.fa" fasta "${Target}_homer_short" \
-fasta ~/NRLB/testing/"${ZeroTag}"/"${ZeroTag}.fa" \
-mknown /data/weirauchlab/databank/appdata/HOMER/customizedTFs/v2.00/human_cisbp200.motif \
-mcheck /data/weirauchlab/databank/appdata/HOMER/customizedTFs/v2.00/human_cisbp200.motif \
-noredun -len 8,10 -p 3

## Run homer for long
findMotifs.pl "${Target}_${ZeroTag}_${Cycle}.fa" fasta "${Target}_homer_long" \
-fasta ~/NRLB/testing/"${ZeroTag}"/"${ZeroTag}.fa" \
-mcheck /data/weirauchlab/databank/appdata/HOMER/customizedTFs/v2.00/human_cisbp200.motif \
-noredun -len 16,18 -p 3

rm "${Target}_${ZeroTag}_${Cycle}.fa"
rm ~/NRLB/testing/"${ZeroTag}"/"${ZeroTag}.fa"