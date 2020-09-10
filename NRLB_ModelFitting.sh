#!/bin/bash

## Author: cainu5
## NRLB analysis 
## 09/01/20

(( TRACE )) && set -x
set -euo pipefail
trap 'Something went wrong. Script error on line #$LINENO.'  ERR


## Bring in arguments
TF=${1:?Expected a TF as argument #1}
LINK1=${2:?Expected download link of 1st cycle as argument #2}
LINK2=${3:?Expected download link of zero cycle as argument #3}

cd NRLBruns

## ACQUIRE DATA FROM ENA DATABASE ##

## Check if folder already exists
if [ -d "$TF" ]
then
	echo "Have you run NRLB for this dataset?"
	exit 1
fi

mkdir "$TF"


## Check if zero cycle has been previously downloaded
ZeroTag=$( echo $LINK2 | cut -d '_' -f 2 )
if [-d "$ZeroTag" ]
then
	#If zero tag was downloaded, only download TF file
	echo "Zero cycle has been downloadeed"
	cd "$TF"
	curl -fLOJ "$LINK1"
else
	mkdir "$ZeroTag"
	#download files into TF folder and then move Zero cycle; so parallel can 
	# be used
	cd "$TF"
	parallel curl -fl0J ::: "$LINK1" "$LINK2"
	mv *_ZeroCycle*.fastq.gz ../"$ZeroTag"
	
fi



cd ..
#check file integrities
if gzip -t "$TF/$TF*.fastq.gz"
then
	echo "Download of TF file was successful"
else
	echo "Download of TF file was unsuccessful"
	exit 1
fi

if gzip -t "$ZeroTag/*$ZeroTag*.fastq.gz"
then
	echo "Download of Zero tag file was successful"
else 
	echo "Download of Zero tag file was unsuccessful"
	exit 1
fi


## RUN NRLB MODEL FITTING (JAVA) ##

module load jdk
module load R/3.5.3-wrl

cd $TF

bash ~/NRLB/NRLB/ModelFittingrunNRLB.sh DataConfig.txt NRLBConfig.txt


## Run R package

R

library("ggseqlogo", lib.loc="/users/cainu5/Rpackages")
library("NRLBtools", lib.loc="/users/cainu5/Rpackages")

Models = load.mo























