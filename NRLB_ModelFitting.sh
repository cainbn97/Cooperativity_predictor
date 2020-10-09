#!/bin/bash

## Author: cainu5
## NRLB analysis 
## 09/01/20

(( TRACE )) && set -x
set -euo pipefail
trap 'Something went wrong. Script error on line #$LINENO.'  ERR


## Bring in arguments
TF=${1:?Expected a TF as argument #1}
# LINK1=${2:?Expected download link of 1st cycle as argument #2}
# LINK2=${3:?Expected download link of zero cycle as argument #3}

# cd ~/NRLB/testing

# ## ACQUIRE DATA FROM ENA DATABASE ##

# ## Check if folder already exists
# if [ -d "$TF" ]
# then
	# echo "Have you run NRLB for this dataset?"
	# # exit 1
# fi

# mkdir "$TF"; cd "$TF"
# curl -flOJ "$LINK1"

# # Check if zero cycle has been previously downloaded
# ZeroTag=$( echo $LINK2 | cut -d '_' -f 2 )
# if [ -d "$ZeroTag" ]
# then
	# #If zero tag was downloaded, only download TF file
	
	# echo "Zero cycle has been downloadeed"
# else
	# cd ~/NRLB/testing
	# mkdir "$ZeroTag"; cd "$ZeroTag"
	
	# #download files into TF folder and then move Zero cycle; so parallel can 
	# # be used
	# curl -flOJ "$LINK2"
	
# fi



# cd ~/NRLB/testing
# #check file integrities
# if [ gzip -t "$TF/$TF*.fastq.gz" ]
# then
	# echo "Download of TF file was successful"
# else
	# echo "Download of TF file was unsuccessful"
	# exit 1
# fi

# if [ gzip -t "$ZeroTag/*$ZeroTag*.fastq.gz" ]
# then
	# echo "Download of Zero tag file was successful"
# else 
	# echo "Download of Zero tag file was unsuccessful"
	# exit 1
# fi


## RUN NRLB MODEL FITTING (JAVA) ##


cd ~/NRLB/NRLB/ModelFitting

bash runNRLB.sh ~/NRLB/testing/"$TF"/Cycle4/DataConfig.txt ~/NRLB/testing/"$TF"/Cycle4/NRLBConfig.txt

mkdir ~/NRLB/testing/"$TF"/Cycle4/NRLBModels
mkdir ~/NRLB/testing/"$TF"/Cycle4/R0Models/

mv ~/NRLB/NRLB/ModelFitting/NRLBModels/"${TF}_NSBinding_4-NRLBConfig.txt" ~/NRLB/testing/"$TF"/Cycle4/NRLBModels/"${TF}_NSBinding_4-NRLBConfig.txt"
mv ~/NRLB/NRLB/ModelFitting/NRLBModels/"${TF}_NSBinding_4-NRLBConfig.list" ~/NRLB/testing/"$TF"/Cycle4/NRLBModels/"${TF}_NSBinding_4-NRLBConfig.list"
mv ~/NRLB/NRLB/ModelFitting/NRLBModels/"${TF}_NSBinding_4-NRLBConfig.dat" ~/NRLB/testing/"$TF"/Cycle4/NRLBModels/"${TF}_NSBinding_4-NRLBConfig.dat"
mv ~/NRLB/NRLB/ModelFitting/NRLBModels/"${TF}_NSBinding_4-NRLBConfig.csv" ~/NRLB/testing/"$TF"/Cycle4/NRLBModels/"${TF}_NSBinding_4-NRLBConfig.csv"


mv ~/NRLB/NRLB/ModelFitting/R0Models/"${TF}_NSBinding_4.txt" ~/NRLB/testing/"$TF"/Cycle4/R0Models/"${TF}_NSBinding_4.txt"
mv ~/NRLB/NRLB/ModelFitting/R0Models/"${TF}_NSBinding_4.dat" ~/NRLB/testing/"$TF"/Cycle4/R0Models/"${TF}_NSBinding_4.dat"

# Run R package

module load R/3.5.3-wrl

cd ~/NRLB/testing/"$TF"/Cycle3/NRLBModels/

Rscript ~/NRLB/code/NRLB_ModelVisualization.R "${TF}_NSBinding_3"

## Move pwms to folder
mkdir pwms
mv *motif* pwms

## Move logos to folder
mkdir motifs_logos
mv *.png motifs_logos




















