#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Generate fastq files with set monomer and dimer contents
This is very much a brute force approach

cainu5
11/03/20
'''

import argparse
from random import choice
from random import shuffle
import os
import re

def rand_nuc(length):
    DNA = ''   
    for count in range(int(length)):
        DNA+=choice('CGTA')
    return DNA
    
def dimer():
    site = choice(monomer_sites)
    DNA = site + rand_nuc(Spacer) + site
    return DNA

def check_counts(sequence):
    nmatches = 0
    for match in re.finditer('TAAT', sequence):
        nmatches +=1
    for match in re.finditer('ATTA', sequence):
        nmatches +=1
    
    ## Take care of overlapping sequences
    for match in re.finditer('TAATAAT', sequence):
        nmatches +=1
    for match in re.finditer('ATTATTA', sequence):
        nmatches +=1
    for match in re.finditer('TAATTA', sequence):
        nmatches +=1
    for match in re.finditer('ATTAAT', sequence):
        nmatches +=1
    return nmatches

## parse arguments from command line
parser = argparse.ArgumentParser(description=
    'Generate fastq files with set monomer and dimer contents')
parser.add_argument('-N','--number', action = 'store', dest = 'Seq_tot', 
    type = int, default = 1000, help='Number of fastq sequences')
parser.add_argument('-M','--monomer', action = 'store', dest = 'Seq_mon', 
    type = int, default = 500, help='Number of monomer sites')
parser.add_argument('-D','--dimer', action = 'store', dest = 'Seq_dim', 
    type = int, default = 100, help='Number of dimer sites')
parser.add_argument('--name', action = 'store', dest = 'name', 
    type = str, default = 'mock_fastq', help = 'Enter the desired name of fastq file')
parser.add_argument('-L','--length', dest = 'Seq_len', 
    type = int, default = 20, help='Length of each sequence')
parser.add_argument('-S','--spacer', action = 'store', dest = 'spacer',
    type = int, default = 7, help = 'Enter desired dimer sequence')

args = parser.parse_args()
Seq_tot = args.Seq_tot
Seq_mon = args.Seq_mon
Seq_dim = args.Seq_dim
name = args.name
Seq_len = args.Seq_len
Spacer = args.spacer

## Output the inputs
print('Generating fastq.gz file with '+ str(Seq_tot) + ' sequences')
print('File will contain ' + str(Seq_mon) + ' monomer sites')
print('File will contain ' + str(Seq_dim) + ' dimer sites')
print('Dimer sites have a spacer of ' + str(Spacer) + ' bp')
print('Each strand will be ' + str(Seq_len) + ' bp')

monomer_sites = ['TAAT','ATTA']

with open('temp.fa', 'w') as fastq:
    ## Write monomer sites dispersed in random sequences
    for m in range(Seq_mon):
        ## Add random sequence with monomer
        ## Keeps shuffling until only single monomer site found
        site_count_match = False
        while site_count_match == False:
            left_flank = choice(range(Seq_len-4))
            right_flank = Seq_len - 4 - left_flank
            sequence = rand_nuc(left_flank)+choice(monomer_sites)+ rand_nuc(right_flank)+'\n'
            nmatches = check_counts(sequence)
            if nmatches == 1:
                site_count_match = True
        fastq.write(sequence)
    
    ## Write dimer sites dispersed in random sequences
    for m in range(Seq_dim):
        ## Add random sequence with dimer
        ## Keeps shuffling until only single dimer site found    
        site_count_match = False
        while site_count_match == False:
            left_flank = choice(range(Seq_len-(8+Spacer)))
            right_flank = Seq_len - (8+Spacer) - left_flank
            sequence = rand_nuc(left_flank)+dimer()+ rand_nuc(right_flank)+'\n'
            nmatches = check_counts(sequence)
            if nmatches == 2:
                site_count_match = True
        fastq.write(sequence)
    
    ## Write random strands that do not contain monomer or dimers    
    for m in range(Seq_tot - Seq_mon - Seq_dim):
        site_count_match = False
        while site_count_match == False:
            sequence = rand_nuc(Seq_len)+'\n'
            nmatches = check_counts(sequence)
            if nmatches == 0:
                site_count_match = True
        fastq.write(sequence)
            
        
## Open file and read it in for randomization
lines = open('temp.fa','r').readlines()
shuffle(lines)

## Write randomized file to fastq format
filename = name + '.fastq'
with open(filename, 'w') as fastq:
    for l in lines:
        
        ## Add dummy header
        fastq.write('@GAIIX:62VYNAAXX:1\n')
        ## Add sequence
        fastq.write(l)
        ## Add dummy quality scores
        fastq.write('IIIIIIIIIIIIIIIIIGII\n+\n')

## Remove temp file
os.remove('temp.fa')

## Confirm output
file = open(filename, encoding = 'utf8')
file = file.read().strip()

## Check total sequences
tot = '[G|A|T|C]{'+str(Seq_len)+'}'
nmatches = 0
for match in re.finditer(tot, file):
    nmatches += 1

if Seq_tot == nmatches:
    print('Total sequences match input')
else:
    print('Total sequences do not match input.')
    print('Input: '+ str(Seq_tot))
    print('Output: '+ str(nmatches))
    
## Check dimer sequences
dim_1 = 'TAAT[G|A|T|C]{'+str(Spacer)+'}TAAT'
dim_2 = 'ATTA[G|A|T|C]{'+str(Spacer)+'}ATTA'
nmatches = 0
for match in re.finditer(dim_1, file):
    nmatches += 1

for match in re.finditer(dim_2, file):
    nmatches +=1

if Seq_dim == nmatches:
    print('Dimer sites match input')
else:
    print('Dimer sites do not match input.')
    print('Input: '+ str(Seq_dim))
    print('Output: '+ str(nmatches))
 
## Check monomer sequences
nmatches = 0
for match in re.finditer('TAAT', file):
    nmatches +=1
    
for match in re.finditer('ATTA', file):
    nmatches +=1
    
if Seq_mon == (nmatches - Seq_dim*2):
    print('Monomer sites match input')
else:
    print('Monomer sites do not match input')
    print('Input: '+ str(Seq_mon))
    print('Output: '+ str(nmatches-Seq_dim*2))
