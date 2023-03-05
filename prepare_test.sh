#!/bin/sh
#This script will download the human genome fasta file 
#necessary for running test example of quaich_aging workflow
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gzip -d hg38.fa.gz
mv hg38.fa ./resources/genome/
