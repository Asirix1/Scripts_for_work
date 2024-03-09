#!/bin/bash
export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.5.0:$PATH
for i in Lamprey Gallus Human_GM12878
do 
mkdir /storage2/asirix/MOTIF/$i
cd /storage2/asirix/MOTIF/$i
wget wget -P ./ http://172.26.5.30/ftp/by_Project/HybridChromatin/$i/ctcf_chip/*.bigwig
wget wget -P ./ http://172.26.5.30/ftp/by_Project/HybridChromatin/$i/genome/*.f*a
bigWigToBedGraph ./*.bigwig ./$i'.bedGraph'
bed2fasta -o ./$1'_chip.fa' ./$i'.bedGraph' ./*.f*a 
meme-chip ./$1'_chip.fa' 
