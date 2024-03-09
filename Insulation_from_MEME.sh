#!/bin/bash
export PATH=$HOME/meme/bin:$HOME/meme/libexec/meme-5.5.0:$PATH
for t in Danio_rerio B_belcheri GalGal Mouse hg38 Lamprey C_punctatum Drosophila Oikopleura_dioica
do
g=$t'.fna'
h=$t'.hic'

samtools faidx /storage3/asirix/HumanDros/genomes/$g
mkdir /storage3/asirix/HumanDros/FIMO_out/$t
cd /storage3/asirix/HumanDros/FIMO_out/$t

fimo --max-stored-scores 25000 --motif 1 /storage3/asirix/HumanDros/Lamprey_MEME/combined.meme /storage3/asirix/HumanDros/genomes/$g
d='/storage3/asirix/HumanDros/FIMO_out/'$t'/fimo_out/fimo.tsv'

r=$(cat $d | wc -l )
let p=($r-4)
head -n $p $d | tail -n +2 > ./$t'_called_motifs'.tsv


awk 'BEGIN {OFS="\t"}; {print $3,$4,$5,$7}' ./$t'_called_motifs'.tsv | sort -k1,1 -k2,2n > ./$t'_called_motifs_sorted'.bedGraph
bedtools merge -i ./$t'_called_motifs_sorted'.bedGraph -c 4 -o mean > $t'_called_motifs_merged'.bedGraph
bedGraphToBigWig ./$t'_called_motifs_merged'.bedGraph /storage3/asirix/HumanDros/genomes/$g.fai ./$t'_called_motifs'.bigwig

cat ./$t'_called_motifs_sorted'.bedGraph > ./'Lamprey_CTCF_on_'$t.bed

mkdir /storage3/asirix/HumanDros/Figs/$t
cd /storage3/asirix/HumanDros/Figs/$t

cooler balance /storage3/asirix/HumanDros/cool/$t'.cool'
coolpup.py /storage3/asirix/HumanDros/cool/$t'.cool' '/storage3/asirix/HumanDros/FIMO_out/'$t'/Lamprey_CTCF_on_'$t.bed --local --nproc 0
plotpup.py --dpi 600 --input_pups ./*.clpy --no_score --output ../'Lamprey_CTCF_on_'$t.png

done

