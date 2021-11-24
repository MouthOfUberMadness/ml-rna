#!/bin/bash
sed -r 's/, complete genome//' $1.fasta >$1_1
sed -e 's/(/_/g' $1_1 > $1_2
rm $1_1
sed -e 's/)/_/g' $1_2 > $1_3
rm $1_2
sed '/^[^>]/s/T/U/g' $1_3 >$1_processed.fasta
rm $1_3
# sed -r 's/, complete genome//' $1_raw.fasta >$1.fasta
mafft --maxiterate 1000 --thread 12 $1_processed.fasta > $1_raw.aln
rm $1_processed.fasta
sed '/^[^>]/s/-/n/g' $1_raw.aln >$1.aln
rm $1_raw.aln