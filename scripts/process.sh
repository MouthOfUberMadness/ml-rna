#!/bin/bash
sed -r 's/, complete genome//' $1.fasta >$1_1
sed '/^[^>]/s/T/U/g' $1_1 >$1_processed.fasta
# sed -r 's/, complete genome//' $1_raw.fasta >$1.fasta
mafft $1_processed.fasta > $1_raw.aln
sed '/^[^>]/s/-/n/g' $1_raw.aln >$1.aln
rm $1_1
rm $1_processed.fasta
rm $1_raw.aln