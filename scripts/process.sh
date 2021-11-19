#!/bin/bash
sed -r 's/, complete genome//' sarbecoviruses_raw.fasta >sarbecoviruses_1
sed '/^[^>]/s/T/U/g' sarbecoviruses_1 >sarbecoviruses.fasta
# sed -r 's/, complete genome//' sarbecoviruses_raw.fasta >sarbecoviruses.fasta
mafft sarbecoviruses.fasta > sarbecoviruses_raw.aln
sed '/^[^>]/s/-/n/g' sarbecoviruses_raw.aln >sarbecoviruses.aln
rm sarbecoviruses_1
rm sarbecoviruses.fasta
rm sarbecoviruses_raw.aln