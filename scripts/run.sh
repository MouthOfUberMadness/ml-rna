#!/bin/bash

echo "Merbecovirus: starting codon (aligned), starting S (aligned), end S1 (aligned), start S2 (aligned), ending S (aligned) on HKU5-1"
./ml-rna 307 22099 24434 24580 26312 --analyse-phyl ../data/merbecoviruses.dnd

echo "Sarbecovirus: starting codon (aligned), starting S (aligned), end S1 (aligned), start S2 (aligned), ending S (aligned) on SarsCov2 Wuhan"
./ml-rna 296 21823 23981 24110 25789 --analyse-phyl ../data/bats_rodents_pangos_civet.dnd