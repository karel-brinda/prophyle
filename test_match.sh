#! /usr/bin/env bash -f -e -o pipefail

fq=tests/simulation_bacteria.1000.fq

echo
echo " =========== 1 ============"
echo

./prophyle/prophyle.py classify -P _index_test $fq | samtools view -b > _test_bam1.bam

echo
echo " =========== 2 ============"
echo

./prophyle/prophyle.py classify -P _index_test $fq $fq | samtools view -b > _test_bam2.bam

