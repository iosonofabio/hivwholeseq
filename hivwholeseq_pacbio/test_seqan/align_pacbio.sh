#!/bin/sh
refseq=$(cat NL4-3.fasta)
readseq=$(head -n 2062 /ebio/ag-neher/share/data/PacBio_HIV_Karolinska/run23/ccs_reads/pb_023_1_ccs_reads.fastq.txt | tail -n 1)

./seqan $refseq $readseq

