#!/bin/sh
################################################################
#
# Map the human reads from the HIV cell samples to the HLA
# cluster to see how many we have and possibly use OptiType to
# get the HLA type of the patients.
#
################################################################
if [ "$#" -ne 2 ]; then
 echo "Please call this script with two parameters, the patient name and either --submit or --run."
 exit 3
fi

pname="$1"

# Submit self to the cluster via qsub if requested
use_submit="$2"
if [ "$use_submit" = --submit ]; then
 HIVWHOLESEQFOLDER="/ebio/ag-neher/share/users/fzanini/phd/sequencing/scripts/mapping/hivwholeseq"
 JOBLOGOUT="$HIVWHOLESEQFOLDER/cluster/logout"
 JOBLOGERR="$HIVWHOLESEQFOLDER/cluster/logerr"
 JOBSCRIPT="$HIVWHOLESEQFOLDER/specific/hla/optitype_preprocess.sh"
 qsub -cwd -b y -S '/bin/bash' -o $JOBLOGOUT -e $JOBLOGERR -N hla$pname -l h_rt=23:59:59 -l h_vmem=4G $JOBSCRIPT $pname --run
 exit 0
fi

datapath="/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/specific/HLA/$pname"

read1path="$datapath/read1.fastq.gz"
read2path="$datapath/read2.fastq.gz"
refpath="/ebio/ag-neher/home/fzanini/programs/OptiType/data/hla_reference_dna.fasta"
outputpath="$datapath/sample_finished.sam"

razers3 -vv -ll 450 -le 100 --percent-identity 90 --max-hits 1 --distance-range 0 --output $outputpath $refpath $read1path $read2path
