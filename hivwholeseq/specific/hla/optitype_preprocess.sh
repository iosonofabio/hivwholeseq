#!/bin/sh
################################################################
#
# Map the human reads from the HIV cell samples to the HLA
# cluster to see how many we have and possibly use OptiType to
# get the HLA type of the patients.
#
################################################################
if [ "$#" -ne 3 ]; then
 echo "Please call this script with 3 parameters: pname, either --submit or --run, and number of threads."
 exit 3
fi

pname="$1"
use_submit="$2"
NCORES="$3"

echo "pname: $pname"
echo "ncores: $NCORES"

# Submit self to the cluster via qsub if requested
if [ "$use_submit" = --submit ]; then
 HIVWHOLESEQFOLDER="/ebio/ag-neher/share/users/fzanini/phd/sequencing/scripts/mapping/hivwholeseq"
 JOBLOGOUT="$HIVWHOLESEQFOLDER/cluster/logout"
 JOBLOGERR="$HIVWHOLESEQFOLDER/cluster/logerr"
 JOBSCRIPT="$HIVWHOLESEQFOLDER/specific/hla/optitype_preprocess.sh"
 VMEM=8G
 RT=23:59:59
 qsub -cwd -b y -S '/bin/bash' -o $JOBLOGOUT -e $JOBLOGERR -N hla$pname -pe parallel "$NCORES" -l h_rt=$RT -l h_vmem=$VMEM $JOBSCRIPT $pname --run $NCORES
 exit 0
fi

datapath="/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/specific/HLA/$pname"

read1path="$datapath/read1.fastq.gz"
read2path="$datapath/read2.fastq.gz"
refpath="/ebio/ag-neher/home/fzanini/programs/OptiType/data/hla_reference_dna.fasta"
outputpath="$datapath/sample_finished.sam"

# NOTE: --unique = -pa -m 1 -dr 0, which is what Oliver suggests
razers3 -v -tc $NCORES -ll 450 -le 100 --percent-identity 90 -m 1 -dr 0 --full-readid --output $outputpath $refpath $read1path $read2path
