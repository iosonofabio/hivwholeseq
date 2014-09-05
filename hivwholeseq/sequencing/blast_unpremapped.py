# vim: fdm=marker
'''
author:     Fabio Zanini
date:       04/11/13
content:    Take the reads that scored unmapped in the premapping and BLAST them
            to get an idea bout contaminants, and make sure we are not throwing
            away useful stuff.
'''
# Modules
import argparse
import pysam
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIWWW, NCBIXML

from hivwholeseq.sequencing.filenames import get_premapped_filename
from hivwholeseq.mapping_utils import pair_generator, reads_to_seqrecord

from hivwholeseq.sequencing.samples import load_sequencing_run


# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Trim and divide reads into fragments')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaID', required=True,
                        help='Adapter ID to analyze (e.g. TS4)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    seq_run = args.run
    adaID = args.adaID
    VERBOSE = args.verbose

    # Specify the dataset
    dataset = load_sequencing_run(seq_run)
    data_folder = dataset.folder

    # Get the BAM filename 
    input_filename = get_premapped_filename(data_folder, adaID, type='bam')

    # Get unmapped reads and BLAST them
    reads_unmapped = []
    n_unmapped = 0
    with pysam.Samfile(input_filename, 'rb') as input_file:
        for reads in pair_generator(input_file):
            if not reads[0].is_unmapped:
                continue

            n_unmapped += 1

            # Take only the first part of read1, to make sure quality is high
            seq = reads[reads[1].is_read1].seq[:200]
            seqb = Seq(seq, IUPAC.ambiguous_dna)

            # Save to file, to test local blast
            reads_unmapped.append(reads[reads[1].is_read1])
            if len(reads_unmapped) >= 100:
                break

            continue

            # BLAST it
            blast_xml = NCBIWWW.qblast("blastn", "nr", seqb)
            blast_record = NCBIXML.read(blast_xml)
            ali = blast_record.alignments
            if len(ali):
                ali = ali[0]
                print ali.title
            else:
                print 'No matches found'

        seqs_unmapped = reads_to_seqrecord(reads_unmapped)

    from Bio import SeqIO
    SeqIO.write(seqs_unmapped, '/ebio/ag-neher/home/fzanini/tmp/seqs_for_blast.fastq', 'fastq')
    SeqIO.write(seqs_unmapped, '/ebio/ag-neher/home/fzanini/tmp/seqs_for_blast.fasta', 'fasta')

