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
from Bio.Blast import NCBIWWW, NCBIXML

from mapping.datasets import MiSeq_runs
from mapping.filenames import get_divided_filenames
from mapping.mapping_utils import pair_generator



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Trim and divide reads into fragments')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--adaID', type=int, required=True,
                        help='Adapter ID to analyze (e.g. 4)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    miseq_run = args.run
    adaID = args.adaID
    VERBOSE = args.verbose

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    # Get the BAM filename 
    input_filename = get_divided_filenames(data_folder, adaID, ['F6'], type='bam')[-2]

    # Get unmapped reads and BLAST them
    with pysam.Samfile(input_filename, 'rb') as input_file:
        for reads in pair_generator(input_file):
            if reads[0].is_unmapped:
                seq = reads[0].seq

                # BLAST it
                seqb = Seq(seq, IUPAC.ambiguous_dna)
                blast_xml = NCBIWWW.qblast("blastn", "nr", seqb)
                blast_record = NCBIXML.read(blast_xml)
                ali = blast_record.alignments
                if len(ali):
                    ali = ali[0]
                    print ali.title
                else:
                    print 'No matches found'
