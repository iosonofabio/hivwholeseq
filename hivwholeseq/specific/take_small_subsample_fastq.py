# vim: fdm=marker
'''
author:     Fabio Zanini
date:       03/12/14
content:    Take a few reads for testing.
'''
# Modules
import gzip
import numpy as np
from Bio import SeqIO

from hivwholeseq.sequencing.filenames import get_read_filenames
from hivwholeseq.sequencing.samples import load_sequencing_run
from hivwholeseq.sequencing.adapter_info import adapters_illumina


# Script
if __name__ == '__main__':

    dataset = load_sequencing_run('Tue28')
   
    adaIDs = ['TS2', 'TS4', 'TS7']

    seqs = []
    for adaID in adaIDs:
        adastr = adapters_illumina[adaID]
        adaqual = list(np.random.randint(30, 39, size=6))

        fns = get_read_filenames(dataset.folder, adaID, gzip=True)

        with gzip.open(fns[0], 'rb') as fr1:
            for i, seq1 in enumerate(SeqIO.parse(fr1, 'fastq')):
                if i == 50:
                    break

                qual = seq1.letter_annotations.pop('phred_quality')
                qual = adaqual + qual
                seq1 = adastr + seq1
                seq1.letter_annotations['phred_quality'] = qual
                seqs.append(seq1)

    np.random.shuffle(seqs)
    SeqIO.write(seqs, '/ebio/ag-neher/home/fzanini/tmp/nextgen_seqs.fastq', 'fastq')
