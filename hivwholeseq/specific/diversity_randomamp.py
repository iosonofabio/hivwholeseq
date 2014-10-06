# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/10/14
content:    Try to get an idea about diversity in a randomly amplified sample as
            opposed to fragment-based PCR.
'''
# Modules
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt

from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file_unfiltered
from hivwholeseq.sequencing.samples import load_sample_sequenced
from hivwholeseq.miseq import alpha
from hivwholeseq.sequencing.primer_info import find_fragment




# Script
if __name__ == '__main__':


    sample = load_sample_sequenced('Chris')
    refseq = SeqIO.read(sample.get_reference_premap_filename(), 'fasta')
    bamfilename = sample.get_premapped_filename()

    counts, inserts = get_allele_counts_insertions_from_file_unfiltered(bamfilename, len(refseq), VERBOSE=2)

    # Average over read types
    counts = counts.sum(axis=0) 
    cov = counts.sum(axis=0)
    consm = alpha[counts.argmax(axis=0)]
    consm[cov < 30] = 'N'
    conss = ''.join(consm)

    pos_fragments = {fr: find_fragment(consm, fr, threshold=0.6)
                     for fr in ('F2o', 'F3bo', 'F4i')}
   
    from hivwholeseq.patients.patients import load_sample_sequenced as lssp
    sample_29184 = lssp('29184')
    cons_29184_F3 = ''.join(sample_29184.get_consensus('F3'))

    from seqanpy import align_overlap
    (score, ali1, ali2) = align_overlap(conss[100: -200], cons_29184_F3)
    start = len(ali2) - len(ali2.lstrip('-'))
    end = len(ali2.rstrip('-'))
    from hivwholeseq.sequence_utils import pretty_print_pairwise_ali
    pretty_print_pairwise_ali([ali1[start: end], ali2[start:end]])


    fig, ax = plt.subplots()
    ax.plot(cov, lw=2)
    ax.set_xlabel('Position in HIV genome [bp]')
    ax.set_ylabel('Coverage')
    ax.grid(True)
    ax.set_title('Chris\'s random amplification')

    plt.ion()
    plt.show()
