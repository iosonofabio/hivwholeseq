# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/10/14
content:    Try to get an idea about diversity in a randomly amplified sample as
            opposed to fragment-based PCR.
'''
# Modules
from collections import defaultdict
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt

from hivwholeseq.utils.one_site_statistics import get_allele_counts_insertions_from_file_unfiltered
from hivwholeseq.sequencing.samples import load_sample_sequenced
from hivwholeseq.miseq import alpha
from hivwholeseq.sequencing.primer_info import find_fragment
from hivwholeseq.utils.sequence import pretty_print_pairwise_ali

from seqanpy import align_overlap


# Globals
coord = {'F3': [3500, 5500], 'F2': [1500, 4000]}




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

    #pos_F3 = find_fragment(consm, 'F3bo', threshold=0.7)
   
    from hivwholeseq.patients.patients import load_samples_sequenced as lssp
    from hivwholeseq.patients.patients import SamplePat
    samples_pat = lssp(patients=['p6'])

    potential_cont = defaultdict(list)
    for samplename_pat, sample_pat in samples_pat.iterrows():
        sample_pat = SamplePat(sample_pat)
        print sample_pat.name

        for fragment, coord_frag in coord.iteritems():
            cons_pat = ''.join(sample_pat.get_consensus(fragment))


            #sample_29184 = lssp('29184')
            #cons_29184_F3 = ''.join(sample_29184.get_consensus('F3'))

            (score, ali1, ali2) = align_overlap(conss[coord_frag[0]: coord_frag[1]], cons_pat)
            start = len(ali2) - len(ali2.lstrip('-'))
            end = len(ali2.rstrip('-'))
            scoremax = 3 * (end - start)
            deltamax = 6 * 30
            
            if score >= scoremax - deltamax:
                potential_cont[fragment].append((samplename_pat,
                                                 (score, ali1[start: end], ali2[start: end])))
                print samplename_pat
                pretty_print_pairwise_ali([ali1[start: end], ali2[start:end]])


    #fig, ax = plt.subplots()
    #ax.plot(cov, lw=2)
    #ax.set_xlabel('Position in HIV genome [bp]')
    #ax.set_ylabel('Coverage')
    #ax.grid(True)
    #ax.set_title('Chris\'s random amplification')

    #plt.ion()
    #plt.show()
