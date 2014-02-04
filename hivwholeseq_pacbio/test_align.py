# vim: fdm=indent
'''
author:     Fabio Zanini
date:       24/01/14
content:    Test parsing the PacBio data.
'''
# Modules
import os
from Bio import SeqIO
import pysam
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt



# Globals
data_folder = '/ebio/ag-neher/share/data/PacBio_HIV_Karolinska/run23/'
samples = [{'name': 'S1', 'description': 'NL4-3', 'filename': 'pb_023_1'},
           {'name': 'S2', 'description': 'Mix1', 'filename': 'pb_023_2'},
           {'name': 'S3', 'description': 'pat', 'filename': 'pb_023_3'},
           ]
samples = pd.DataFrame(samples)


# Functions
def trim_reference(refm, readm, band=100, VERBOSE=0):
    '''Trim the reference around a seed of the read'''

    seed_len = 30
    matches_min = 25
    rl = len(readm)

    # Try out a few seeds: middle of the read, a bit earlier, a bit later
    seed_starts = [max(0, rl / 2 - seed_len / 2),
                   max(0, rl / 2 - 5 * seed_len),
                   min(rl - seed_len, rl / 2 + 5 * seed_len)]

    for seed_start in seed_starts:
        seed = readm[seed_start: seed_start + seed_len]
        mismatches = 0

        # Get scores for the seed
        scores = [(seed == refm[i: i + seed_len]).sum()
                  for i in xrange(len(refm) - (rl - seed_start) + band / 2)]
        pos_seed = np.argmax(scores)

        if scores[pos_seed] >= matches_min:
            if VERBOSE >= 1:
                print 'Position:', pos_seed, 'score:', scores[pos_seed]
        
            pos_trim_start = pos_seed - seed_start - band / 2
            pos_trim_end = pos_seed + (rl - seed_start) + band / 2
            ref_trim = refm[max(0, pos_trim_start): min(len(refm), pos_trim_end)]

            if pos_trim_end > len(ref_trim):
                ref_trim = np.concatenate([ref_trim, np.repeat('-', pos_trim_end - len(ref_trim))])
            if pos_trim_start < 0:
                ref_trim = np.concatenate([np.repeat('-', -pos_trim_start), ref_trim])

            return (pos_trim_start, ref_trim)
    
    raise ValueError('Seed not found in reference')


def align_overlap_seqan(refm, readm, band=100, VERBOSE=0, raw_output=False):
    '''Global alignment via SeqAn'''
    from Bio.Seq import reverse_complement as revcom

    if len(refm) > len(readm) + 2 * band:
        try:
            (pos, refm) = trim_reference(refm, readm, band=band, VERBOSE=VERBOSE)
        except ValueError:
            try:
                readm = np.fromstring(revcom(readm.tostring()), 'S1')
                (pos, refm) = trim_reference(refm, readm, band=band, VERBOSE=VERBOSE)
            except ValueError:
                raise ValueError('Seed not found in reference on any strand')

    # Call SeqAn
    # OLD WAY: pipes
    if False:
        import subprocess as sp
        output = sp.check_output(['./test_seqan/align_global',
                                  refm.tostring(),
                                  readm.tostring(),
                                  str(band)])
        output = output.split('\n')
        if len(output) < 10:
            raise ValueError('=Overlap alignment failed')
    
        # Parse the horrible output
        ali1 = ''.join([output[5 * i + 3].strip(' ') for i in xrange(len(output[1:]) // 5)])
        ali2 = ''.join([output[5 * i + 5].strip(' ') for i in xrange(len(output[1:]) // 5)])

    # NEW WAY: wrapper
    else:
        import seqan_module.seqanpy as sap
        (score, ali1, ali2) = sap.align_overlap(refm.tostring(),
                                                readm.tostring(),
                                                band=band,
                                                trim=True)
        pos += score

    if raw_output:
        return (pos, ali1, ali2)

    # Make a Biopython MSA
    from Bio.Align import MultipleSeqAlignment as MSA
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet.IUPAC import ambiguous_dna
    ali = MSA([SeqRecord(Seq(ali1, alphabet=ambiguous_dna), id='ref'),
               SeqRecord(Seq(ali2, alphabet=ambiguous_dna), id='read')])

    return (pos, ali)




# Script
if __name__ == '__main__':

    samplename = 'S1'

    # Parse consensus reads
    sample = samples.set_index('name').loc[samplename]

    reads_iter = SeqIO.parse(data_folder+'ccs_reads/'+sample['filename']+\
                             '_ccs_reads.fastq.txt', 'fastq')
    
    reads_all = [reads_iter.next() for i in xrange(2100)]
    reads = reads_all[2000:2100]
    

    # Get NL4-3 reference
    from hivwholeseq.reference import load_NL43  
    refseq = load_NL43()
    refm = np.array(refseq)

    # Align a few reads
    from hivwholeseq.mapping_utils import align_muscle
    poss = []
    alis = []
    for i, read in enumerate(reads_iter):

        print i,

        readm = np.array(read)
        try:
            (pos, ali) = align_overlap_seqan(refm, readm)
        except ValueError:
            pos = ali = None

        print ali

        poss.append(pos)
        alis.append(ali)

