# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/12/14
content:    Check the V3 alignment in p4/15313, as it looks weird.
'''
# Modules
import os
from operator import itemgetter
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import Phylo

from hivwholeseq.patients.patients import load_patient
from hivwholeseq.sequence_utils import align_muscle
from hivwholeseq.tree_utils import build_tree_fasttree



# Functions
def load_alignments(filename):
    '''Load alignments from website file'''
    import zipfile, zlib
    from Bio import AlignIO
    import StringIO

    alis = []
    with zipfile.ZipFile(filename, 'r') as zf:
        for fn in zf.namelist():
            f = StringIO.StringIO(zf.read(fn))
            ali = {'time': float(fn.split('_')[0]),
                   'ali': AlignIO.read(f, 'fasta')}
            alis.append(ali)

    return alis


def get_region_count_trajectories(patient, region, VERBOSE=0, countmin=5):
    '''Get haplotype trajectories in a region (from the website alignments)'''
    import numpy as np
    from hivwholeseq.website.filenames import get_precompiled_alignments_filename

    filename = get_precompiled_alignments_filename(patient.code, region)
    alis = load_alignments(filename)
    
    seqs_set = set()
    for ali in alis:
        seqs_set |= set([''.join(seq).replace('-', '')
                         for seq in ali['ali']
                         if int(seq.name.split('_')[1]) >= countmin])
    seqs_set = list(seqs_set)
    
    hct = np.zeros((len(seqs_set), len(alis)), int)
    for it, ali in enumerate(alis):
        for seq in ali['ali']:
            s = ''.join(seq).replace('-', '')
            count = int(seq.name.split('_')[1])
            if count < countmin:
                continue
            iseq = seqs_set.index(s)
            hct[iseq, it] = count

    seqs_set = np.array(seqs_set, 'S'+str(np.max(map(len, seqs_set))))
    times = np.array(map(itemgetter('time'), alis))
    ind = np.array([i for i, t in enumerate(patient.times) if t in times])
    
    return (hct.T, ind, seqs_set)



# Script
if __name__ == '__main__':

    VERBOSE = 2
    pname = '15313'
    regions = ['psi', 'PR', 'V3']
    countmin = 5

    patient = load_patient(pname)

    for region in regions:
    
        (hct, ind, seqs) = get_region_count_trajectories(patient, region,
                                                         countmin=countmin,
                                                         VERBOSE=VERBOSE)

        times = patient.times[ind]

        hct = hct.T
        hft = 1.0 * hct / hct.sum(axis=0)

        seqsali = align_muscle(*seqs, sort=True)
        tree = build_tree_fasttree(seqsali, VERBOSE=VERBOSE)
        tree.ladderize()
    
        fig, ax = plt.subplots()
        Phylo.draw(tree, axes=ax)
        ax.set_title(patient.code+', '+region)
    
        fig, ax = plt.subplots()
        for hf in hft:
            ax.plot(times, hf, lw=2)
    
        ax.set_xlabel('Time [days]')
        ax.set_ylabel('Haplotype frequency')
    
        plt.tight_layout()

        import sys
        sys.exit()

    plt.ion()
    plt.show()




