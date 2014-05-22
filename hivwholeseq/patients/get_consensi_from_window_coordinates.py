# vim: fdm=marker
'''
author:     Fabio Zanini
date:       22/05/14
content:    Get an alignment of consensi from all time points of a patient,
            from a window selected by coordinates in a reference.
'''
# Modules
import os
import argparse
from itertools import izip
import numpy as np
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import ambiguous_dna

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.patients import patients as patients_all
from hivwholeseq.patients.filenames import get_initial_consensus_filename, \
        get_foldername, get_coordinate_map_filename, \
        get_allele_count_trajectories_filename



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Update initial consensus')
    parser.add_argument('--reference', default='HXB2',
                        help='Select reference strain to align against')
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--coordinates', required=True, type=int, nargs=2,
                        help='Coordinates of window (BOTH extremes includes)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--save', action='store_true',
                        help='Save the allele frequency trajectories to file')

    args = parser.parse_args()
    refname = args.reference
    pnames = args.patients
    coords = args.coordinates
    VERBOSE = args.verbose
    save_to_file = args.save
    
    # use Python coordinates
    coords[1] += 1

    if pnames is None:
        patients = [p for p in patients_all]
    else:
        patients = [p for p in patients_all if p.id in pnames]    
    if not len(patients):
        raise ValueError('No patients found!')

    if VERBOSE >= 3:
        print 'patients', map(attrgetter('id'), patients)

    for patient in patients:
        pname = patient.id

        found = False
        for fragment in ['F'+str(i+1) for i in xrange(7)] + ['genomewide']:
            map_fn = get_coordinate_map_filename(pname, fragment, refname=refname)
            mapco = np.loadtxt(map_fn, dtype=int, unpack=True)
            frag_start = mapco[0, 0]
            frag_end = mapco[0, -1] + 1

            if (frag_start <= coords[0]) and (frag_end + 1 >= coords[1]):
                found = True
                break

        if not found:
            raise ValueError('Window not found in any fragment or genomewide')

        if coords[0] in mapco[0]:
            win_start_init = mapco[1, mapco[0] == coords[0]][0]
        else:
            found = False
            coord_start = coords[0] - 1
            while coord_start >= frag_start:
                if coord_start in mapco[0]:
                    win_start_init = mapco[1, mapco[0] == coord_start][0]
                    found = True
                    break
                coord_start -= 1

            if not found:
                raise ValueError('Start not found (this should NOT happen)?!')

        if (coords[1] - 1) in mapco[0]:
            win_end_init = mapco[1, mapco[0] == (coords[1] - 1)][0] + 1
        else:
            found = False
            coord_end = coords[1] + 1
            while coord_end <= frag_end:
                if coord_end - 1 in mapco[0]:
                    win_end_init = mapco[1, mapco[0] == (coord_end - 1)][0] + 1
                    found = True
                    break
                coord_end += 1

            if not found:
                raise ValueError('End not found (this should NOT happen)?!')

        act_filename = get_allele_count_trajectories_filename(pname, fragment)
        act = np.load(act_filename)
        consensi = []
        samplenames = patient.samples
        for (samplename, ac) in izip(samplenames, act):
            ac = ac[:, win_start_init: win_end_init]
            cons = alpha[ac.argmax(axis=0)]
            cons[ac.sum(axis=0) == 0] = 'N'
            cons_rec = SeqRecord(Seq(''.join(cons), ambiguous_dna),
                                 id=pname+'_'+samplename+'_'+refname+'_'+'-'.join(map(str, coords)),
                                 name=pname+'_'+samplename+'_'+refname+'_'+'-'.join(map(str, coords)),
                                 description=pname+'_'+samplename+'_'+refname+'_'+'-'.join(map(str, coords))+', consensus')
            consensi.append(cons_rec)
        
        if save_to_file:
            SeqIO.write(consensi,
                        '/ebio/ag-neher/home/fzanini/tmp/seqs_'+pname+'_'+'-'.join(map(str, coords))+'.fasta',
                        'fasta')

        




