# vim: fdm=marker
'''
author:     Fabio Zanini
date:       24/03/14
content:    Study the frequency of known drug resistance mutations, e.g. in pol.
'''
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import ambiguous_dna

from hivwholeseq.miseq import alpha
from hivwholeseq.patients.patients import get_patient
from hivwholeseq.patients.filenames import get_initial_reference_filename, \
        get_mapped_to_initial_filename, get_allele_frequency_trajectories_filename, \
        get_allele_count_trajectories_filename


# Globals
prot_gene_dict = {'RT': 'pol',
                  'PR': 'pol'}

# FIXME: do better than this!
prot_pos = {'PR': [55 * 3, (55 + 300) * 3]}

mutations_DR = {'RT': ['K70R', 'T215Y'],
                'PR': ['L10I', 'K20R', 'L33F', 'M36I', 'M46I', 'I54V', 'A71V', 'V77I', 'V82A', 'L90M']}



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Study drug resistance mutations')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--protein', required=True,
                        help='Protein to analyze (e.g. RT)')

    args = parser.parse_args()
    pname = args.patient
    VERBOSE = args.verbose
    protein = args.protein

    patient = get_patient(pname)

    gene = prot_gene_dict[protein]

    # HIV makes polyproteins: each protein has coordinates within its gene 
    (start, end) = prot_pos[protein]

    aft = np.load(get_allele_frequency_trajectories_filename(pname, gene))[:, :, start: end]
    consm = alpha[aft[0].argmax(axis=0)]
    conss = consm.tostring()
    prot = Seq(conss, ambiguous_dna).translate()

    for mutstr in mutations_DR[protein]:
        al_anc = mutstr[0]
        al_mut = mutstr[-1]
        pos_prot = int(mutstr[1:-1])
        pos_gene = 3 * pos_prot

        # Check whether the consensus of the patient starts mutated already
        cod_anc_pat = Seq(conss[pos_gene: pos_gene + 3], ambiguous_dna)
        al_anc_pat = cod_anc_pat.translate()[0]
        if al_anc_pat != al_anc:
            print mutstr, 'the patient initial consensus has a mutated allele already:', al_anc_pat
            continue

        for pos_in_cod in xrange(3):
            for ia, nuc in enumerate(alpha[:4]):
                if nuc == cod_anc_pat[pos_in_cod]:
                    continue
                cod_tmp = list(cod_anc_pat)
                cod_tmp[pos_in_cod] = nuc
                cod_tmp = ''.join(cod_tmp)
                al_tmp = Seq(cod_tmp, ambiguous_dna).translate()[0]
                if al_tmp == al_mut:
                    aft_tmp = aft[:, ia, pos_gene + pos_in_cod]
                    aft_tmp[aft_tmp < 1e-3] = 0
                    if not (aft_tmp > 1e-3).any():
                        aft_tmp = 'ABSENT'
                    print mutstr, ''.join(cod_anc_pat), '->', cod_tmp, aft_tmp
