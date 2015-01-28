# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/01/15
content:    Check the annotated genomewide reference for problems in the proteins,
            et al.
'''
# Modules
import argparse

from hivwholeseq.patients.patients import load_patients, Patient


# Functions
def check_protein(fea, seqgw, VERBOSE=0, delta_pos=2.5):
    '''Check a protein annotation'''
    seq = fea.extract(seqgw).seq

    if len(seq) % 3:
        raise ValueError('The length of '+fea.id+' is not a multiple of 3')

    if 'N' in seq:
        raise ValueError('N nucleotides found in '+fea.id)

    if '-' in seq:
        raise ValueError('Gaps found in '+fea.id)

    prot = seq.translate()

    if ('*' in prot) and (prot.find('*') != len(prot) - 1):
        raise ValueError('Premature stops found in '+fea.id)

    if 'X' in prot:
        raise ValueError('X amino acids found in '+fea.id)

    # Compare to HXB2
    from hivwholeseq.reference import load_custom_reference
    ref = load_custom_reference('HXB2', region=fea.id)

    from seqanpy import align_global
    (score, alis, alir) = align_global(seq, ref, score_gapopen=-20)
    if VERBOSE >= 3:
        from hivwholeseq.utils.sequence import pretty_print_pairwise_ali
        pretty_print_pairwise_ali((alir, alis), name1='HXB2', name2='seq',
                                  width=100)

    scoremax = 3 * len(alis)
    delta = scoremax - score
    if delta > delta_pos * len(alis):
        raise ValueError('The sequence of '+fea.id+' looks different from HXB2')



# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Check patient samples')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--force', action='store_true',
                        help='Do not stop for errors')

    args = parser.parse_args()
    VERBOSE = args.verbose
    pnames = args.patients
    use_force = args.force

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[pnames]


    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        if VERBOSE >= 1:
            print 'Patient:', patient.name

        ref = patient.get_reference('genomewide', 'gb')

        for fea in ref.features:
            if fea.type == 'protein':
                if VERBOSE >= 2:
                    print 'Checking', fea.id
                try:
                    check_protein(fea, ref, VERBOSE=VERBOSE)
                except ValueError:
                    if use_force:
                        print 'ERROR!'
                    else:
                        raise

