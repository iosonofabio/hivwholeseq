# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/09/14
content:    Extract features from annotated initial patient refernces, and
            align them across patients to have a look.
'''
# Modules
import argparse
from operator import attrgetter, itemgetter
from Bio import SeqIO
from Bio.Alphabet.IUPAC import ambiguous_dna

from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.sequence_utils import correct_genbank_features_load



# Functions
def collect_features():
    '''Collect feature names'''
    features = []
    from hivwholeseq.genome_info import genes
    features.extend(genes)

    return features





# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Annotate initial reference')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')

    args = parser.parse_args()
    VERBOSE = args.verbose
    pnames = args.patients

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[patients.index.isin(pnames)]

    refseqs = {}
    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        if VERBOSE:
            print 'Patient:', patient.name

        fn = patient.get_reference_filename('genomewide', 'gb')
        refseq = SeqIO.read(fn, 'gb', alphabet=ambiguous_dna)
        correct_genbank_features_load(refseq)
        refseqs[pname] = refseq

    features = collect_features()
    for featname in features:
        print featname
        for pname, refseq in refseqs.iteritems():
            for feature in refseq.features:
                if feature.id == featname:
                    break
            else:
                continue

            print '{:>5s}'.format(pname), feature.extract(refseq).seq.translate()

        print ''
