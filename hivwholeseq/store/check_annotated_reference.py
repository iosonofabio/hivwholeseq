# vim: fdm=marker
'''
author:     Fabio Zanini
date:       27/01/15
content:    Check the annotated genomewide reference for problems in the proteins,
            et al.
'''
# Modules
def check_protein(fea, ref, VERBOSE=0):
    '''Check a protein annotation'''
    #TODO



# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Check patient samples')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')

    args = parser.parse_args()
    VERBOSE = args.verbose
    pnames = args.patients

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
                check_protein(fea, ref, VERBOSE=VERBOSE)
