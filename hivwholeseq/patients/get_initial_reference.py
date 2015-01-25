# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/03/14
content:    Load the initial reference for manual analysis.
'''
# Modules
import argparse

from hivwholeseq.patients.patients import load_patient



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get initial patient reference')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--region', required=True,
                        help='Regions to analyze (e.g. V3 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    pname = args.patient
    region = args.region
    VERBOSE = args.verbose

    patient = load_patient(pname)
    ref = patient.get_reference(region)
