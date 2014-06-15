# vim: fdm=marker
'''
author:     Fabio Zanini
date:       15/06/14
content:    Basic script to obtain the data from a seq run, e.g. adapters, sample
            names, etc.
'''
# Modules
import argparse

from hivwholeseq.samples import SampleSeq, load_sequencing_run



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Demultiplex HIV reads')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28, test_tiny)')
    
    args = parser.parse_args()
    seq_run = args.run

    dataset = load_sequencing_run(seq_run)
