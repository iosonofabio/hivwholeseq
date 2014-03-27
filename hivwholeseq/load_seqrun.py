# vim: fdm=marker
'''
author:     Fabio Zanini
date:       26/03/14
content:    Load sequencing run data for manual inspection.
'''
# Modules
import argparse

from hivwholeseq.datasets import MiSeq_runs


# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Load sequencing run')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')

    args = parser.parse_args()
    seq_run = args.run

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

