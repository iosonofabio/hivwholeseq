# vim: fdm=marker
'''
author:     Fabio Zanini
date:       19/05/14
content:    Map the emPCR reads, they use totally different primers etc.
'''
# Modules
import os
import gzip
import numpy as np
from itertools import izip
import argparse
import subprocess as sp
from Bio.SeqIO.QualityIO import FastqGeneralIterator as FGI

from hivwholeseq.filenames import get_read_filenames
from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.mapping_utils import stampy_bin, get_number_reads_fastq_open
from hivwholeseq.adapter_info import foldername_adapter



# Globals
sd = {'mixconv1542': {'adaID': 'TS2'},
      'mixconv1968': {'adaID': 'TS5'},
      'mixem1542': {'adaID': 'TS6'},
      'mixem1968': {'adaID': 'TS12'},
     }
# Stampy parameters
stampy_gapopen = 60	        # Default: 40
stampy_gapextend = 5 	    # Default: 3
stampy_sensitive = True     # Default: False
subsrate = '0.05'



# Script
if __name__ == '__main__':

    # Parse input args: this is used to call itself recursively
    parser = argparse.ArgumentParser(description='Errors in HIV pure strains')
    parser.add_argument('--maxreads', type=int, default=100,
                        help='Number of reads analyzed')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--PCRtype', default='1542',
                        help='Mix to study (1542, 1968)')

    args = parser.parse_args()
    maxreads = args.maxreads
    VERBOSE = args.verbose
    PCRtype = args.PCRtype
    mixconvname = 'mixconv'+PCRtype
    mixemname = 'mixem'+PCRtype
    seqrun = 'Tuen6'

    dataset = MiSeq_runs[seqrun]
    data_folder = dataset['folder']

    for samplename in (mixconvname, mixemname):
        if VERBOSE >= 1:
            print samplename

        adaID = sd[samplename]['adaID']
        reads_filenames = get_read_filenames(data_folder, adaID, gzip=True)

        if maxreads != -1:
            read_subset_filenames = [r.replace('.fastq.gz', '_subset.fastq.gz') for r in reads_filenames]

            with gzip.open(reads_filenames[0], 'rb') as f1, \
                 gzip.open(reads_filenames[1], 'rb') as f2, \
                 gzip.open(read_subset_filenames[0], 'wb') as fo1, \
                 gzip.open(read_subset_filenames[1], 'wb') as fo2:

                if VERBOSE >= 1:
                    print 'Counting reads'
                #n_reads = get_number_reads_fastq_open(f1)
                n_reads = 500000
    
                if VERBOSE >= 1:
                    print 'Getting random indices'
                ind = np.arange(n_reads)
                np.random.shuffle(ind)
                ind = ind[:maxreads]
                ind.sort()

                ri1 = FGI(f1)
                ri2 = FGI(f2)

                if VERBOSE >= 1:
                    print 'Iterating over read pairs'

                ind_i = 0
                for irp, (read1, read2) in enumerate(izip(ri1, ri2)):

                    if irp == ind[ind_i]:
                        if VERBOSE >= 3:
                            print irp
                        elif VERBOSE == 2:
                            if not ((ind_i + 1) % 100):
                                print ind_i + 1

                        fo1.write("@%s\n%s\n+\n%s\n" % read1)
                        fo2.write("@%s\n%s\n+\n%s\n" % read2)
                        ind_i += 1
                        if ind_i == maxreads:
                            break

            reads_filenames = read_subset_filenames

            output_filename = data_folder+foldername_adapter(adaID)+'mapped.sam'
            # FIXME
            output_filename = output_filename.replace('.sam', '_subset.sam')

            # Map
            call_list = [stampy_bin,
                         '-g', data_folder+'emPCR_results/NL43_pol_emPCR_'+PCRtype,
                         '-h', data_folder+'emPCR_results/NL43_pol_emPCR_'+PCRtype,
                         '-o', output_filename,
                         '--overwrite',
                         '--substitutionrate='+subsrate,
                         '--gapopen', stampy_gapopen,
                         '--gapextend', stampy_gapextend,
                         '--sensitive']

            ## Take only a (random) subsample: stampy uses the fraction of reads
            ## intead of the number
            #if maxreads > 0:
            #    n_pairs_tot = get_number_reads(input_filename, 'bam') / 2
            #    frac_pairs = 1.0 * maxreads / n_pairs_tot
            #    random_seed = np.random.randint(1e5)
            #    call_list.append('-s', frac_pairs + random_seed)

            call_list = call_list + ['-M'] + reads_filenames
            call_list = map(str, call_list)
            if VERBOSE >=2:
                print ' '.join(call_list)
            sp.call(call_list)


