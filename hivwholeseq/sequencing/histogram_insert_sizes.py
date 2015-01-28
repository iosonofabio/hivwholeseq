# vim: fdm=indent
'''
author:     Fabio Zanini
date:       06/08/13
content:    Get the histogram of insert sizes.

            This is particularly important because some of the maps contain
            contaminants (fosmid plasmid chunks), which give rise to spuriously
            short inserts.
'''
# TODO: finish up!
# Modules
import os
import pysam
import numpy as np
from Bio import SeqIO

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.miseq import alpha
from hivwholeseq.sequencing.filenames import get_consensus_filename, get_mapped_filename
from hivwholeseq.utils.mapping import convert_sam_to_bam



# Globals
maxreads = 5000
match_len_min = 30
trim_bad_cigars = 3



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Demultiplex HIV reads')
    parser.add_argument('--run', required=True,
                        help='Seq run to analyze (e.g. Tue28)')
    parser.add_argument('--adaIDs', nargs='*',
                        help='Adapter IDs to analyze (e.g. TS2)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    seq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    VERBOSE = args.verbose

    # Specify the dataset
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = load_adapter_table(data_folder)['ID']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Iterate over all requested samples
    for adaID in adaIDs:
        for fragment in fragments:

            # Read reference
            reffilename = get_consensus_filename(data_folder, adaID, fragment)
            refseq = SeqIO.read(reffilename, 'fasta')
            ref = np.array(refseq)
        
            # Open BAM
            bamfilename = get_mapped_filename(data_folder, adaID, fragment,
                                              filtered=False)
            if not os.path.isfile(bamfilename):
                convert_sam_to_bam(bamfilename)
            with pysam.Samfile(bamfilename, 'rb') as bamfile:
        
                # Iterate through reads
                for i, read in enumerate(bamfile):
                
                    # Limit to the first reads
                    if i >= maxreads: break
            
                    # Print output
                    if VERBOSE and not ((i +1) % 10000):
                        print (i+1)
                
                    # Ignore unmapped reads
                    if read.is_unmapped or not read.is_proper_pair:
                        continue
                
                    # Divide by read 1/2 and forward/reverse
                    js = 2 * read.is_read2 + read.is_reverse
                
                    # Sequence and position
                    # Note: stampy takes the reverse complement already
                    seq = read.seq
                    pos = read.pos
                
                    # Read CIGAR code for indels, and anayze each block separately
                    # Note: some reads are weird ligations of HIV and contaminants
                    # (e.g. fosmid plasmids). Those always have crazy CIGARs, because
                    # only the HIV part maps. We hence trim away short indels from the
                    # end of reads (this is still unbiased).
                    cigar = read.cigar
                    len_cig = len(cigar)
                    good_cigars = np.array(map(lambda x: (x[0] == 0) and (x[1] >= match_len_min), cigar), bool, ndmin=1)
                    # If no long match, skip read
                    # FIXME: we must skip the mate pair also? But what if it's gone already?
                    # Note: they come in pairs: read1 first, read2 next, so we can just read two at a time
                    if not (good_cigars).any():
                        continue
                    # If only one long match, no problem
                    if (good_cigars).sum() == 1:
                        first_good_cigar = last_good_cigar = good_cigars.nonzero()[0][0]
                    # else include stuff between the extremes
                    else:
                        tmp = good_cigars.nonzero()[0]
                        first_good_cigar = tmp[0]
                        last_good_cigar = tmp[-1]
                        good_cigars[first_good_cigar: last_good_cigar + 1] = True
            
                    # Iterate over CIGARs
                    for ic, block in enumerate(cigar):
            
                        # Inline block
                        if block[0] == 0:
                            # Exclude bad CIGARs
                            if good_cigars[ic]: 
            
                                # The first and last good CIGARs are matches: trim them (unless they end the read)
                                if (ic == first_good_cigar) and (ic != 0): trim_left = trim_bad_cigars
                                else: trim_left = 0
                                if (ic == last_good_cigar) and (ic != len_cig - 1): trim_right = trim_bad_cigars
                                else: trim_right = 0
                
                                seqb = np.array(list(seq[trim_left:block[1] - trim_right]), 'S1')
                                # Increment counts
                                for j, a in enumerate(alpha):
                                    posa = (seqb == a).nonzero()[0]
                                    if len(posa):
                                        counts[js, j, pos + trim_left + posa] += 1
                
                            # Chop off this block
                            if ic != len_cig - 1:
                                seq = seq[block[1]:]
                                pos += block[1]
                
                        # Deletion
                        elif block[0] == 2:
                            # Exclude bad CIGARs
                            if good_cigars[ic]: 
                                # Increment gap counts
                                counts[js, 4, pos:pos + block[1]] += 1
                
                            # Chop off pos, but not sequence
                            pos += block[1]
                
                        # Insertion
                        elif block[0] == 1:
                            # Exclude bad CIGARs
                            if good_cigars[ic]: 
                                seqb = seq[:block[1]]
                                inserts[js][(pos, seqb)] += 1
                
                            # Chop off seq, but not pos
                            if ic != len_cig - 1:
                                seq = seq[block[1]:]
                
                        # Other types of cigar?
                        else:
                            raise ValueError('CIGAR type '+str(block[0])+' not recognized')
            
            
