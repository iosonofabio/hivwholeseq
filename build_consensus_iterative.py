#!/ebio/ag-neher/share/programs/EPD/bin/python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       12/08/13
content:    Map a subsample to HXB2 first, then build a consensus, then map
            recursively againt consensi until the latter do not change anymore.

            This elaborate scheme is necessary to cope with indel-rich regions
            (e.g. V loops) in variants away from HXB2 (e.g. F10).

            Note that mismapping is a serious problem for us, so maximal care
            must be taken here.

            Note: this script can call itself in a parallel fashion with the
            flag --submit as a cluster job.
'''
# Modules
import os
import sys
import subprocess as sp
import time
import argparse
import re
from operator import *
from itertools import izip
from collections import defaultdict
from collections import Counter
import numpy as np
import pysam
import Bio.SeqIO as SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Horizontal import of modules from this folder
from mapping.adapter_info import load_adapter_table, foldername_adapter
from mapping.miseq import alpha, read_types, pair_generator
from mapping.mapping_utils import stampy_bin, subsrate, convert_sam_to_bam
from mapping.filenames import get_HXB2_fragmented, get_read_filenames,\
        get_HXB2_index_file, get_HXB2_hash_file, get_consensus_filename


# Globals
VERBOSE = 3
# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']

# Consensus building
maxreads = 100000
match_len_min = 30
trim_bad_cigars = 3

# Cluster submit
import mapping
JOBDIR = mapping.__path__[0].rstrip('/')+'/'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
JOBSCRIPT = JOBDIR+'mapping_subsample_recursive.py'
cluster_time = '0:59:59'
vmem = '8G'
interval_check = 10



# Functions
def get_ref_file(data_folder, adaID, fragment, n_iter, ext=True):
    '''Get the reference filename'''
    if n_iter == 0:
        fn = get_HXB2_fragmented(data_folder, fragment=fragment)
        if not ext:
            # Strip the '.fasta' part
            fn = fn[:-6]
    else:
        fn = '_'.join(['consensus', str(n_iter-1), fragment])
        fn = data_folder+'subsample/'+foldername_adapter(adaID)+'map_iter/'+fn
        if ext:
            fn = fn+'.fasta'
    return fn


def get_index_file(data_folder, adaID, fragment, n_iter, ext=True):
    '''Get the index filename, with or w/o extension'''
    if n_iter == 0:
        return get_HXB2_index_file(data_folder, fragment, ext=ext)
    else:
        fn = get_ref_file(data_folder, adaID, fragment, n_iter, ext=False)
        if ext:
            fn = fn+'.stidx'
        return fn


def get_hash_file(data_folder, adaID, fragment, n_iter, ext=True):
    '''Get the index filename, with or w/o extension'''
    if n_iter == 0:
        return get_HXB2_hash_file(data_folder, fragment, ext=ext)
    else:
        fn = get_ref_file(data_folder, adaID, fragment, n_iter, ext=False)
        if ext:
            fn = fn+'.sthash'
        return fn


def get_mapped_filename(data_folder, adaID, fragment, n_iter, type='bam'):
    '''Get the mapped filenames'''
    filename = 'mapped_to_'
    if n_iter == 0:
        filename = filename + 'HXB2'
    else:
        filename = filename + 'consensus_'+str(n_iter - 1)
    filename = filename+'_'+fragment
    if type == 'sam':
        filename = filename + '.sam'
    elif type == 'bam':
        filename = filename + '.bam'
    else:
        raise ValueError('Type of mapped reads file not recognized')
 
    return data_folder+'subsample/'+foldername_adapter(adaID)+'map_iter/'+filename


def make_index_and_hash(data_folder, adaID, fragment, n_iter, VERBOSE=0):
    '''Make index and hash files for reference or consensus'''
    if VERBOSE:
        print 'Build stampy hashes: '+'{:02d}'.format(adaID)+' '+fragment\
                +' iteration '+str(n_iter)

    # Make folder if necessary
    dirname = os.path.dirname(get_index_file(data_folder, adaID, fragment, n_iter))
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    # 1. Make genome index file for 6 fragments (chromosomes)
    if not os.path.isfile(get_index_file(data_folder, adaID, fragment, n_iter, ext=True)):
        sp.call([stampy_bin,
                 '--species="HIV adaID '+'{:02d}'.format(adaID)+' fragment '+fragment+'"',
                 '-G', get_index_file(data_folder, adaID, fragment, n_iter, ext=False),
                 get_ref_file(data_folder, adaID, fragment, n_iter),
                 ])
    
    # 2. Build a hash file for 6 fragments
    if not os.path.isfile(get_hash_file(data_folder, adaID, fragment, n_iter, ext=True)):
        sp.call([stampy_bin,
                 '-g', get_index_file(data_folder, adaID, fragment, n_iter, ext=False),
                 '-H', get_hash_file(data_folder, adaID, fragment, n_iter, ext=False),
                 ])


def map_stampy(data_folder, adaID, fragment, n_iter, VERBOSE=0):
    '''Map using stampy'''
    if VERBOSE:
        print 'Map via stampy: '+'{:02d}'.format(adaID)+' '+fragment\
                +' iteration '+str(n_iter)

    read_filenames = get_read_filenames(data_folder, adaID, subsample=True,
                                        premapped=True, fragment=fragment)

    mapped_filename = get_mapped_filename(data_folder, adaID, fragment,
                                          n_iter, type='sam')

    # Make folder if necessary
    dirname = os.path.dirname(mapped_filename)
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    # Map
    call_list = [stampy_bin,
                 '-g', get_index_file(data_folder, adaID, fragment,
                                      n_iter, ext=False),
                 '-h', get_hash_file(data_folder, adaID, fragment,
                                     n_iter, ext=False), 
                 '-o', mapped_filename,
                 '--substitutionrate='+subsrate,
                 '-M'] + read_filenames
    if VERBOSE >=2:
        print ' '.join(call_list)
    sp.call(call_list)


def fork_self(data_folder, adaID, fragment, iterations_max=0, VERBOSE=0):
    '''Fork self for each adapter ID and fragment'''
    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'mit '+'{:02d}'.format(adaID)+' '+fragment,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--adaID', adaID,
                 '--fragment', fragment,
                 '--verbose', VERBOSE,
                 '--iterations', iterations_max,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def make_consensus(data_folder, adaID, fragment, n_iter, VERBOSE=0):
    '''Make consensus sequence from the mapped reads'''
    if VERBOSE:
        print 'Build consensus: '+'{:02d}'.format(adaID)+' '+fragment\
                +' iteration '+str(n_iter)
    
    # Read reference
    reffilename = get_ref_file(data_folder, adaID, fragment, n_iter)
    refseq = SeqIO.read(reffilename, 'fasta')
    ref = np.array(refseq)
    
    # Allele counts and inserts
    counts = np.zeros((len(read_types), len(alpha), len(ref)), int)
    inserts = [defaultdict(int) for k in read_types]

    # Open BAM file
    bamfilename = get_mapped_filename(data_folder, adaID, fragment, n_iter)
    if not os.path.isfile(bamfilename):
        convert_sam_to_bam(bamfilename)
    with pysam.Samfile(bamfilename, 'rb') as bamfile:

        # Iterate over pairs of reads
        for i_pairs, reads in enumerate(pair_generator(bamfile)):

            # Set read1 and read2
            (read1, read2) = reads

            # Check a few things to make sure we are looking at paired reads
            if read1.qname != read2.qname:
                raise ValueError('Read pair '+str(i_pairs)+': reads have different names!')

            # Ignore unmapped reads
            if (read1.is_unmapped or (not read1.is_proper_pair)
                or read2.is_unmapped or (not read2.is_proper_pair)):
                continue

            # Limit to the first reads
            if 2 * i_pairs >= maxreads: break

            # What to do with bad CIGARs? Shall we skip both?

            # Cover both reads of the pair
            for i_within_pair, read in enumerate(reads):
            
                # Divide by read 1/2 and forward/reverse
                if read.is_read1: js = 0
                else: js = 2
                if read.is_reverse: js += 1
            
                # FIXME: if the insert size is less than 250, we read over into the adapters
                # TODO: should not this be taken care of by our block-by-block strategy below?
                # Note: skip also the mate
                if read.isize <= 250:
                    pass
            
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
            
                # Sequence and position
                # Note: stampy takes the reverse complement already
                seq = read.seq
                pos = read.pos

                # Iterate over CIGARs
                for ic, (block_type, block_len) in enumerate(cigar):

                    # Check for pos: it should never exceed the length of the fragment
                    if (block_type in [0, 1, 2]) and (pos > len(ref)):
                        raise ValueError('Pos exceeded the length of the fragment')
            
                    # Inline block
                    if block_type == 0:
                        # Exclude bad CIGARs
                        if good_cigars[ic]: 
            
                            # The first and last good CIGARs are matches: trim them (unless they end the read)
                            if (ic == first_good_cigar) and (ic != 0): trim_left = trim_bad_cigars
                            else: trim_left = 0
                            if (ic == last_good_cigar) and (ic != len_cig - 1): trim_right = trim_bad_cigars
                            else: trim_right = 0
            
                            seqb = np.array(list(seq[trim_left:block_len - trim_right]), 'S1')
                            # Increment counts
                            for j, a in enumerate(alpha):
                                posa = (seqb == a).nonzero()[0]
                                if len(posa):
                                    counts[js, j, pos + trim_left + posa] += 1
            
                        # Chop off this block
                        if ic != len_cig - 1:
                            seq = seq[block_len:]
                            pos += block_len
            
                    # Deletion
                    elif block_type == 2:
                        # Exclude bad CIGARs
                        if good_cigars[ic]: 
                            # Increment gap counts
                            counts[js, 4, pos:pos + block_len] += 1
            
                        # Chop off pos, but not sequence
                        pos += block_len
            
                    # Insertion
                    # an insert @ pos 391 means that seq[:391] is BEFORE the insert,
                    # THEN the insert, FINALLY comes seq[391:]
                    elif block_type == 1:
                        # Exclude bad CIGARs
                        if good_cigars[ic]: 
                            seqb = seq[:block_len]
                            inserts[js][(pos, seqb)] += 1
            
                        # Chop off seq, but not pos
                        if ic != len_cig - 1:
                            seq = seq[block_len:]
            
                    # Other types of cigar?
                    else:
                        raise ValueError('CIGAR type '+str(block_type)+' not recognized')

    
    # Build consensus
    # Make consensi for each of the four categories
    consensi = np.zeros((len(read_types), len(refseq)), 'S1')
    for js, count in enumerate(counts):
        # Positions without reads are considered N
        # (this should happen only at the ends)
        count.T[(count).sum(axis=0) == 0] = np.array([0, 0, 0, 0, 0, 1])
        consensi[js] = alpha[count.argmax(axis=0)]
    
    # Add inserts that are present in more than half the reads
    inserts_consensi = []
    for js, insert in enumerate(inserts):
        insert_consensus = []
        for (k, i) in insert.iteritems():
            pos_pre = max(0, k[0]-1)
            pos_post = min(len(refseq), k[0]+1)
            if pos_post <= pos_pre:
                print 'Strange insert:', pos_pre, pos_post, k, len(refseq), counts[js, :, pos_pre:pos_post].sum(axis=0)
            threshold = 0.5 * counts[js, :, pos_pre:pos_post].sum(axis=0).mean()
            if (i > 10) and (i > threshold):
                insert_consensus.append(k)
        inserts_consensi.append(insert_consensus)
    
    # Make final consensus
    # This has two steps: 1. counts; 2. inserts
    # We should check that different read types agree (e.g. forward/reverse) on
    # both counts and inserts if they are covered
    # 1.0: Prepare everything as 'N': ambiguous
    consensus = np.repeat('N', len(refseq))
    # 1.1: Put safe stuff
    ind_agree = (consensi == consensi[0]).all(axis=0)
    consensus[ind_agree] = consensi[0, ind_agree]
    # 1.2: Ambiguous stuff requires more care
    polymorphic = []
    # 1.2.1: If is covered by some read types only, look among those
    amb_pos = (-ind_agree).nonzero()[0]
    for pos in amb_pos:
        cons_pos = consensi[:, pos]
        cons_pos = cons_pos[cons_pos != 'N']
        # If there is unanimity, ok
        if (cons_pos == cons_pos[0]).all():
            consensus[pos] = cons_pos[0]
        # else, assign only likely mismapped deletions
        else:
            # Restrict to nongapped things (caused by bad mapping)
            cons_pos = cons_pos[cons_pos != '-']
            if (cons_pos == cons_pos[0]).all():
                consensus[pos] = cons_pos[0]
            # In case of polymorphisms, take any most abundant nucleotide
            else:
                polymorphic.append(pos)
                tmp = zip(*Counter(cons_pos).items())
                consensus[pos] = tmp[0][np.argmax(tmp[1])]
    
    # 2.0 get all inserts
    insert_names = set()
    for insert in inserts:
        insert_names |= set(insert.keys())
    # 2.1 iterate over inserts and check all read types
    insert_consensus = []
    for insert_name in insert_names:
        # Get counts for the insert and local coverage
        ins_counts = np.array([insert[insert_name] for insert in inserts])
        # FIXME: is there a division by zero here?
        cov_loc = counts[:, :, insert_name[0]: min(len(refseq), insert_name[0] + 2)].sum(axis=1).mean(axis=1)
        # Get read types in which the insert is called
        types_ins_called = ins_counts > 0.5 * cov_loc
        # Get read types which actually cover this region
        types_cov = cov_loc > 3
        # Check the insert counts compared to the local coverage
        # if any read type has coverage and all read types with coverage agree, ok
        if types_cov.sum() and (types_ins_called[types_cov]).all():
            insert_consensus.append(insert_name)
    insert_consensus.sort(key=itemgetter(0))
    # 3. put inserts in
    consensus_final = []
    pos = 0
    for insert_name in insert_consensus:
        # Indices should be fine...
        # an insert @ pos 391 means that seq[:391] is BEFORE the insert,
        # THEN the insert, FINALLY comes seq[391:]
        consensus_final.append(''.join(consensus[pos:insert_name[0]]))
        consensus_final.append(insert_name[1])
        pos = insert_name[0]
    consensus_final.append(''.join(consensus[pos:]))
    # Strip initial and final Ns and gaps
    consensus_final = ''.join(consensus_final).strip('N')
    consensus_final = re.sub('-', '', consensus_final)

    return refseq, consensus_final


def check_new_old_consensi(refseq, consensus):
    '''Check the old and the new consensi, whether they match'''
    # exact match
    match = str(refseq.seq) == consensus
    return match


def write_consensus(data_folder, adaID, fragment, n_iter, consensus, final=False):
    '''Write consensus sequences into FASTA'''
    name = 'adaID_'+'{:02d}'.format(adaID)+'_'+fragment+'_consensus'
    consensusseq = SeqRecord(Seq(consensus),
                             id=name, name=name)
    if final:
        outfile = get_consensus_filename(data_folder, adaID, fragment)
    else:
        outfile = get_ref_file(data_folder, adaID, fragment, n_iter+1)
    SeqIO.write(consensusseq, outfile, 'fasta')
             




# Script
if __name__ == '__main__':

    # Parse input args: this is used to call itself recursively
    parser = argparse.ArgumentParser(description='Map HIV reads recursively')
    parser.add_argument('--adaID', metavar='00', type=int, nargs='?',
                        default=0,
                        help='Adapter ID sample to map')
    parser.add_argument('--fragment', metavar='F0', nargs='?',
                        default='F0',
                        help='Fragment to map (F1-F6)')
    parser.add_argument('--iterations', type=int, default=0,
                        help=('Maximal number of map/consensus iterations'))
    parser.add_argument('--verbose', type=int, default=0,
                        help=('Verbosity level [0-3]'))
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')

    args = parser.parse_args()
    adaID = args.adaID
    fragment = args.fragment
    iterations_max = args.iterations
    VERBOSE = args.verbose
    submit = args.submit

    # If the script is called with no adaID, iterate over all
    if adaID == 0:
        adaIDs = load_adapter_table(data_folder)['ID']
    else:
        adaIDs = [adaID]
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # If the script is called with no fragment, iterate over all
    if fragment == 'F0':
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    else:
        fragments = [fragment]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Iterate over all requested samples
    for frag in fragments:

        # Make hashes for the reference, they are shared across adaIDs
        make_index_and_hash(data_folder, adaID, frag, 0, VERBOSE=VERBOSE)

        for adaID in adaIDs:

            # Submit to the cluster self if requested
            if submit:
                fork_self(data_folder, adaID, frag, iterations_max, VERBOSE=VERBOSE)

            # or else, build the consensus
            else:

                # Iterate the consensus building until convergence
                n_iter = 0
                done = False
                while not done:
                
                    # Make hashes of consensus (the reference is hashed earlier)
                    if (n_iter > 0):
                        make_index_and_hash(data_folder, adaID, frag, n_iter, VERBOSE=VERBOSE)

                    # Map against reference or consensus
                    map_stampy(data_folder, adaID, frag, n_iter, VERBOSE=VERBOSE)

                    # Build consensus
                    refseq, consensus = make_consensus(data_folder, adaID, fragment, n_iter, VERBOSE=VERBOSE)

                    # Save consensus to file
                    write_consensus(data_folder, adaID, fragment, n_iter, consensus, final=False)

                    # Check whether the old and new consensi match
                    match = check_new_old_consensi(refseq, consensus)

                    # Start a new round
                    if (not match) and (n_iter + 1 < iterations_max):
                            n_iter += 1
                            if VERBOSE:
                                print 'Starting again for iteration '+str(n_iter) 

                    # or terminate
                    else:
                        write_consensus(data_folder, adaID, fragment, n_iter, consensus, final=True)
                        done = True

                        if VERBOSE:
                            if match:
                                print 'Consensus converged.' 
                            else:
                                print 'Maximal number of iterations reached.'
                        
