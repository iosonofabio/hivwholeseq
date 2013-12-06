#!/usr/bin/env python
# vim: fdm=indent
'''
author:     Fabio Zanini
date:       12/08/13
content:    Take reads mapped to HXB2, then build a consensus, then map
            recursively againt consensi until the latter do not change anymore.

            This elaborate scheme is necessary to cope with indel-rich regions
            (e.g. V loops) in variants away from HXB2 (e.g. F10).

            Note that mismapping is a serious problem for us, so maximal care
            must be taken here.
'''
#TODO: tidy up the F5a/b mess.
# Modules
import os
import subprocess as sp
import argparse
import re
from operator import itemgetter
from collections import defaultdict
from collections import Counter
import numpy as np
import pysam
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Horizontal import of modules from this folder
from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.adapter_info import load_adapter_table, foldername_adapter
from hivwholeseq.miseq import alpha, read_types
from hivwholeseq.mapping_utils import stampy_bin, subsrate, convert_sam_to_bam,\
        pair_generator, align_muscle
from hivwholeseq.filenames import get_HXB2_fragmented, \
        get_HXB2_index_file, get_HXB2_hash_file, get_consensus_filename, \
        get_divided_filenames
from hivwholeseq.one_site_statistics import get_allele_counts_insertions_from_file_unfiltered


# Globals
# Stampy parameters
stampy_gapopen = 60	        # Default: 40
stampy_gapextend = 5 	    # Default: 3
stampy_sensitive = True     # Default: False

# Consensus building
maxreads = 1e5
match_len_min = 10
trim_bad_cigars = 3
coverage_min = 10
# Minimal quality required for a base to be considered trustful (i.e. added to 
# the allele counts), in phred score. Too high: lose coverage, too low: seq errors.
# Reasonable numbers are between 30 and 36.
qual_min = 30

# Cluster submit
import hivwholeseq
JOBDIR = hivwholeseq.__path__[0].rstrip('/')+'/'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
JOBSCRIPT = JOBDIR+'build_consensus_iterative.py'
cluster_time = '0:59:59'
vmem = '2G'
interval_check = 10



# Functions
def fork_self(miseq_run, adaID, fragment, n_reads=1000, iterations_max=0, VERBOSE=0):
    '''Fork self for each adapter ID and fragment'''
    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'bci '+'{:02d}'.format(adaID)+' '+fragment,
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--run', miseq_run,
                 '--adaIDs', adaID,
                 '--fragments', fragment,
                 '-n', n_reads,
                 '--iterations', iterations_max,
                 '--verbose', VERBOSE,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def get_ref_file(data_folder, adaID, fragment, n_iter, ext=True):
    '''Get the reference filename'''
    if n_iter == 0:
        fn = get_HXB2_fragmented(fragment, ext=ext, trim_primers=True)
    else:
        # There are two sets of primers for fragment F5
        if 'F5' in fragment:
            fragment = 'F5'
        fn = '_'.join(['consensus', str(n_iter-1), fragment])
        fn = data_folder+foldername_adapter(adaID)+'map_iter/'+fn
        if ext:
            fn = fn+'.fasta'
    return fn


def get_refcum_file(data_folder, adaID, fragment, ext=True):
    '''Get the file with the cumulated consensi'''
    # There are two sets of primers for fragment F5
    if 'F5' in fragment:
        fragment = 'F5'
    fn = '_'.join(['consensus', 'alliters', fragment])
    fn = data_folder+foldername_adapter(adaID)+'map_iter/'+fn
    if ext:
        fn = fn+'.fasta'
    return fn


def get_index_file(data_folder, adaID, fragment, n_iter, ext=True):
    '''Get the index filename, with or w/o extension'''
    if n_iter == 0:
        return get_HXB2_index_file(fragment, ext=ext)
    else:
        # There are two sets of primers for fragment F5
        if 'F5' in fragment:
            fragment = 'F5'
        fn = get_ref_file(data_folder, adaID, fragment, n_iter, ext=False)
        if ext:
            fn = fn+'.stidx'
        return fn


def get_hash_file(data_folder, adaID, fragment, n_iter, ext=True):
    '''Get the index filename, with or w/o extension'''
    if n_iter == 0:
        return get_HXB2_hash_file(fragment, ext=ext)
    else:
        # There are two sets of primers for fragment F5
        if 'F5' in fragment:
            fragment = 'F5'
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
 
    return data_folder+foldername_adapter(adaID)+'map_iter/'+filename


def extract_reads_subsample(data_folder, adaID, fragment, n_reads, VERBOSE=0):
    '''Extract a subsample of reads from the initial sample mapped to HXB2'''
    # There are two sets of primers for fragment F5
    if 'F5' in fragment:
        frag_out = 'F5'
    else:
        frag_out = fragment

    # Count the number of reads first
    input_filename = get_divided_filenames(data_folder, adaID, [fragment], type='bam')[0]
    with pysam.Samfile(input_filename, 'rb') as bamfile_in:
        n_reads_tot = sum(1 for read in bamfile_in) / 2

    # Pick random numbers among those
    # Get the random indices of the reads to store
    # (only from the middle of the file)
    ind_store = np.arange(int(0.00 * n_reads_tot), int(1 * n_reads_tot))
    np.random.shuffle(ind_store)
    ind_store = ind_store[:n_reads]
    ind_store.sort()

    if VERBOSE >= 2:
        print 'Random indices between '+str(ind_store[0])+' and '+str(ind_store[-1])

    # Output file
    output_filename = get_mapped_filename(data_folder, adaID, frag_out, 0, type='bam')
    # Make folder if necessary
    dirname = os.path.dirname(output_filename)
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    # Copy those to subsample
    with pysam.Samfile(input_filename, 'rb') as bamfile_in:
        with pysam.Samfile(output_filename, 'wb', template=bamfile_in) as bamfile_out:

            n_written = 0
            for i, (read1, read2) in enumerate(pair_generator(bamfile_in)):
                # Log activity
                if VERBOSE >= 3:
                    if not ((i+1) % 10000):
                        print i+1, n_written, ind_store[n_written]
    
                # If you hit a read, write them
                if i == ind_store[n_written]:
                    bamfile_out.write(read1)
                    bamfile_out.write(read2)
                    n_written += 1
    
                # Break after the last one
                if n_written >= n_reads:
                    break


def make_index_and_hash(data_folder, adaID, fragment, n_iter, VERBOSE=0):
    '''Make index and hash files for reference or consensus'''
    # There are two sets of primers for fragment F5
    if (n_iter > 0) and ('F5' in fragment):
         fragment = 'F5'

    if VERBOSE:
        print 'Build stampy hashes: '+'{:02d}'.format(adaID)+' '+fragment\
                +' iteration '+str(n_iter)

    # Make folder if necessary
    dirname = os.path.dirname(get_index_file(data_folder, adaID, fragment, n_iter))
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    # 1. Make genome index file for 6 fragments (chromosomes)
    if not os.path.isfile(get_index_file(data_folder, adaID, fragment, n_iter, ext=True)):
        call_list = [stampy_bin,
                     '--species="HIV adaID '+'{:02d}'.format(adaID)+' fragment '+fragment+'"',
                     '-G', get_index_file(data_folder, adaID, fragment, n_iter, ext=False),
                     get_ref_file(data_folder, adaID, fragment, n_iter),
                    ]
        if VERBOSE >= 3:
            print ' '.join(call_list)
        sp.call(call_list)
    
    # 2. Build a hash file for 6 fragments
    if not os.path.isfile(get_hash_file(data_folder, adaID, fragment, n_iter, ext=True)):
        call_list = [stampy_bin,
                     '-g', get_index_file(data_folder, adaID, fragment, n_iter, ext=False),
                     '-H', get_hash_file(data_folder, adaID, fragment, n_iter, ext=False),
                    ]
        if VERBOSE >= 3:
            print ' '.join(call_list)
        sp.call(call_list)


def map_stampy(data_folder, adaID, fragment, n_iter, VERBOSE=0):
    '''Map using stampy'''
    # There are two sets of primers for fragment F5
    if 'F5' in fragment:
        frag_out = 'F5'
    else:
        frag_out = fragment

    if VERBOSE:
        print 'Map via stampy: '+'{:02d}'.format(adaID)+' '+fragment\
                +' iteration '+str(n_iter)

    # Input and output files
    input_filename = get_mapped_filename(data_folder, adaID, frag_out,
                                          n_iter - 1, type='bam')
    output_filename = get_mapped_filename(data_folder, adaID, frag_out,
                                          n_iter, type='sam')
    # Make folder if necessary
    dirname = os.path.dirname(output_filename)
    if not os.path.isdir(dirname):
        os.mkdir(dirname)

    # Map
    call_list = [stampy_bin,
                 '-g', get_index_file(data_folder, adaID, fragment,
                                      n_iter, ext=False),
                 '-h', get_hash_file(data_folder, adaID, fragment,
                                     n_iter, ext=False), 
                 '-o', output_filename,
                 '--substitutionrate='+subsrate,
                 '--gapopen', stampy_gapopen,
                 '--gapextend', stampy_gapextend]
    if stampy_sensitive:
        call_list.append('--sensitive')
    call_list = call_list + ['-M', input_filename]
    call_list = map(str, call_list)
    if VERBOSE >=2:
        print ' '.join(call_list)
    sp.call(call_list)


def make_consensus(data_folder, adaID, fragment, n_iter, qual_min=20, VERBOSE=0):
    '''Make consensus sequence from the mapped reads'''
    # There are two sets of primers for fragment F5
    if 'F5' in fragment:
        frag_out = 'F5'
    else:
        frag_out = fragment

    if VERBOSE:
        print 'Build consensus: '+'{:02d}'.format(adaID)+' '+frag_out\
                +' iteration '+str(n_iter)
    
    # Read reference
    reffilename = get_ref_file(data_folder, adaID, fragment, n_iter)
    refseq = SeqIO.read(reffilename, 'fasta')
    ref = np.array(refseq)

    # Open BAM file
    bamfilename = get_mapped_filename(data_folder, adaID, frag_out, n_iter)
    if not os.path.isfile(bamfilename):
        convert_sam_to_bam(bamfilename)

    # Call lower-level function for allele counts and inserts
    (counts, inserts) = get_allele_counts_insertions_from_file_unfiltered(bamfilename,\
                                len(refseq), qual_min=qual_min,
                                match_len_min=match_len_min)

    # Build consensus
    # Make allele count consensi for each of the four categories
    consensi = np.zeros((len(read_types), len(refseq)), 'S1')
    for js, count in enumerate(counts):
        # Positions without reads are considered N
        # (this should happen only at the ends)
        count.T[(count).sum(axis=0) == 0] = np.array([0, 0, 0, 0, 0, 1])
        consensi[js] = alpha[count.argmax(axis=0)]

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
    
    # 2. Inserts
    # This is a bit tricky, because we could have a variety of inserts at the
    # same position because of PCR or sequencing errors (what we do with that
    # stuff, is another question). So, assuming the mapping program maps
    # everything decently at the same starting position, we have to look for
    # prefix families
    insert_consensus = []

    # Iterate over all insert positions
    for pos, insert in inserts.iteritems():

        # Get the coverage around the insert in all four read types
        # Note: the sum is over the four nucleotides, the mean over the two positions
        cov = counts[:, :, pos: min(len(refseq), pos + 2)].sum(axis=1).mean(axis=1)

        # Determine what read types have coverage
        is_cov = cov > coverage_min
        covcs = cov[is_cov]

        # A single covered read type is sufficient: we are going to do little
        # with that, but as a consensus it's fine
        if not is_cov.any():
            continue

        # Use prefixes/suffixes, and come down from highest frequencies
        # (averaging fractions is not a great idea, but -- oh, well).
        # Implement it as a count table: this is not very efficient a priori,
        # but Python is lame anyway if we start making nested loops
        len_max = max(map(len, insert.keys()))
        align_pre = np.tile('-', (len(insert), len_max))
        align_suf = np.tile('-', (len(insert), len_max))
        counts_loc = np.zeros((len(insert), len(read_types)), int)
        for i, (s, cs) in enumerate(insert.iteritems()):
            align_pre[i, :len(s)] = list(s)
            align_suf[i, -len(s):] = list(s)
            counts_loc[i] = cs

        # Make the allele count tables
        counts_pre_table = np.zeros((len(read_types), len(alpha), len_max), int)
        counts_suf_table = np.zeros((len(read_types), len(alpha), len_max), int)
        for j, a in enumerate(alpha):
            counts_pre_table[:, j] = np.dot(counts_loc.T, (align_pre == a))
            counts_suf_table[:, j] = np.dot(counts_loc.T, (align_suf == a))

        # Look whether any position of the insertion has more than 50% freq
        # The prefix goes: ----> x
        ins_pre = []
        for cs in counts_pre_table.swapaxes(0, 2):
            freq = 1.0 * (cs[:, is_cov] / covcs).mean(axis=1)
            iaM = freq.argmax()
            if (alpha[iaM] != '-') and (freq[iaM] > 0.5):
                ins_pre.append(alpha[iaM])
            else:
                break
        # The suffix goes: x <----
        ins_suf = []
        for cs in counts_suf_table.swapaxes(0, 2)[::-1]:
            freq = 1.0 * (cs[:, is_cov] / covcs).mean(axis=1)
            iaM = freq.argmax()
            if (alpha[iaM] != '-') and (freq[iaM] > 0.5):
                ins_suf.append(alpha[iaM])
            else:
                break
        ins_suf.reverse()

        if VERBOSE >= 4:
            if ins_pre or ins_suf:
                print ''.join(ins_pre)
                print ''.join(ins_suf)

        # Add the insertion to the list (they should agree in most cases)
        if ins_pre:
            insert_consensus.append((pos, ''.join(ins_pre)))
        elif ins_suf:
            insert_consensus.append((pos, ''.join(ins_suf)))

    # 3. put inserts in
    insert_consensus.sort(key=itemgetter(0))
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


    ## FIXME
    #pos = 979
    #print n_iter, pos
    #print counts[:, :, pos]
    #print consensus[979]


    ## FIXME
    #motif = 'TGCCTTGGAA'
    #pos = consensi[0].tostring().find(motif) + len(motif)
    #print n_iter
    #print motif
    #print consensi[0].tostring()[pos-len(motif): pos + 1]
    #print counts[:, :, pos]
    #import ipdb; ipdb.set_trace()

    return refseq, consensus_final


def check_new_old_consensi(refseq, consensus):
    '''Check the old and the new consensi, whether they match'''
    # exact match
    match = str(refseq.seq) == consensus
    return match


def write_consensus_intermediate(data_folder, adaID, fragment, n_iter, consensus):
    '''Write consensus sequences into FASTA'''
    # There are two sets of primers for fragment F5
    if 'F5' in fragment:
        frag_out = 'F5'
    else:
        frag_out = fragment

    name = 'adaID_'+'{:02d}'.format(adaID)+'_'+frag_out+'_consensus'
    consensusseq = SeqRecord(Seq(consensus),
                             id=name, name=name)
    outfile = get_ref_file(data_folder, adaID, frag_out, n_iter+1)
    SeqIO.write(consensusseq, outfile, 'fasta')

    # In addition, cumulate all consensi in one file for checking purposes,
    outfile = get_refcum_file(data_folder, adaID, frag_out)
    # If it's the first iteration, write the reference too
    if n_iter == 1:
        SeqIO.write(SeqIO.read(get_ref_file(data_folder, adaID, fragment, 0), 'fasta'),
                    outfile, 'fasta')
    consensusseq.name = consensusseq.id = consensusseq.name+'_'+str(n_iter)
    with open(outfile, 'a') as f:
        SeqIO.write(consensusseq, f, 'fasta')


def write_consensus_final(data_folder, adaID, fragment, consensus):
    '''Write the final consensus'''
    # There are two sets of primers for fragment F5
    if 'F5' in fragment:
        frag_out = 'F5'
    else:
        frag_out = fragment

    # Get the consensus
    name = 'adaID_'+'{:02d}'.format(adaID)+'_'+frag_out+'_consensus'
    consensusseq = SeqRecord(Seq(consensus), id=name, name=name)

    # Output file
    outfile = get_consensus_filename(data_folder, adaID, frag_out, trim_primers=True)

    # Write output
    SeqIO.write(consensusseq, outfile, 'fasta')

    # Align all consensi via muscle and store
    seqs = list(SeqIO.parse(get_refcum_file(data_folder, adaID, frag_out), 'fasta'))
    ali = align_muscle(*seqs)
    AlignIO.write(ali, get_refcum_file(data_folder, adaID, frag_out), 'fasta')




# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Build consensus, iteratively')
    parser.add_argument('--run', type=int, required=True,
                        help='MiSeq run to analyze (e.g. 28, 37)')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('-n', type=int, default=1000,
                        help='Number of (random) reads used for the consensus')
    parser.add_argument('--iterations', type=int, default=5,
                        help=('Maximal number of map/consensus iterations'))
    parser.add_argument('--verbose', type=int, default=0,
                        help=('Verbosity level [0-3]'))
    parser.add_argument('--submit', action='store_true',
                        help='Execute the script in parallel on the cluster')

    args = parser.parse_args()
    miseq_run = args.run
    adaIDs = args.adaIDs
    fragments = args.fragments
    n_reads = args.n
    iterations_max = args.iterations
    VERBOSE = args.verbose
    submit = args.submit

    # Specify the dataset
    dataset = MiSeq_runs[miseq_run]
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

            # There are two sets of primers for fragment F5
            if fragment == 'F5':
                fragment = dataset['primerF5'][dataset['adapters'].index(adaID)]

            # Submit to the cluster self if requested
            if submit:
                fork_self(miseq_run, adaID, fragment, n_reads, iterations_max,
                          VERBOSE=VERBOSE)
                continue

            # Extract subsample of reads mapped to HXB2
            extract_reads_subsample(data_folder, adaID, fragment, n_reads,
                                    VERBOSE=VERBOSE)
            # Build first consensus and save to file
            refseq, consensus = make_consensus(data_folder, adaID, fragment, 0,
                                               VERBOSE=VERBOSE)
            write_consensus_intermediate(data_folder, adaID, fragment, 0, consensus)

            # NOTE: we assume the consensus is different from HXB2, so we do not
            # check any further and do at least one more iteration

            # Iterate the consensus building until convergence
            n_iter = 1
            while True:
            
                # Make hashes of consensus (the reference is hashed earlier)
                make_index_and_hash(data_folder, adaID, fragment, n_iter,
                                    VERBOSE=VERBOSE)

                # Map against reference or consensus
                map_stampy(data_folder, adaID, fragment, n_iter, VERBOSE=VERBOSE)

                # Build consensus
                refseq, consensus = make_consensus(data_folder, adaID, fragment,
                                                   n_iter, VERBOSE=VERBOSE)

                # Save consensus to file
                write_consensus_intermediate(data_folder, adaID, fragment, n_iter, consensus)

                # Check whether the old and new consensi match
                match = check_new_old_consensi(refseq, consensus)

                # Start a new round
                if (not match) and (n_iter + 1 < iterations_max):
                    n_iter += 1
                    if VERBOSE:
                        print 'Starting again for iteration '+str(n_iter) 

                # or terminate
                else:
                    write_consensus_final(data_folder, adaID, fragment, consensus)
                    if VERBOSE:
                        if match:
                            print 'Consensus converged at iteration '+str(n_iter)+\
                                    ': adaID', adaID, fragment 
                        else:
                            print 'Maximal number of iterations reached: adaID', \
                                    adaID, fragment
                    
                    break 
