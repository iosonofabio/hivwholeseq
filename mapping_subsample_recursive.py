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

            Note: this script forks in a tree-like fashion. The first call might
            refer to a single adapter ID, in which case the script controls one
            single child, the stampy process, and waits until that's done. This
            script might also get called for all adapter IDs, in which case it
            proceeds linearly to produce hash files of HXB2, then forks into
            subprocesses, each of which tends to its own stampy child. This is
            not most efficient in terms of cluster jobs, but much easier to
            control and debug than a single master script (unless we wait for all
            IDs to finish at each iteration).
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
import mapping
from mapping.adapter_info import load_adapter_table, foldername_adapter
from mapping.miseq import alpha, read_types, pair_generator


# Globals
VERBOSE = 3
# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']
HXB2_fragmented_file = 'HXB2_fragmented.fasta'

# Mapping
stampy_bin = '/ebio/ag-neher/share/programs/bundles/stampy-1.0.22/stampy.py'
subsrate = '0.05'

# Consensus building
maxreads = 100000
match_len_min = 30
trim_bad_cigars = 3

# Cluster submit
JOBDIR = mapping.__path__[0].rstrip('/')+'/'
JOBLOGERR = JOBDIR+'logerr'
JOBLOGOUT = JOBDIR+'logout'
JOBSCRIPT = JOBDIR+'mapping_subsample_recursive.py'
cluster_time = '0:59:59'
vmem = '8G'
interval_check = 10



# Functions
def get_ref_file(data_folder, adaID=0, refname='HXB2'):
    '''Get the reference filename'''
    if refname == 'HXB2':
        reffile = data_folder+HXB2_fragmented_file
    elif adaID == 0:
        raise ValueError('Adapter ID not specified.')
    else:
        dirname = foldername_adapter(adaID)
        reffile = data_folder+'subsample/'+dirname+refname+'_fragmented.fasta'
    return reffile


def get_index_file(data_folder, adaID=0, refname='HXB2', ext=True):
    '''Get the index filename, with or w/o extension'''
    if refname == 'HXB2':
        fn = data_folder+'subsample/HXB2'
    else:
        dirname = foldername_adapter(adaID)
        fn = data_folder+'subsample/'+dirname+refname+'_fragmented'
    if ext:
        fn = fn+'.stidx'
    return fn


def get_hash_file(data_folder, adaID=0, refname='HXB2', ext=True):
    '''Get the index filename, with or w/o extension'''
    if refname == 'HXB2':
        fn = data_folder+'subsample/HXB2'
    else:
        dirname = foldername_adapter(adaID)
        fn = data_folder+'subsample/'+dirname+refname+'_fragmented'
    if ext:
        fn = fn+'.sthash'
    return fn


def get_mapped_file(data_folder, adaID, refname):
    '''Get the mapped filename'''
    dirname = foldername_adapter(adaID)
    return data_folder+'subsample/'+dirname+'mapped_to_'+refname+'_fragmented.sam'


def get_sam_file(data_folder, adaID, refname):
    '''Get the SAM filename'''
    dirname = foldername_adapter(adaID)
    return data_folder+'subsample/'+dirname+'mapped_to_'+refname+'_fragmented.sam'


def get_bam_file(data_folder, adaID, refname):
    '''Get the SAM filename'''
    dirname = foldername_adapter(adaID)
    return data_folder+'subsample/'+dirname+'mapped_to_'+refname+'_fragmented.bam'


def get_iteration(refname):
    '''Find the iteration number'''
    if refname == 'HXB2':
        return -1
    try:
        refnum = int(refname.split('_')[-1])
    except ValueError:
        raise ValueError('Reference not understood')
    return refnum


def transform_reference(reference):
    '''Add some text to the reference without affecting the command line'''
    if reference == 'HXB2':
        return reference

    if (len(reference) < 2) or (reference[0] != 'c'):
        raise ValueError('Reference not understood')
    try:
        refnum = int(reference[1:])
    except ValueError:
        raise ValueError('Reference not understood')

    return 'consensus_'+str(refnum)


def detransform_reference(refname):
    '''Delete some test from reference'''
    if refname == 'HXB2':
        return refname
    refnum = get_iteration(refname)
    return 'c'+str(refnum)


def next_refname(refname, shortversion=False):
    '''Find the next reference'''
    refnum = get_iteration(refname)
    nextref = 'c'+str(refnum+1)
    if shortversion:
        return nextref
    else:
        return transform_reference(nextref)


def make_index_and_hash(data_folder, refname, adaID=0):
    '''Make index and hash files for reference or consensus'''
    # 1. Make genome index file for 6 fragments (chromosomes)
    if not os.path.isfile(get_index_file(data_folder, adaID, refname, ext=True)):
        if VERBOSE:
            print 'Make '+refname+' index'
        sp.call([stampy_bin,
                 '--species="HIV 6 fragments"',
                 '-G', get_index_file(data_folder, adaID, refname, ext=False),
                 get_ref_file(data_folder, adaID, refname),
                 ])
    
    # 2. Build a hash file for 6 fragments
    if not os.path.isfile(get_hash_file(data_folder, adaID, refname, ext=True)):
        if VERBOSE:
            print 'Make '+refname+' hash'
        sp.call([stampy_bin,
                 '-g', get_index_file(data_folder, adaID, refname, ext=False),
                 '-H', get_hash_file(data_folder, adaID, refname, ext=False),
                 ])


def call_stampy_for_mapping(data_folder, refname, adaID):
    '''Call stampy for actual mapping'''
    # Directory to read
    dirname = foldername_adapter(adaID)
    if not os.path.isdir(data_folder+'subsample/'+dirname):
        raise ValueError('Adapter ID '+str(adaID)+': folder not found')
    
    if VERBOSE:
        print 'Map adaID '+str(adaID)+' to '+refname
    
    # Stampy command line
    # Note: no --solexa option as of 2013 (illumina v1.8)
    qsub_list = ['qsub','-cwd',
                 '-o', JOBLOGOUT,
                 '-e', JOBLOGERR,
                 '-N', 'stampy_'+'{:02d}'.format(adaID),
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 stampy_bin,
                 '-g', get_index_file(data_folder, adaID, refname, ext=False),
                 '-h', get_hash_file(data_folder, adaID, refname, ext=False), 
                 '-o', get_mapped_file(data_folder, adaID, refname),
                 '--substitutionrate='+subsrate,
                 '-M',
                 data_folder+'subsample/'+dirname+'read1_filtered_trimmed.fastq',
                 data_folder+'subsample/'+dirname+'read2_filtered_trimmed.fastq']
    qsub_list = map(str,qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)

    # Call stampy and check output
    qsub_output = sp.check_output(qsub_list).rstrip('\n')
    if VERBOSE:
        print qsub_output
    jobID = qsub_output.split(' ')[2]
    if VERBOSE:
        print jobID

    return jobID


def fork_self(data_folder, refname, adaID, stage=3, iterations_max=0):
    '''Fork self for each adapter ID'''
    qsub_list = ['qsub','-cwd',
                 '-b', 'y',
                 '-S', '/bin/bash',
                 '-o',JOBLOGOUT,
                 '-e',JOBLOGERR,
                 '-N', 'map_'+'{:02d}'.format(adaID),
                 '-l', 'h_rt='+cluster_time,
                 '-l', 'h_vmem='+vmem,
                 JOBSCRIPT,
                 '--adaID', adaID,
                '--reference', detransform_reference(refname),
                 '--stage', stage,
                 '--verbose', 3,
                 '--iterations', iterations_max,
                ]
    qsub_list = map(str, qsub_list)
    if VERBOSE:
        print ' '.join(qsub_list)
    sp.call(qsub_list)


def make_consensus(data_folder, refname, adaID):
    '''Make consensus sequence from the mapped reads'''
    # Directory to read
    dirname = foldername_adapter(adaID)
    if not os.path.isdir(data_folder+'subsample/'+dirname):
        raise ValueError('Adapter ID '+str(adaID)+': folder not found')
    
    if VERBOSE:
        print 'Building consensus after mapping against '+refname
    
    # Convert SAM (output of stampy) to BAM (more efficient)
    if not os.path.isfile(get_bam_file(data_folder, adaID, refname)):
        if VERBOSE:
            print 'Converting SAM to BAM...',
        samfile = pysam.Samfile(get_sam_file(data_folder, adaID, refname), 'r')
        bamfile = pysam.Samfile(get_bam_file(data_folder, adaID, refname), 'wb',
                                template=samfile)
        for s in samfile:
            bamfile.write(s)
        if VERBOSE:
            print 'DONE.'
    
    # Open BAM
    bamfile = pysam.Samfile(get_bam_file(data_folder, adaID, refname), 'rb')
    
    # Chromosome list
    chromosomes = bamfile.references
    
    # Read reference (fragmented)
    refseqs_raw = list(SeqIO.parse(get_ref_file(data_folder, adaID, refname),
                                   'fasta'))
    # Sort according to the chromosomal ordering
    refseqs = []
    for chromosome in chromosomes:
        for seq in refseqs_raw:
            if chromosome == seq.id:
                refseqs.append(seq)
                break
    refs = [np.array(refseq) for refseq in refseqs]
    
    # Allele counts and inserts
    counts_all = [np.zeros((len(read_types), len(alpha), len(s)), int)
                  for s in refseqs]
    inserts_all = [[defaultdict(int) for k in read_types]
                   for s in refseqs]

    # Count alleles
    # Iterate over pairs of reads
    for i_pairs, (read1, read2) in enumerate(pair_generator(bamfile)):

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
        for i_within_pair, read in enumerate([read1, read2]):
        
            # Find out on what chromosome the read has been mapped
            fragment = read.tid
            counts = counts_all[fragment]
            inserts = inserts_all[fragment]
            refseq = refseqs[fragment]
            ref = refs[fragment]
        
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

    
    # Build consensus for each fragment
    consensus_all = {}
    for chromosome, counts, inserts, refseq in \
                            izip(chromosomes, counts_all, inserts_all, refseqs):
    
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
                    print chromosome, pos_pre, pos_post, k, len(refseq), counts[js, :, pos_pre:pos_post].sum(axis=0)
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
    
        # Store
        consensus_all[chromosome] = consensus_final

    return refseqs, consensus_all


def check_new_old_consensi(refseqs, consensus_all):
    '''Check the old and the new consensi, whether they match'''
    
    # check whether the old consensus is exactly equal to the newly built one
    all_match = True
    for refseq in refseqs:
        if str(refseq.seq) != consensus_all[refseq.id]:
            all_match = False
            break

    return all_match


def write_consensus(data_folder, refname, adaID, conensus_all):
    '''Write consensus sequences into FASTA'''
    consensusseqs = [SeqRecord(Seq(consensus),
                               id=fname+' '+'{:02d}'.format(adaID)+'_consensus',
                               name=fname+' '+'{:02d}'.format(adaID)+'_consensus')
                     for fname, consensus in consensus_all.iteritems()]
    SeqIO.write(consensusseqs,
                get_ref_file(data_folder, adaID, next_refname(refname)),
                'fasta')



# Script
if __name__ == '__main__':

    # Parse input args: this is used to call itself recursively
    parser = argparse.ArgumentParser(description='Map HIV reads recursively')
    parser.add_argument('--stage', metavar='N', type=int, nargs='?',
                        default=1,
                        help=('1: initialize, 2: map, 3: build consensus'))
    parser.add_argument('--adaID', metavar='00', type=int, nargs='?',
                        default=0,
                        help='Adapter ID sample to map')
    parser.add_argument('--reference', metavar='REF', default='HXB2',
                        help=('What reference to use for mapping/consensus '+
                              'building:\n'+
                              '- HXB2: use that reference\n'+
                              '- c0: use the raw consensus\n'+
                              '- c1: use the 1st refined consensus\n'+
                              '- etc.'))
    parser.add_argument('--iterations', type=int, default=0,
                        help=('Maximal number of map/consensus iterations'))
    parser.add_argument('--verbose', type=int, default=0,
                        help=('Verbosity level [0-3]'))

    args = parser.parse_args()
    stage = args.stage
    adaID = args.adaID
    refname = transform_reference(args.reference)
    iterations_max = args.iterations
    VERBOSE = args.verbose

    # MAIN RECURSION LOOP (IT'S BROKEN AT THE END OF THE SCRIPT)
    # Note: exceptions also stop the loop.
    done = False
    while not done:
    
        ###########################################################################
        # 1. PREPARE HXB2 OR CONSENSUS FOR MAPPING
        ###########################################################################
        if stage == 1:

            if VERBOSE >= 2:
                print str(get_iteration(refname))+': entered stage 1...'

            # Make index and hash files out of reference/consensus sequence
            make_index_and_hash(data_folder, refname, adaID=adaID)

            # This was easy, move to the mapping part
            stage = 2

            if VERBOSE >= 2:
                print 'Exit stage 1.'
    
        ###########################################################################
        # 2. MAP AGAINST HXB2 OR CONSENSUS
        ###########################################################################
        if stage == 2:

            if VERBOSE >= 2:
                print str(get_iteration(refname))+': entered stage 2...'
    
            # If the script is called with no adaID, make it a master script for all
            # adapters: otherwise, keep it a master script for its own stampy only
            if adaID == 0:
                adaIDs = load_adapter_table(data_folder)['ID']
            else:
                adaIDs = [adaID]
    
            # Call stampy for mapping and get the jobIDs
            jobIDs = [call_stampy_for_mapping(data_folder, refname, adaID) for adaID in adaIDs]
    
            # Wait for all stampy children to finish, then either proceed to
            # consensus building (if only one child is left) or fork an instance
            # of self to that stage 3
            mapping_is_done = np.zeros(len(jobIDs), bool)
            while not mapping_is_done.all():
                time.sleep(interval_check)
    
                if VERBOSE >= 3:
                    print 'Mapping progress: ', list(map(int, mapping_is_done))
    
                # Ask qstat about our stampy jobs
                qstat_output = sp.check_output(['qstat'])
                if VERBOSE >= 3:
                    print qstat_output

                # Check all unfinished mappings
                for i, adaID in enumerate(adaIDs):
                    if mapping_is_done[i]:
                        continue

                    # If that mapping is done, proceed to consensus building
                    # FIXME: this check is rather lousy
                    is_in_qstat_output = jobIDs[i] in qstat_output
                    if not is_in_qstat_output:
                        mapping_is_done[i] = True

                        # Fork self to a child unless there is only one process
                        if len(adaIDs) > 1:
                            fork_self(data_folder, refname, adaID, stage=3,
                                      iterations_max=iterations_max)
                            done = True
                        else:
                            stage = 3

            if VERBOSE >= 2:
                print 'Exit stage 2.'

    
        ###########################################################################
        # 3. MAKE CONSENSUS AFTER MAPPING
        ###########################################################################
        if stage == 3:

            if VERBOSE >= 2:
                print str(adaID)+', '+str(get_iteration(refname))+': entered stage 3...'

            # Make consensus
            refseqs, consensus_all = make_consensus(data_folder, refname, adaID)

            # Check whether the old and new consensi match
            all_match = check_new_old_consensi(refseqs, consensus_all)
            if all_match:
                print 'Consensus converged.'
                done = True

            else:
                # Save consensi to file
                write_consensus(data_folder, refname, adaID, consensus_all)

                # Check whether we are at the max number of iterations
                if (get_iteration(refname) + 1 >= iterations_max):
                    print 'Maximal number of iterations reached.'
                    done = True
                else:
                    if VERBOSE:
                        print 'Starting again for iteration '+str(get_iteration(refname)+1)
                    refname = next_refname(refname)
                    stage = 1

            if VERBOSE >= 2:
                print 'Exit stage 3.'
