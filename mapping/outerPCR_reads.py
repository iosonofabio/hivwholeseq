# vim: fdm=marker
'''
author:     Fabio Zanini
date:       16/09/13
content:    Figure out what those errors at 10-15 and 235-240 from the fragment
            ends are.
'''
# TODO: convert this script into a general script to check the abundance of outer
# PCR-amplified reads
# Modules
import os
import argparse
import pysam
import numpy as np
from Bio import SeqIO

# Horizontal import of modules from this folder
from mapping.adapter_info import load_adapter_table
from mapping.miseq import read_types
from mapping.filenames import get_mapped_filename, get_allele_counts_filename, \
        get_coverage_filename, get_consensus_filename
from mapping.mapping_utils import sort_bam, index_bam
from mapping.minor_allele_frequency import get_minor_allele_counts



# Globals
# FIXME
from mapping.datasets import dataset_testmiseq as dataset
data_folder = dataset['folder']



# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Get allele counts')
    parser.add_argument('--adaIDs', nargs='*', type=int,
                        help='Adapter IDs to analyze (e.g. 2 16)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--errind', type=int, default=None, nargs='+',
                        help='Index of the errors to focus on')

    args = parser.parse_args()
    adaIDs = args.adaIDs
    VERBOSE = args.verbose
    errind = args.errind

    # If the script is called with no adaID, iterate over all
    if not adaIDs:
        adaIDs = load_adapter_table(data_folder)['ID']
    if VERBOSE >= 3:
        print 'adaIDs', adaIDs

    # Restrict to fragment 4, we know the problem is there
    fragment = 'F4'

    # Iterate over all requested samples
    table = load_adapter_table(data_folder)
    for adaID in adaIDs:
        samplename = table['sample'][table['ID'] == adaID][0]

        counts = np.load(get_allele_counts_filename(data_folder, adaID, fragment))
        coverage = np.load(get_coverage_filename(data_folder, adaID, fragment))

        (counts_major,
         counts_minor,
         counts_minor2) = get_minor_allele_counts(counts, n_minor=2)

        # Get minor allele frequencies and identities
        nu_minor = 1.0 * counts_minor[:, :, 1] / (coverage + 1e-6)
        all_minor = counts_minor[:, :, 0]

        # Find the problems
        pos_errs = (nu_minor > 0.3e-2).all(axis=0).nonzero()[0]
        all_errs = all_minor[:, pos_errs]

        # Read reference
        reffilename = get_consensus_filename(data_folder, adaID, fragment)
        refseq = SeqIO.read(reffilename, 'fasta')
        ref = np.array(refseq)
 
        # Find the reads overlapping with that position in the reference
        bamfilename = get_mapped_filename(data_folder, adaID, fragment, type='bam',
                                          filtered=True, sort=True)

        # Sort BAM if needed
        if not os.path.isfile(bamfilename):
            sort_bam(bamfilename)
            index_bam(bamfilename)

        # The indices of errors to study (manual)
        if errind is None:
            errind = np.arange(len(pos_errs))
        for indi in errind:
            pos_err = pos_errs[indi]
            all_err = all_errs.T[indi]
            #with pysam.Samfile(bamfilename, 'rb') as bamfile:
            bamfile = pysam.Samfile(bamfilename, 'rb')
            if True:
                print bamfile.references
                reference = '_'.join(['adaID', '{:02d}'.format(adaID),
                                  fragment, 'consensus'])

                # Fetch the reads that cover the focal site
                read_it = bamfile.fetch(reference, pos_err, pos_err+1)

                # Divide reads by allele at the error position
                if VERBOSE >= 1:
                    print 'Classifying read into types...'
                reads = {k: [[] for _ in read_types] for k in ['WT', 'mut']}
                for read in read_it:

                    # Divide by read 1/2 and forward/reverse
                    js = 2 * read.is_read2 + read.is_reverse

                    # Get over the CIGARs until that position
                    seq = read.seq
                    pos = read.pos

                    for (block_type, block_len) in read.cigar:
                        # Deletions are boring, skip
                        if block_type == 2:
                            pos += block_len

                        # Inserts slightly less boring, but still
                        elif block_type == 1:
                            seq = seq[block_len:]

                        # Matches are cool
                        elif block_type == 0:
                            if pos + block_len > pos_err:
                                # Look at the allele at that position
                                al = seq[pos_err - pos]
                                if al != ref[pos_err]:#FIXME == alpha[all_err[js]]:
                                    reads['mut'][js].append(read)
                                else:
                                    reads['WT'][js].append(read)
                                break

                            pos += block_len
                            seq = seq[block_len:]

                        # The rest is suspicious
                        else:
                            raise ValueError('CIGAR type '+str(block_type)+' not recognized')

                ## Get the mate reads (probably getting names and reiterating is faster...)
                #if VERBOSE >= 1:
                #    print 'Finding mate reads...'
                #reads_mate = {k: [[] for _ in read_types] for k in ['WT', 'mut']}
                #for key in ['mut']:#reads:
                #    readst = reads[key]
                #    reads_matet = reads_mate[key]
                #    for js, read_type in enumerate(read_types):
                #        for i, read in enumerate(readst[js]):
                #            if VERBOSE >= 3:
                #                if not ((i+1) % 10):
                #                    print i+1
                #            reads_matet[js].append(bamfile.mate(read))

                # Get the insert sizes
                from operator import attrgetter
                isizes_err = [map(attrgetter('isize'), reads['mut'][js])
                              for js, read_type in enumerate(read_types)]

                # Plot the histograms
                import matplotlib.pyplot as plt
                import matplotlib.cm as cm
                fig, ax = plt.subplots(1, 1)
                colors = [cm.jet(int(255.0 * js / len(read_types)))
                          for js, read_type in enumerate(read_types)]
                ax.hist(map(np.abs, isizes_err), bins=np.linspace(0, 600, 30),
                        normed=True,
                        color=colors, alpha=1, label=read_types)

                # Plot the distance of the primer from the edge of the fragment
                from sequence_utils.annotate_HXB2 import load_HXB2
                HXB2 = load_HXB2()
                for fea in HXB2.features:
                    if 'inner PCR' not in fea.type:
                        continue
                    if ((pos > 1000) and \
                        (fea.type[-1] == 'S') and \
                        ('F'+str(int(fragment[-1]) + 1) in fea.type)) or\
                       ((pos < 1000) and \
                        (fea.type[-1] == 'R') and \
                        ('F'+str(int(fragment[-1]) - 1) in fea.type)):
                            break
                primer = str(fea.extract(HXB2).seq)
                streches = primer.split(ref[pos_err])
                strech = streches[np.argmax(map(len, streches))]
                primer_pos = str(refseq.seq)[pos_err - 20: pos_err + 20].find(strech) - primer.find(strech) + pos_err - 20
                if pos_err < 1000:
                    line_x = [primer_pos, primer_pos + len(primer)]
                else:
                    line_x = [len(ref) - (primer_pos + len(primer)), len(ref) - primer_pos]
                plt.axvline(line_x[0], lw=2, c='r', ls='--')
                plt.axvline(line_x[1], lw=2, c='r', ls='--')

                ax.set_xlabel('Insert size of wrong reads')
                ax.set_ylabel('Distribution')
                ax.legend(loc=1)

                ax.set_title(samplename+' '+fragment+', error position: '+str(pos_err),
                             fontsize=20)
                plt.tight_layout(rect=(0, 0, 1, 0.97))
                plt.ion()
                plt.show()

                # Pick the strange ones
                rm = reads['mut']
                rmfoc = rm[1]
                rmg = [r for r in rmfoc if np.abs(r.isize) > len(ref) - pos_err + 10]
                ind = np.arange(len(rmg))
                np.random.shuffle(ind)
                rms = np.array(rmg)[ind[:10]]
                pai = [(r, bamfile.mate(r)) for r in rms]
                print [[[p.pos, p.pos + sum(bl for bt, bl in p.cigar if bt in [0, 2])] for p in pa] for pa in pai]




                ## Analyze the reads
                #start_poss = []
                #from collections import Counter
                #for key in reads.iterkeys():
                #    print key
                #    # Extract the start positions of the reads (only frequenct ones)
                #    for i, read_type in enumerate(read_types):
                #        start_pos_hist = Counter([r.pos for r in reads[key][i]])
                #        if key == 'mut':
                #            start_poss.append(start_pos_hist.most_common(1)[0][0])
                #        print start_pos_hist.most_common(5)
                #    print

                #    for i, read_type in enumerate(read_types):
                #        cigar_hist = Counter([tuple(r.cigar) for r in reads[key][i]])
                #        print cigar_hist.most_common(5)
                #    print
                #    print    

                ## Align the error reads via Waterman-Smith
                #import re
                #from operator import attrgetter
                #from Bio.Seq import Seq
                #from Bio.SeqRecord import SeqRecord
                #from Bio.Alphabet import IUPAC
                #import Bio.SeqIO as SeqIO
                #from Bio.Emboss.Applications import WaterCommandline

                #cwd = os.path.abspath(os.curdir)+'/'
                #if not os.path.isdir(cwd+'tmp'):
                #    os.mkdir(cwd+'tmp')

                #seqs_err = []
                #for js, read_type in enumerate(read_types):
                #    seq_err = reads_to_seqrecord(reads['mut'][js])
                #    seqs_err.append(seq_err)

                #    #rt = re.sub(r' ', r'_', read_type)
                #    #refrfilename = 'tmp/adaID_'+'{:02d}'.format(adaID)+'_'+fragment+'_'+str(pos_err)+'_'+rt+'_ref.fa'
                #    #SeqIO.write(refseq[start_poss[js]:min(len(refseq), start_poss[js] + 250)],
                #    #            cwd+refrfilename, 'fasta')

                #    #readsfilename = 'tmp/adaID_'+'{:02d}'.format(adaID)+'_'+fragment+'_'+str(pos_err)+'_'+rt+'.fa'
                #    #alifilename = 'tmp/adaID_'+'{:02d}'.format(adaID)+'_'+fragment+'_'+str(pos_err)+'_'+rt+'_aligned.txt'
                #    #SeqIO.write(seq_err, readsfilename, 'fasta')
                #    #cmd = WaterCommandline(asequence=cwd+refrfilename,
                #    #                       bsequence=cwd+readsfilename,
                #    #                       outfile=cwd+alifilename,
                #    #                       gapopen=10, gapextend=3)
                #    #print cmd
                #    #cmd()
                #    
                ## Check phred scores along the read
                #seqs_right = []
                #for js, read_type in enumerate(read_types):
                #    seq_right = reads_to_seqrecord(reads['WT'][js])
                #    seqs_right.append(seq_right)

                #scores_err = [np.zeros(250) for _ in read_types]
                #for js, seq_err in enumerate(seqs_err):
                #    cov = np.zeros(250)
                #    for seq in seq_err:
                #        scores_err[js][:len(seq)] += seq.letter_annotations['phred_quality']
                #        cov[:len(seq)] += 1
                #    scores_err[js] = 1.0 * scores_err[js] / cov
                #scores_right = [np.zeros(250) for _ in read_types]
                #for js, seq_right in enumerate(seqs_right):
                #    cov = np.zeros(250)
                #    for seq in seq_right:
                #        scores_right[js][:len(seq)] += seq.letter_annotations['phred_quality']
                #        cov[:len(seq)] += 1
                #    scores_right[js] = 1.0 * scores_right[js] / cov

                #import matplotlib.pyplot as plt
                #fig, axs = plt.subplots(2, 2, figsize=(14, 12))
                #axs = axs.ravel()
                #for js, read_type in enumerate(read_types):
                #    ax = axs[js]
                #    ax.plot(np.arange(250), scores_err[js], label='mut', lw=2)
                #    ax.plot(np.arange(250), scores_right[js], label='WT', lw=2)
                #    if js > 1:
                #        ax.set_xlabel('Position in read')
                #    if js in [0, 2]:
                #        ax.set_ylabel('mean phred')
                #    ax.set_title(read_type)
                #    ax.grid()
                #    ax.axvline(pos_err - start_poss[js], lw=2, c='r', ls='--')

                #ax.legend(loc=3)

                #fig.suptitle(samplename+' '+fragment+', error position: '+str(pos_err),
                #             fontsize=20)
                #plt.tight_layout(rect=(0, 0, 1, 0.97))
                #plt.ion()
                #plt.show()

                ## Check positions of errors in reads
                ## Note: this must be done with the mapped reads
                #errs_read = {k: [np.zeros(250) for _ in read_types] for k in ['mut', 'WT']}
                #covs_read = {k: [np.zeros(250) for _ in read_types] for k in ['mut', 'WT']}
                #for (key, errs_readm) in errs_read.iteritems():
                #    covs_readm = covs_read[key]
                #    for (err_read, cov_read, readst) in izip(errs_readm, covs_readm, reads[key]):
                #        for read in readst:
                #            pos = read.pos
                #            seq = read.seq
                #            pos_read = 0
        
                #            for (block_type, block_len) in read.cigar:
        
                #                # Deletions are boring, skip
                #                if block_type == 2:
                #                    pos += block_len
        
                #                # Inserts slightly less boring, but still
                #                elif block_type == 1:
                #                    seq = seq[block_len:]
                #                    pos_read += block_len
        
                #                # Matches are cool
                #                elif block_type == 0:
                #                    # Find mismatches and score them
                #                    seqb = np.array(list(seq[:block_len]))
                #                    seqr = ref[pos: pos + block_len]
                #                    err_read[(seqb != seqr).nonzero()[0] + pos_read] += 1
                #                    cov_read[pos_read: pos_read + block_len] += 1
        
                #                    pos += block_len
                #                    seq = seq[block_len:]
                #                    pos_read += block_len
        
                #                # The rest is suspicious
                #                else:
                #                    raise ValueError('CIGAR type '+str(block_type)+' not recognized')

                ## Plot it
                #import matplotlib.pyplot as plt
                #colors = {'mut': 'b', 'WT': 'g'}
                #fig, axs = plt.subplots(2, 2, figsize=(14, 12))
                #axs = axs.ravel()
                #for js, read_type in enumerate(read_types):
                #    ax = axs[js]
                #    for key in errs_read:
                #        # Normalize by the mut/WT coverage
                #        ax.plot(np.arange(250), 1.0 * errs_read[key][js] / covs_read[key][js], label=key, lw=2, c=colors[key])
                #        ax.scatter(np.arange(250), 1.0 * errs_read[key][js] / covs_read[key][js], lw=2, s=50,
                #                   color=colors[key], edgecolor='none')
                #    if js > 1:
                #        ax.set_xlabel('Position in read')
                #    if js in [0, 2]:
                #        ax.set_ylabel('number of errors')
                #    ax.set_title(read_type)
                #    ax.set_yscale('log')
                #    ax.grid()
                #    ax.axvline(pos_err - start_poss[js], lw=2, c='r', ls='--')

                #ax.legend(loc=3)

                #fig.suptitle(samplename+' '+fragment+', error position: '+str(pos_err),
                #             fontsize=20)
                #plt.tight_layout(rect=(0, 0, 1, 0.97))
                #plt.ion()
                #plt.show()

                ### Plot histograms (because WT is samples much better)
                ##fig, axs = plt.subplots(2, 2, figsize=(14, 12))
                ##axs = axs.ravel()
                ##for js, read_type in enumerate(read_types):
                ##    ax = axs[js]
                ##    for key in errs_read:
                ##        # Fit best binomial
                ##        from scipy.optimize import curve_fit
                ##        from scipy.stats import binom
                ##        f = lambda x, p: binom.pmf(x, len(errs_read[key][js]), p)
                ##        p = curve_fit(f, np.arange(19), [1e-3])[0][0]

                ##        ax.hist(errs_read[key][js],
                ##                bins=np.arange(20), normed=True, alpha=0.5,
                ##                label=key+', P = '+'{:4f}'.format(p))

                ##        ax.legend(loc=1)
                ##    if js > 1:
                ##        ax.set_xlabel('Position in read')
                ##    if js in [0, 2]:
                ##        ax.set_ylabel('number of errors')
                ##    ax.set_title(read_type)
                ##    ax.grid()


                ##fig.suptitle(samplename+' '+fragment+', error position: '+str(pos_err),
                ##             fontsize=20)
                ##plt.tight_layout(rect=(0, 0, 1, 0.97))
                ##plt.ion()
                ##plt.show()

                #
                ## Check motifs

    
                #    # TODO: 1. check phred scores...OK
                #    #       2. check error positions in the read...
                #    #       3. check DNA motifs...
