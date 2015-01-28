# vim: fdm=marker
'''
author:     Fabio Zanini
date:       02/04/14
content:    Get reads from the V3 loop from any one sample for Thomas Leitner,
            with some data processing.
'''
# Modules
import argparse
import sys
import os
from itertools import izip
import numpy as np
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet.IUPAC import ambiguous_dna

from seqanpy import align_global

from hivwholeseq.patients.patients import get_patient
from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.patients.filenames import get_mapped_to_initial_filename, \
        get_initial_consensus_filename
from hivwholeseq.reference import load_HXB2
from hivwholeseq.utils.mapping import pair_generator




# Script
if __name__ == '__main__':

    # Input arguments
    parser = argparse.ArgumentParser(description='Extract seqs for Thomas')
    parser.add_argument('--PCRtype', default=['PCR1', 'PCR2'], nargs='+',
                        help='PCR1 or PCR2')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    PCRtypes = args.PCRtype
    VERBOSE = args.verbose

    VERBOSE = 1
    pname = '20097'
    time = 336
    fragment = 'F5'
    maxgap = 0

    patient = get_patient(pname)

    cos = {}
    for PCRtype in PCRtypes:
        if VERBOSE:
            print PCRtype
        samplename = [s for (s, t) in izip(patient.samples, patient.times()) if \
                      (t == time) and PCRtype in s][0]
    
    
        # Make sure it's a multiple of three, in case they translate
        # Plus the coordinates of Jan are a bit fuzzy, because the inner fwd primer
        # length is only 25, not 30...
        pos_V3_HXB2 = [6983, 7352]
    
        ref_rec = load_HXB2()
        refm = np.array(ref_rec)
        V3ref = ref_rec[pos_V3_HXB2[0]: pos_V3_HXB2[1]].seq
    
        cons_rec = SeqIO.read(get_initial_consensus_filename(pname, fragment), 'fasta')
        cons = cons_rec.seq
        consm = np.array(cons_rec)
    
        ## Find the coordinates in the consensus
        V3primer_inner_fwd = np.fromstring('ACAATGYACACATGGAATTARGCCA', 'S1')
        seed = np.ma.array(V3primer_inner_fwd)
        seed[(seed == 'Y') | (seed == 'R')] = np.ma.masked
        sl = len(seed)
        n_matches = np.array([(consm[i: i + sl] == seed).sum() for i in xrange(len(consm) - sl)])
        start = (n_matches).argmax()
        if n_matches[start] < (sl - seed.mask.sum()) * 0.75:
            raise ValueError('Seed not found reliably')
        start += sl
        print 'V3 primer: fwd', start
    
        from Bio.Seq import reverse_complement as rc
        V3primer_inner_rev = np.fromstring(rc('AGAAAAATTCYCCTCYACAATTAAA'), 'S1')
        seed = np.ma.array(V3primer_inner_rev)
        seed[(seed == 'Y') | (seed == 'R')] = np.ma.masked
        sl = len(seed)
        n_matches = np.array([(consm[i: i + sl] == seed).sum() for i in xrange(len(consm) - sl)])
        end = (n_matches).argmax()
        if n_matches[end] < (sl - seed.mask.sum()) * 0.75:
            raise ValueError('Seed not found reliably')
        print 'V3 primer: rev', end
    
        V3con = cons[start: end]
        V3s = start
        V3e = end
        V3l = V3e - V3s
        print V3con.translate()
    
        # Go down to the reads (unfiltered, do a liberal filtering here)
        bamfilename = get_mapped_to_initial_filename(pname, samplename, fragment, filtered=False)
        with pysam.Samfile(bamfilename, 'rb') as bamfile:
            reads_V3 = []
            
            from hivwholeseq.sequencing.filter_mapped_reads import check_overhanging_reads, \
                    get_distance_from_consensus, trim_bad_cigar

            # Iterate over all pairs
            n_good = 0
            n_wrongname = 0
            n_unmapped = 0
            n_unpaired = 0
            n_mutator = 0
            n_suspect = 0
            n_mismapped_edge = 0
            n_badcigar = 0
            for read_pair in pair_generator(bamfile):

                #####################################
                # SOFT FILTER 
                #####################################
                # Assign names
                (read1, read2) = read_pair
                i_fwd = read_pair[0].is_reverse
                i_rev = not i_fwd

                # Check a few things to make sure we are looking at paired reads
                if read1.qname != read2.qname:
                    n_wrongname += 1
                    raise ValueError('Read pair '+str(irp)+': reads have different names!')

                # Ignore unmapped reads
                if read1.is_unmapped or read2.is_unmapped:
                    if VERBOSE >= 2:
                        print 'Read pair '+read1.qname+': unmapped'
                    n_unmapped += 1
                    continue
            
                # Ignore not properly paired reads (this includes mates sitting on
                # different fragments)
                if (not read1.is_proper_pair) or (not read2.is_proper_pair):
                    if VERBOSE >= 2:
                        print 'Read pair '+read1.qname+': not properly paired'
                    n_unpaired += 1
                    continue

                # Mismappings are sometimes at fragment edges:
                # Check for overhangs beyond the edge
                skip = check_overhanging_reads(read_pair, len(consm))
                if skip:
                    n_mismapped_edge += 1
                    continue

                # Trim the bad CIGARs from the sides, if there are any good ones
                match_len_min = 30
                trim_bad_cigars = 3
                skip = trim_bad_cigar(read_pair, match_len_min=match_len_min,
                                       trim_left=trim_bad_cigars,
                                       trim_right=trim_bad_cigars,
                                       cons=str(cons))
                if skip:
                    n_badcigar += 1
                    continue

                # Mismappings are often characterized by many mutations:
                # check the number of mismatches and skip reads with too many
                max_mismatches = 100
                dc = get_distance_from_consensus(consm, read_pair,
                                                 threshold=20,
                                                 VERBOSE=VERBOSE)
                if (dc > max_mismatches).any():
                    if VERBOSE >= 2:
                        print 'Read pair '+read1.qname+': too many mismatches '+\
                                '('+str(dc[0])+' + '+str(dc[1])+')'
                    n_mutator += 1
                    #import ipdb; ipdb.set_trace()
                    continue

                #####################################

                # Cover both edges
                start_fwd = read_pair[i_fwd].pos
                if start_fwd > V3s:
                    continue
    
                end_fwd = start_fwd + sum(bl for (bt, bl) in read_pair[i_fwd].cigar if bt in (0, 2))
                if end_fwd <= V3s:
                    continue
    
                start_rev = read_pair[i_rev].pos
                if start_rev >= V3e:
                    continue
    
                end_rev = start_rev + sum(bl for (bt, bl) in read_pair[i_rev].cigar if bt in (0, 2))
                if end_rev < V3e:
                    continue
    
                len_gap = start_rev - end_fwd
                if len_gap > maxgap:
                    continue
    
                reads_V3.append(read_pair)
            
            print '# of reads covering the segment:', len(reads_V3)
        
            # Precompute conversion table
            SANGER_SCORE_OFFSET = ord("!")
            q_mapping = dict()
            for letter in range(0, 255):
                q_mapping[chr(letter)] = letter - SANGER_SCORE_OFFSET
        
            reads_V3_fastq = []
            for read_pair in reads_V3:
                i_fwd = read_pair[0].is_reverse
                i_rev = not i_fwd
                name = read_pair[i_fwd].qname
                id = name
                descr = name
    
                start_fwd = read_pair[i_fwd].pos
                end_fwd = start_fwd + sum(bl for (bt, bl) in read_pair[i_fwd].cigar if bt in (0, 2))
                start_rev = read_pair[i_rev].pos
                end_rev = start_rev + sum(bl for (bt, bl) in read_pair[i_rev].cigar if bt in (0, 2))
                len_gap = start_rev - end_fwd
        
                if read_pair[i_fwd].isize < 0:
                    raise ValueError('Wrong isize?!')
        
                # FWD read
                # The insert is more than 300 + tolerance long, so we take the whole fwd read
                pos_ref = read_pair[i_fwd].pos
                pos_ref_start_rev = read_pair[i_rev].pos
                pos_read_start = None
                pos_read_start_rev = None
                pos_read = 0
                for (bt, bl) in read_pair[i_fwd].cigar:
                    if bt == 1:
                        pos_read += bl
        
                    elif bt == 2:
                        if (pos_read_start is None) and (pos_ref + bl > V3s):
                            pos_read_start = pos_read
                                
                        if (pos_read_start_rev is None) and (pos_ref + bl >= pos_ref_start_rev):
                            pos_read_start_rev = pos_read
        
                        if (pos_read_start is not None) and (pos_read_start_rev is not None):
                            break
        
                        pos_ref += bl
        
                    elif bt == 0:
                        if (pos_read_start is None) and (pos_ref + bl > V3s):
                            pos_read_start = pos_read + V3s - pos_ref
                                
                        if (pos_read_start_rev is None) and (pos_ref + bl >= pos_ref_start_rev):
                            pos_read_start_rev = pos_read + pos_ref_start_rev - pos_ref
        
                        if (pos_read_start is not None) and (pos_read_start_rev is not None):
                            break
        
                        pos_ref += bl
                        pos_read += bl
        
                if pos_read_start is None:
                    raise ValueError('Start not found!?')
        
                seq_fwd = read_pair[i_fwd].seq[pos_read_start:]
                qual_fwd = read_pair[i_fwd].qual[pos_read_start:]
        
                # REV read
                # Analogously, because the insert is big enough, we do not care about the start
                pos_ref = read_pair[i_rev].pos
                pos_read_end = None
                pos_read = 0
                for (bt, bl) in read_pair[i_rev].cigar:
                    if bt == 1:
                        pos_read += bl
                    elif bt == 2:
                        if pos_ref + bl >= V3e:
                            pos_read_end = pos_read
                            break
                        pos_ref += bl
        
                    elif bt == 0:
                        if pos_ref + bl >= V3e:
                            pos_read_end = pos_read + V3e - pos_ref
                            break
                        pos_ref += bl
                        pos_read += bl
        
                if pos_read_end is None:
                    raise ValueError('End not found!?')
        
                seq_rev = read_pair[i_rev].seq[:pos_read_end]
                qual_rev = read_pair[i_rev].qual[:pos_read_end]
        
                # Merge the two reads
                if pos_read_start_rev is None:
                    gap_len = V3l - len(seq_fwd) - len(seq_rev)
                    #if gap_len != len_gap:
                    #    import ipdb; ipdb.set_trace()
                    seq = seq_fwd+('N' * gap_len)+seq_rev
                    qual = qual_fwd+(chr(SANGER_SCORE_OFFSET) * gap_len)+qual_rev
                    if len(seq) != len(qual):
                        raise ValueError('Seq NO OVERLAP has a different length from qual!?')
        
                else:
                    #if VERBOSE >= 3:
                    #    print seq_fwd[-pos_read_start + pos_read_start_rev - 10: -pos_read_start + pos_read_start_rev + 10];
                    #    print (' ' * 10)+seq_rev[:10]
        
                    len_overlap = len(seq_fwd) - (pos_read_start_rev - pos_read_start)
                    overlap_ali = align_global(seq_fwd[-pos_read_start + pos_read_start_rev:],
                                               seq_rev[: len_overlap], band=10)
        
                    if VERBOSE >= 3:
                        print overlap_ali[1][:50]; print overlap_ali[1][:50]
        
                    # Check whether the two reads agree
                    seq_over = []
                    qual_over = []
                    pos_qual_fwd = pos_read_start_rev - pos_read_start
                    pos_qual_rev = 0
                    for i in xrange(len_overlap):
                        nuc_fwd = overlap_ali[1][i]
                        nuc_rev = overlap_ali[2][i]
                        if (nuc_fwd == '-') and (nuc_rev != '-'):
                            pos_qual_rev += 1
                            continue
                        elif (nuc_fwd != '-') and (nuc_rev == '-'):
                            pos_qual_fwd += 1
                            continue
                        elif (nuc_fwd == '-') and (nuc_rev == '-'):
                            continue
        
                        if nuc_fwd != nuc_rev:
                            seq_over.append('N')
                            qual_over.append('!')
        
                        else:
                            seq_over.append(nuc_fwd)
                            qual_over.append(chr(ord(qual_fwd[pos_qual_fwd]) + \
                                                 ord(qual_rev[pos_qual_rev]) - \
                                                 SANGER_SCORE_OFFSET))
        
                        pos_qual_fwd += 1
                        pos_qual_rev += 1
        
                    seq_over = ''.join(seq_over)
                    qual_over = ''.join(qual_over)
                    if len(seq_over) != len(qual_over):
                        raise ValueError('Seq OVER has a different length from qual!?')
        
                    seq_fwd = seq_fwd[: pos_read_start_rev - pos_read_start]
                    qual_fwd = qual_fwd[: pos_read_start_rev - pos_read_start]
                    if len(seq_fwd) != len(qual_fwd):
                        raise ValueError('Seq FWD has a different length from qual!?')
        
                    seq_rev = seq_rev[len_overlap:]
                    qual_rev = qual_rev[len_overlap:]
                    if len(seq_rev) != len(qual_rev):
                        raise ValueError('Seq REV has a different length from qual!?')
        
        
                    # Join FWD, OVERLAP, REV
                    seq = seq_fwd+seq_over+seq_rev
                    qual = qual_fwd+qual_over+qual_rev
                    if len(seq) != len(qual):
                        raise ValueError('Seq MERGED has a different length from qual!?')
                
                # Convert to FASTQ
                read_V3_fastq = SeqRecord(Seq(seq, ambiguous_dna),
                                          id=id, name=name, description=descr)
            
                # Get the qualities second
                qualities = [q_mapping[letter] for letter in qual]
                if qualities and (min(qualities) < 0 or max(qualities) > 93):
                    raise ValueError("Invalid character in quality string")
        
                #For speed, will now use a dirty trick to speed up assigning the
                #qualities. We do this to bypass the length check imposed by the
                #per-letter-annotations restricted dict. This is equivalent to:
                #record.letter_annotations["phred_quality"] = qualities
                dict.__setitem__(read_V3_fastq._per_letter_annotations,
                                 "phred_quality", qualities)
        
                reads_V3_fastq.append(read_V3_fastq)
        
        SeqIO.write(reads_V3_fastq, '/ebio/ag-neher/home/fzanini/tmp/seqs_V3_Thomas_'+PCRtype+'_'+str(time)+'.fastq', 'fastq')
        SeqIO.write(reads_V3_fastq, '/ebio/ag-neher/home/fzanini/tmp/seqs_V3_Thomas_'+PCRtype+'_'+str(time)+'.fasta', 'fasta')
    
        from collections import Counter
        seqs = [str(r.seq) for r in reads_V3_fastq]
        co = Counter(seqs)
        cos[PCRtype] = co

    ## Plot cumulative spectrum
    #import matplotlib.pyplot as plt
    #plt.figure()
    #for PCRtype, co in cos.iteritems():
    #    x = (np.arange(len(co)) + 1.0)
    #    y = 1.0 * np.sort(co.values())[::-1]
    #    y /= y.sum()
    #    plt.plot(x, y, lw=2, label=PCRtype)

    #plt.xscale('log')
    #plt.yscale('log')
    #plt.xlabel('rank')
    #plt.ylabel('fraction of seqs')

    #plt.title('time = '+str(time)+' days')
    #plt.legend(loc=3)
    #plt.ion()
    #plt.show()

    ## Make rank correlation
    #from operator import itemgetter
    #seqs_common_PCR1 = map(itemgetter(0), cos['PCR1'].most_common())
    #seqs_common_PCR2 = map(itemgetter(0), cos['PCR2'].most_common())
    #seqs_common = list(set(cos['PCR1'].keys() + cos['PCR2'].keys()))
    #for seq in seqs_common:
    #    if seq not in seqs_common_PCR1:
    #        seqs_common_PCR1.append(seq)
    #    if seq not in seqs_common_PCR2:
    #        seqs_common_PCR2.append(seq)

    ## Make frequency correlation
    #freqs = []
    #for seq in seqs_common:
    #    if seq in cos['PCR1']:
    #        freq = [cos['PCR1'][seq]]
    #    else:
    #        freq = [0]
    #    if seq in cos['PCR2']:
    #        freq.append(cos['PCR2'][seq])
    #    else:
    #        freq.append(0)
    #    freqs.append(freq)

    #freqs = np.array(freqs, float)
    #freqs /= freqs.sum(axis=0)
    #freqs = freqs.T

    #plt.figure()
    #plt.scatter(freqs[0] + 1e-5, freqs[1] + 1e-5)
    #plt.xlabel('Frequency in PCR1')
    #plt.ylabel('Frequency in PCR2')
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.xlim(0.9e-5, 1)
    #plt.ylim(0.9e-5, 1)
    #plt.title('time = '+str(time)+' days')

    #plt.ion()
    #plt.show()

