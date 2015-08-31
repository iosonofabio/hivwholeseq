# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/11/13
content:    Find out what's wrong with fragment F6.
'''
# Modules
import os
import numpy as np
import pysam
from Bio import SeqIO
from Bio.Seq import reverse_complement

from mapping.miseq import alpha
from mapping.datasets import MiSeq_runs
from mapping.filenames import get_mapped_filename, get_consensus_filename, \
        get_divided_filenames
from mapping.mapping_utils import sort_bam, index_bam, pair_generator, convert_sam_to_bam
from mapping.one_site_statistics import get_allele_counts_insertions_from_file_unfiltered, \
        filter_nus



# Script
if __name__ == '__main__':

    # The dip in coverage happens around position 977
    pos_dip = 979

    # General parameters
    miseq_run = 28
    adaID = 4
    fragment = 'F6'
    dataset = MiSeq_runs[miseq_run]
    data_folder = dataset['folder']

    ## Input read file
    #bamfilename = get_mapped_filename(data_folder, adaID, fragment, type='bam',
    #                                  filtered=True, sort=True)

    ## Sort and index file
    #if not os.path.isfile(bamfilename):
    #    sort_bam(bamfilename)
    #    index_bam(bamfilename)

    ## Get all the reads covering that position first
    #bamfile = pysam.Samfile(bamfilename, 'rb')
    #read_it = bamfile.fetch(bamfile.references[0], pos_dip, pos_dip+1)
    #for read in read_it:
    #    if read.is_reverse:
    #        print 'rev',
    #    else:
    #        print 'fwd',

    #    # Get the allele at position 979
    #    pos_all = 979
    #    pos_ref = read.pos
    #    pos_read = 0
    #    for (bt, bl) in read.cigar:
    #        if bt == 1:
    #            pos_read += bl
    #        elif bt == 2:
    #            if pos_ref + bl > pos_all:
    #                al = '-'
    #                break
    #            pos_ref += bl
    #        elif bt == 0:
    #            if pos_ref + bl > pos_all:
    #                al = read.seq[pos_read + pos_all - pos_ref]
    #                if (pos_ref < pos_all - 5) and (pos_ref + bl > pos_all + 3):
    #                    motif = read.seq[pos_read + pos_all - pos_ref - 4: pos_read + pos_all - pos_ref + 4]
    #                else:
    #                    motif = ''
    #                break
    #            pos_read += bl
    #            pos_ref += bl

    #    print read.pos, read.pos + sum(bl for (bt, bl) in read.cigar if bt in (0, 2)), al, motif, read.isize

    # Look for HIV signatures in the trashed reads
    consensus_filename = get_consensus_filename(data_folder, adaID, fragment)
    consseq = SeqIO.read(consensus_filename, 'fasta')
    conss = str(consseq.seq)

    # Get the sequence of the suspicious plasmid
    ### Lentiviral vector
    ##vector_name = 'KC262216.1'
    ## Plant vector
    #vector_name = 'KF499077.1'
    # Plant vector 2
    vector_name = 'KF206147.1'
    vector_filename = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/reference/'+\
            'vectors/'+vector_name+'.fasta'
    vecseq = SeqIO.read(vector_filename, 'fasta')
    vecmat = np.array(vecseq)
    vecs = str(vecseq.seq)
    vecsrev = reverse_complement(str(vecseq.seq))

    ## Take a few random signatures to the left and right of the site
    #signatures_left = [conss[600: 615],
    #                   conss[700: 715],
    #                   conss[800: 815],
    #                   conss[900: 915],
    #                   conss[920: 935],
    #                   conss[940: 955]]
    #signatures_right = [conss[1000: 1015],
    #                    conss[1020: 1035],
    #                    conss[1040: 1055],
    #                    conss[1100: 1115],
    #                    conss[1200: 1215]]
    #signatures_vector = [vecs[i: i+15] for i in xrange(0, len(vecs), 20)]
    #signatures_vectorrev = [vecsrev[i: i+15] for i in xrange(0, len(vecsrev), 20)]
    #trashed_filename = get_divided_filenames(data_folder, adaID, 'F6')[-2]
    #vector_reads_filename = os.path.dirname(trashed_filename)+'/vector_reads.bam'
    #with pysam.Samfile(trashed_filename, 'rb') as trashed_file:
    #    with pysam.Samfile(vector_reads_filename, 'wb', template=trashed_file) as vector_file:

    #        for reads in pair_generator(trashed_file):
    #            if reads[0].is_unmapped:
    #                possvs = []
    #                possvsrev = []
    #                foundvs = []
    #                foundvsrev = []
    #                stop = False
    #                for read in reads:
    #                    # Look only for right or only for left
    #                    poss = map(read.seq.find, signatures_left)
    #                    found = np.array([p != -1 for p in poss], int)
    #                    if found.any():
    #                        print ' '.join(map(str, found)),
    #                        print '*'
    #     
    #                    # Signatures of HIV found, now look for signatures of vector
    #                    possv = map(read.seq.find, signatures_vector)
    #                    foundv = np.array([p != -1 for p in possv], int)
    #                    possvs.append(possv)
    #                    foundvs.append(foundv)
    #                    possvrev = map(read.seq.find, signatures_vectorrev)
    #                    foundvrev = np.array([p != -1 for p in possvrev], int)
    #                    possvsrev.append(possvrev)
    #                    foundvsrev.append(foundvrev)

    #                    if (foundv.sum() >= 3) or (foundvrev.sum() >= 3):
    #                        stop = True
    #
    #                if stop:
    #                    vector_file.write(reads[0])
    #                    vector_file.write(reads[1])
    #                    #import ipdb; ipdb.set_trace()

    ## Get the mapped reads and make a consensus
    #mapped_reads_filename = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/'+\
    #        'run28_test_samples/adapterID_04/contamination_vector/mapped_reads_to_'+vector_name+'.bam'
    #if not os.path.isfile(mapped_reads_filename):
    #    convert_sam_to_bam(mapped_reads_filename)
    #(counts, inserts) = get_allele_counts_insertions_from_file_unfiltered(mapped_reads_filename, len(vecs),
    #                                                                      qual_min=25, maxreads=10000)
    #nus = filter_nus(counts, counts.sum(axis=1))
    #consmat = alpha[nus.argmax(axis=0)]
    #conss = ''.join(consmat)

    #covered = counts.sum(axis=1).sum(axis=0) > 0

    # Get a bunch of reads from fragment F6, save them to interleaved fastq
    # (for spades to assemble)
    maxreads = 1000
    from mapping.mapping_utils import reads_to_seqrecord
    ftmp = get_divided_filenames(data_folder, adaID, ['F6'])
    (F6_filename, trashed_filename) = (ftmp[0], ftmp[2])
    with pysam.Samfile(F6_filename, 'rb') as F6f:
        reads_F6 = []
        for irp, reads in enumerate(pair_generator(F6f)):
            if irp > maxreads:
                break
            reads_F6.append(reads[0])
            reads_F6.append(reads[1])
        seqs_F6 = reads_to_seqrecord(reads_F6)

    # Save some seqs from the unmapped
    maxreads = 10000
    with pysam.Samfile(trashed_filename, 'rb') as trash_f:
        reads_trash = []
        for irp, reads in enumerate(pair_generator(trash_f)):
            if irp > maxreads:
                break
            reads_trash.append(reads[0])
            reads_trash.append(reads[1])
        seqs_trash = reads_to_seqrecord(reads_trash)

    # Write all reads out for assembly
    assembly_file = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/'+\
            'run28_test_samples/adapterID_04/contamination_vector/F6_trash_ass.fastq'
    SeqIO.write(seqs_F6+seqs_trash, assembly_file, 'fastq')

    ## Design fwd primer for region ~700-800 to Sanger sequence the stuff
    #cmat_other = [np.array(SeqIO.read(get_consensus_filename(data_folder, adaID, F), 'fasta'))
    #              for F in ['F1', 'F2', 'F3', 'F4', 'F5', 'F6']] + [vecmat]
    #cmat = np.array(consseq) 
    #primers = []
    #prlen = 30
    #for pos in xrange(800, 700, -10):
    #    pr = cmat[pos: pos + prlen]
    #    # Look for approximate matches in all the consensi
    #    from collections import Counter
    #    n_matches = Counter()
    #    for cmato in cmat_other:
    #        for poso in xrange(len(cmato) - prlen):
    #            nm = (cmato[poso: poso + prlen] == pr).sum()
    #            if nm > prlen / 2:
    #                n_matches[nm] += 1
    #    from operator import itemgetter
    #    print ''.join(pr), pos, sorted(n_matches.items(), key=itemgetter(0), reverse=True)
