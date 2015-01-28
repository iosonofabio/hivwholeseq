# vim: fdm=marker
'''
author:     Fabio Zanini
date:       16/11/14
content:    Collect trajectories of RNA structures, using local haplotypes.
'''
# Modules
import os
import argparse
from operator import itemgetter
from itertools import izip
import numpy as np
from Bio import SeqIO, AlignIO
from seqanpy import align_global

from hivwholeseq.patients.patients import load_patient
from hivwholeseq.utils.sequence import pretty_print_pairwise_ali



# Functions
def get_fragment_roi(seq, start, end):
    '''Get the best fragment for a genomic region'''
    for fea in seq.features:
        if fea.type == 'fragment':
            for locpart in fea.location.parts:
                if (locpart.nofuzzy_start <= start) and \
                   (locpart.nofuzzy_end >= end):
                    return fea


def parse_ct_file_multiple(filename):
    '''Parse CT file with RNA structures from RNAstructure'''
    structs = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            fields = line.split()
            if 'ENERGY' in line:
                structs.append({'energy': float(fields[-2]),
                                'filename': filename,
                                'index': len(structs),
                                'pairs': {}})

            else:
                pos1 = int(fields[0]) - 1
                pos2 = int(fields[-2]) - 1
                if pos2 != -1:
                    structs[-1]['pairs'][pos1] = pos2
                    structs[-1]['pairs'][pos2] = pos1

    return structs
    

def predict_RNA_structure(seq, label='seq', maxstructs=1, VERBOSE=0):
    '''Predict RNA secondary structures using RNAstructure'''
    import os
    import subprocess as sp
    from hivwholeseq.utils.generic import mkdirs
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Alphabet.IUPAC import ambiguous_dna

    rna_fold_bin = '/ebio/ag-neher/home/fzanini/programs/RNAstructure_cli/exe/Fold'

    # Make tmp input file
    tmp_file_in = '/ebio/ag-neher/home/fzanini/tmp/RNAfold/'+label+'.fasta'
    tmp_file_out = '/ebio/ag-neher/home/fzanini/tmp/RNAfold/'+label+'.ct'
    mkdirs(os.path.dirname(tmp_file_in))
    seqrec = SeqRecord(Seq(seq, ambiguous_dna),
                       id=label,
                       name=label,
                       description='')
    SeqIO.write(seqrec, tmp_file_in, 'fasta')

    # Call RNAStructure with all the crap (env vars, etc)
    rna_tables = '/ebio/ag-neher/home/fzanini/programs/RNAstructure_cli/data_tables/'
    env = os.environ.copy()
    env['DATAPATH'] = rna_tables
    call_list = [rna_fold_bin, '-m', str(maxstructs), tmp_file_in, tmp_file_out]
    if VERBOSE >= 2:
        print ' '.join(call_list)
    output = sp.check_output(call_list, shell=False)
    if VERBOSE >= 3:
        print output

    if 'Writing output ct file...done.' in output:
        structs = parse_ct_file_multiple(tmp_file_out)
    else:
        IOError('RNAstructure had problems predicting the structure')

    return structs


def plot_circle_compare(filename_ct1, filename_ct2, VERBOSE=0):
    '''Plot circle comparison of two structures'''
    import os
    import subprocess as sp

    cc_bin = '/ebio/ag-neher/home/fzanini/programs/RNAstructure_cli/exe/CircleCompare'
    filename_out = filename_ct1.replace('.ct', '.svg')
    call_list = [cc_bin, '--svg', '-n', '1', filename_ct1, filename_ct2, filename_out]
    if VERBOSE >= 2:
        print ' '.join(call_list)
    output = sp.check_output(call_list, shell=False)
    if VERBOSE >= 3:
        print output
    return filename_out


# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Get local haplotypes',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--structure', required=True,
                        help='RNA structure of interest')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-4]')
    parser.add_argument('--maxreads', type=int, default=-1,
                        help='Number of reads analyzed per sample')
    parser.add_argument('--plot', action='store_true',
                        help='Plot local haplotype trajectories')
    #parser.add_argument('--save', default=None,
    #                    help='Save to this filename')

    args = parser.parse_args()
    pname = args.patient
    stname = args.structure
    VERBOSE = args.verbose
    maxreads = args.maxreads
    use_plot = args.plot
    #save_path = args.save

    patient = load_patient(pname)
    patient.discard_nonsequenced_samples()
    refseq = patient.get_reference('genomewide', 'gb')

    for fea in refseq.features:
        if fea.id == stname:
            break

    start = fea.location.nofuzzy_start
    end = fea.location.nofuzzy_end

    fea_fragment = get_fragment_roi(refseq, start, end)
    fragment = fea_fragment.id
    start -= fea_fragment.location.nofuzzy_start
    end -= fea_fragment.location.nofuzzy_start

    (hct, ind, seqs) = patient.get_local_haplotype_count_trajectories(fragment,
                                                                      start, end,
                                                                      filters=['noN'],
                                                                      VERBOSE=VERBOSE)

    # Select only haplotypes that are seen at least a few times
    ind_abu = (hct > 5).any(axis=0)
    hct = hct[:, ind_abu]
    seqs = seqs[ind_abu]
    
    hft = (1.0 * hct.T / hct.sum(axis=1)).T

    # Initial consensus
    i0 = hct[0].argmax()
    seq0 = seqs[i0]
    
    if use_plot:
        import matplotlib.pyplot as plt
        from hivwholeseq.patients.get_local_haplotypes import plot_haplotype_frequencies

        def on_click(event):
            '''Print sequence on click'''
            mouseevent = event.mouseevent
            artist = event.artist
            i_clicked = int(artist.get_label())
            (score, ali1, ali2) = align_global(seq0, seqs[i_clicked], score_gapopen=-20)
            pretty_print_pairwise_ali((ali1, ali2), name1='cons0', name2='clicked', width=120)


        (fig, ax) = plot_haplotype_frequencies(patient.times[ind], hft,
                                               title=patient.name+', '+stname,
                                               picker=0.1)

        fig.canvas.mpl_connect('pick_event', on_click)


        plt.ion()
        plt.show()


    # Predict RNA structures
    structs = []
    for i, seq in enumerate(seqs):
        if VERBOSE >= 2:
            print 'Predicting structure n', i+1, 'of', len(seqs)
        label = stname+'_'+str(i)
        struct = predict_RNA_structure(seqs[i], label, maxstructs=1, VERBOSE=VERBOSE)[0]
        structs.append(struct)

    # Plot circle comparisons
    struct_cons0 = structs[i0]
    ind_minor = sorted(set((hft > 0.01).any(axis=0).nonzero()[0]) - set([i0]))

    filename_ct0 = struct_cons0['filename']
    for i in ind_minor:
        filename_ct = structs[i]['filename']
        fn_out = plot_circle_compare(filename_ct, filename_ct0, VERBOSE=VERBOSE)

