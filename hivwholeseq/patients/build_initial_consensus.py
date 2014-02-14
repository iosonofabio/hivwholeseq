# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/10/13
content:    Update the initial consensus for the paitent (when we get more data).
'''
# TODO: if the first sequenced sample has no fragment 5, we are in trouble
# In principle, we should use the first *sequenced* F5, but in practice that is
# not in line with the other fragments... wait until that happens before getting
# dirty
# Modules
import os
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.miseq import alpha
from hivwholeseq.filenames import get_consensus_filename
from hivwholeseq.patients.filenames import get_initial_consensus_filename, \
        get_foldername
from hivwholeseq.patients.patients import get_patient
from hivwholeseq.samples import samples, date_to_integer
from hivwholeseq.primer_info import primers_inner
from hivwholeseq.primer_info import primers_coordinates_HXB2_inner as pci



# Functions
def expand_consensus_F5(patient, seq, sample_init,
                        data_folder_init, adaID_init,
                        mismatches_max=2,
                        VERBOSE=0):
    '''Expand consensus of fragment F5 to the union of all primer sets (trimmed)'''
    # Check what primers are used at all in this patient
    samples_seq = get_sequenced_samples(patient)
    primer_set = set()
    for sample in samples_seq:
        dataset = MiSeq_runs[samples[sample]['run']]
        prF5 = dataset['primerF5'][dataset['samples'].index(sample)]
        primer_set.add(prF5)

    # If only one set of primers is used throughout, we do nothing
    if len(primer_set) == 1:
        return seq

    # else, we must pick from neighboring F4 and F6 the missing edges
    # NOTE: as long as only two sets of primers are present, only one of the two
    # edges must be tampered, because they are not nested
    dataset = MiSeq_runs[samples[sample_init]['run']]
    prF5 = dataset['primerF5'][dataset['samples'].index(sample_init)]

    # F4: find the first F5 primer, the start of seq, and fill the gap
    prF5_first = min(primer_set, key=lambda x: pci[x][0][1])
    if prF5_first != prF5:
        seq_F4 = SeqIO.read(get_consensus_filename(data_folder_init, adaID_init, 'F4'),
                            'fasta')
        seq_F4 = np.array(seq_F4[-500:])

        # Find the fwd inner primer end
        prF5_fwd = primers_inner[prF5_first][0]
        prF5_fwd = np.ma.array(np.fromstring(prF5_fwd, 'S1'),
                               mask=[a not in alpha[:4] for a in prF5_fwd])
        for pos in xrange(len(seq_F4) - len(prF5_fwd)):
            n_matches = (seq_F4[pos: pos + len(prF5_fwd)] == prF5_fwd).sum()
            if n_matches >= (-prF5_fwd.mask).sum() - mismatches_max:
                break
        pos_F5_fwd = pos + len(prF5_fwd)

        # Find the beginnning of seq
        seed_len = 20
        seed = np.array(seq[:seed_len], 'S1')
        for pos in xrange(len(seq_F4) - len(seed)):
            n_matches = (seq_F4[pos: pos + len(seed)] == seed).sum()
            if n_matches >= seed_len - mismatches_max:
                break
        pos_seed = pos

        seq = SeqRecord(Seq(seq_F4[pos_F5_fwd: pos_seed].tostring()+str(seq.seq),
                            seq.seq.alphabet),
                        id=seq.id, name=seq.name,
                        description=seq.description)

        if VERBOSE >= 1:
            print 'Expanded F5 from F4'
    
    # F6: same difference
    # FIXME: this parti is still untested
    prF5_last = max(primer_set, key=lambda x: pci[x][1][0])
    if prF5_last != prF5:
        seq_F6 = SeqIO.read(get_consensus_filename(data_folder_init, adaID_init, 'F6'),
                            'fasta')
        seq_F6 = np.array(seq_F6[:500])

        # Find the rev inner primer start
        prF5_rev = primers_inner[prF5_last][1]	# already rev comp --->
        prF5_rev = np.ma.array(np.fromstring(prF5_rev, 'S1'),
                               mask=[a not in alpha[:4] for a in prF5_rev])
        for pos in xrange(len(seq_F6) - len(prF5_rev)):
            n_matches = (seq_F6[pos: pos + len(prF5_rev)] == prF5_rev).sum()
            if n_matches >= (-prF5_rev.mask).sum() - mismatches_max:
                break
        pos_F5_rev = pos

        # Find the end of seq
        seed_len = 20
        seed = np.array(seq[-seed_len:], 'S1')
        for pos in xrange(len(seq_F6) - len(seed)):
            n_matches = (seq_F6[pos: pos + len(seed)] == seed).sum()
            if n_matches >= seed_len - mismatches_max:
                break
        pos_seed = pos

        seq = SeqRecord(Seq(str(seq.seq)+seq_F6[pos_seed: pos_F5_rev].tostring(),
                            seq.seq.alphabet),
                        id=seq.id, name=seq.name,
                        description=seq.description)

        if VERBOSE >= 1:
            print 'Expanded F5 from F6'
 
    return seq



# Script
if __name__ == '__main__':

    # Parse input args
    parser = argparse.ArgumentParser(description='Update initial consensus')
    parser.add_argument('--patient', required=True,
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='*',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')

    args = parser.parse_args()
    pname = args.patient
    fragments = args.fragments
    VERBOSE = args.verbose

    # Get the patient
    patient = get_patient(pname)

    # If the script is called with no fragment, iterate over all
    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    # Make dir for the patient if absent
    pfolder = get_foldername(pname)
    if not os.path.isdir(pfolder):
        os.mkdir(pfolder)
        if VERBOSE >= 1:
            print pname+': folder created.'
    
    # Get the first sequenced sample
    sample_init = patient.initial_sample
    seq_run = sample_init['run']
    adaID = sample_init['adaID']
    dataset = MiSeq_runs[seq_run]
    data_folder = dataset['folder']

    # Check for the existence of an initial consensus already
    for fragment in fragments:

        # Read the new consensus
        input_filename = get_consensus_filename(data_folder, adaID, fragment)
        seq_in = SeqIO.read(input_filename, 'fasta')

        # FIXME: IGNORE F5 PECULIARITY FOR THE TIME BEING
        ## Fragment F5 needs a special treatment, because of the double primers
        #if fragment == 'F5':
        #    seq_in = expand_consensus_F5(patient, seq_in, sample_init,
        #                                 data_folder, adaID,
        #                                 VERBOSE=VERBOSE)

        # Write output
        output_filename = get_initial_consensus_filename(pname, fragment)

        # If absent, just copy the thing over
        if not os.path.isfile(output_filename):
            SeqIO.write(seq_in, output_filename, 'fasta')
            if VERBOSE >= 1:
                print pname+': initial consensus file created for sample', sample_init

        # if present, check whether the sequences are the same (if so, no remapping
        # is needed!). Overwrite the file anyway, because single samples carry
        # their consensus (mapping reference) with them in the folder (not much
        # overhead and MUUUCH cleaner than otherwise).
        else:
            seq_out = SeqIO.read(output_filename, 'fasta')
            SeqIO.write(seq_in, output_filename, 'fasta')
            if str(seq_in.seq) != str(seq_out.seq):
                print 'NOTE: initial consensus updated to '+sample_init+', remap!'
            
