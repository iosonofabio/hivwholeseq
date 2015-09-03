# vim: fdm=marker
'''
author:     Fabio Zanini
date:       03/09/15
content:    Make symlinks to the sequencing samples.
'''
# Modules
import sys
import os

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.samples import itersample
from hivwholeseq.sequencing.samples import load_samples_sequenced as lss
from hivwholeseq.patients.samples import load_samples_sequenced as lssp
from hivwholeseq.sequencing.filenames import get_sample_foldername



# Script
if __name__ == '__main__':

    samples_pat = lssp()
    samples_seq = lss()

    for samplename, sample in itersample(samples_pat):
        root_foldername = sample.get_foldername()+'samples_sequencing/'
        mkdirs(root_foldername)

        for samplenameseq, sampleseq in samples_seq.iterrows():
            if sampleseq['patient sample'] == samplename:
                src_folder = get_sample_foldername(samplenameseq)
                dst_folder = root_foldername+samplenameseq
                if not os.path.islink(dst_folder):
                    os.symlink(src_folder, dst_folder)
                    print 'Symlink:', src_folder, dst_folder
                else:
                    print 'Esists:', dst_folder
