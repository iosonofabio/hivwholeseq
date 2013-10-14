# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/08/13
content:    Check sequencing/PCR erros in pure plasmid HIV samples (they appear
            as minor variants).
'''
# FIXME: modernize!
# TODO: check allele frequencies in overlapping regions
# Modules
import sys
from itertools import izip
import numpy as np
import matplotlib.pyplot as plt
from map_HIV_HXB2 import load_adapter_table
from minor_allele_frequency import get_minor_allele_counts


# Globals
VERBOSE = 1
alpha = np.array(list('ACGT-N'), 'S1')
read_types = ['read1 f', 'read1 r', 'read2 f', 'read2 r']
ref_filename = 'consensus_filtered_trimmed.fasta'
allele_count_filename = 'allele_counts.npy'



# Script
if __name__ == '__main__':

    # Input args
    args = sys.argv
    if len(args) < 2:
        raise ValueError('Please select an adapterID')
    adaID = int(args[1])

    # Get data
    data_folder = ('/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/')
    data_adaID_folder = data_folder+'adapterID_'+'{:02d}'.format(adaID)+'/'

    # Get adapter table
    adapter_table = load_adapter_table(data_folder)
    sample_name = adapter_table['sample'][adapter_table['ID'] == adaID][0]

    ## Read reference
    #ref_file = data_adaID_folder+ref_filename
    #if os.path.isfile(data_adaID_folder+ref_filename): ref_file = data_adaID_folder+ref_filename
    #else: ref_file = '/'.join(data_adaID_folder.split('/')[:-2]+['subsample/']+data_adaID_folder.split('/')[-2:])+ref_filename
    #refseq = SeqIO.read(ref_file, 'fasta')
    #ref = np.array(refseq)

    # Get counts
    # The first dimension is read type, the second alphabet, the third position
    counts = np.load(data_adaID_folder+allele_count_filename)
    
    # Define sequenced region
    sequenced_region = [(counts.sum(axis=1)[2] > 9000).nonzero()[0][0] - 1,
                        (counts.sum(axis=1)[3] > 9000).nonzero()[0][-1] + 2]
    counts = counts[:, :, sequenced_region[0]: sequenced_region[1]]

    # Get minor allele frequencies
    coverage = counts.sum(axis=1)
    major = counts.argmax(axis=1)
    nu_major = 1.0 * counts.max(axis=1) / (coverage + 0.0000001)
    # Top minor allele
    minor, counts_minor = get_minor_allele_counts(counts)
    nu_minor = 1.0 * counts_minor / (coverage + 0.0000001)

    # Plot coverage, genome profile, cumulative error histogram
    import matplotlib.cm as cm
    fig, axs = plt.subplots(3, 1, figsize=(14, 14))
    for i, read_type, nu_m in izip(xrange(len(read_types)), read_types, nu_minor):
        axs[0].plot(sequenced_region[0] + np.arange(coverage.shape[1]), coverage[i],
                    c=cm.jet(int(255.0 * i / len(read_types))))
        axs[1].plot(sequenced_region[0] + np.arange(coverage.shape[1]), nu_m + 0.0000001,
                    lw=1.5, c=cm.jet(int(255.0 * i / len(read_types))))
        ys = np.sort(nu_m[nu_m > 0])
        axs[2].plot(ys, 1 + 1.0 / len(ys) - np.linspace(0, 1, len(ys)), lw=2,
                    label=read_type,
                    c=cm.jet(int(255.0 * i / len(read_types))))

    depth = ys[int(0.99 * len(ys))]
    #depth = 2e-3
    axs[2].plot([depth] * 2, [1e-4, 1.15], lw=2, ls='--', c='k')
    axs[2].text(depth * 1.5, 0.15, 'approximate\nsequencing depth')

    axs[0].set_ylabel('Coverage [# reads]')
    axs[1].set_ylim(1e-5, 1e0)
    axs[1].set_yscale('log')
    axs[2].set_xscale('log')
    axs[2].set_yscale('log')
    axs[2].set_ylim(1e-4, 1.15)
    axs[1].set_xlabel('HIV genome [b.p.]')
    axs[1].set_ylabel('Minor allele frequency')
    axs[2].set_xlabel('Minor allele frequency')
    axs[2].set_ylabel('fraction of alleles with\nfrequency more than x')
    axs[2].legend(loc='lower left')
    axs[0].set_title('Sample: '+sample_name)

    plt.tight_layout()

    plt.savefig(sample_name+'_minor_alleles.pdf')

    plt.ion()
    plt.show()
