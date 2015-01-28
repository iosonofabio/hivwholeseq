# vim: fdm=marker
'''
author:     Fabio Zanini
date:       09/05/14
content:    The run Tuen3 was a 14 plex. We want to resequence the same with
            different sample proportions in order to fill gaps.
'''
# Modules
from itertools import izip

from hivwholeseq.datasets import MiSeq_runs
from hivwholeseq.utils.mapping import get_number_mapped_reads
from hivwholeseq.sequencing.filenames import get_premapped_filename


# Globals
vols = [14.8,
        7.27,
        5.88,
        5.84,
        14.68,
        13.01,
        15.24,
        8.89,
        7.48,
        6.61,
        5.80,
        14.95,
        17.02,
        9.58]

concs = [8.0 / vol for vol in vols]



# Script
if __name__ == '__main__':

    dataset = MiSeq_runs['Tuen3']
    data_folder = dataset['folder']
    adaIDs = dataset['adapters']
    samplenames = dataset['samples']

    n_mapped_reads = []
    for (adaID, samplename) in izip(adaIDs, samplenames):
        bamfilename = get_premapped_filename(data_folder, adaID)
        n = get_number_mapped_reads(bamfilename)
        print adaID, samplename, n
        n_mapped_reads.append(n)

    n_missing_reads = [max(n_mapped_reads) - n for n in n_mapped_reads]
    reads_per_ul = [1.0 * n / vol for (n, vol) in izip(n_mapped_reads, vols)]
    ul_missing = [n / nul for (n, nul) in izip(n_missing_reads, reads_per_ul)]

    ul_left = [35 - vol for vol in vols]
    prop = min([ul / (um + 0.001) for (ul, um) in izip(ul_left, ul_missing)])
    ul_missing_capped = [um * prop for um in ul_missing]
    
    ng_tot = sum(conc * umc for (conc, umc) in izip(concs, ul_missing_capped))

    conc_tot = ng_tot / 30
    conc_tot_nM = conc_tot * 1e6 / 660 / 650

    print 'Sample name\tvol [ul]'
    for (samplename, umc) in izip(samplenames, ul_missing_capped):
        print samplename, '\t', umc
    print 'Final library molality:', conc_tot_nM, 'nM'
