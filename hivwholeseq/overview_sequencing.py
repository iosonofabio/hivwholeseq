# vim: fdm=marker
'''
author:     Fabio Zanini
date:       29/08/14
content:    Command-line overview of the whole project, the sequencing part.
'''
# Modules
from collections import Counter
from itertools import izip
from operator import itemgetter
from hivwholeseq.samples import load_sequencing_runs, SequencingRun, \
        load_samples_sequenced

# Functions
def line_format(seq_run, n_samples_pat):
    '''Format description line for sequencing run'''
    data = [seq_run.name,
            seq_run.library,
            seq_run.plex,
            n_samples_pat[seq_run.name]]
    
    fields = []
    for (datum, form) in izip(data, col_formats):
        fields.append(('{:>'+str(form[1])+form[2]+'}').format(data[len(fields)]))
    
    line = ' | '.join(fields)
    return line



# Script
if __name__ == '__main__':

    print '----------------------------------------'
    print 'HIV whole genome longitudinal sequencing'
    print '----------------------------------------'
    
    seq_runs = load_sequencing_runs()
    samples = load_samples_sequenced()
    print 'Sequencing runs:', ', '.join(seq_runs.index)

    n_samples_pat = Counter(samples.loc[samples['patient sample'] != 'nan']['seq run'])
    linef = lambda x: line_format(x, n_samples_pat)
    col_formats = [('run', 15, 's'),
                   ('library', 20, 's'),
                   ('# samples', 10, 'd'),
                   ('# pat samples', 10, 'd')]

    header = ' | '.join([('{:^'+str(cf)+'s}').format(cn) for (cn, cf, _) in col_formats])
    head_line = '-' * len(header)

    print head_line
    print header
    print head_line
    for seq_runname, seq_run in seq_runs.iterrows():
        print linef(seq_run)

