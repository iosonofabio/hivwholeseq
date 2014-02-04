# vim: fdm=marker
'''
author:     Fabio Zanini
date:       30/01/14
content:    Functions for getting standard filenames.
'''
# Functions
def get_reference_premap_filename(data_folder, samplename):
    '''Get the filename of the reference used from premapping'''
    fn = 'reference.fasta'
    fn = data_folder+samplename+'/premapped/'+fn
    return fn


def get_premapped_file(data_folder, samplename, type='bam'):
    '''Get the filename of the readed mapped to reference to split into fragments'''
    filename = 'premapped'
    filename = data_folder+samplename+'/premapped/'+filename+'.'+type
    return filename

