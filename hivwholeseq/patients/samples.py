# vim: fdm=marker
'''
author:     Fabio Zanini
date:       07/10/14
content:    Description module for HIV patient samples.
'''
# Modules
import numpy as np
import pandas as pd

from hivwholeseq.sequencing.filenames import table_filename



# Classes
class SamplePat(pd.Series):
    '''Patient sample'''


    def __init__(self, *args, **kwargs):
        '''Initialize a patient sample'''
        super(SamplePat, self).__init__(*args, **kwargs)
        self._sequenced_samples = None


    @property
    def _constructor(self):
        return SamplePat


    def get_n_templates_dilutions(self):
        '''Get the time course of the number of templates to PCR, limiting depth'''
        from hivwholeseq.patients.get_template_number import get_template_number
        return get_template_number(self.dilutions)


    def get_foldername(self, PCR=1):
        '''Get the name of the folder with the data for this sample'''
        from hivwholeseq.patients.filenames import get_sample_foldername
        return get_sample_foldername(self.patient, self.name, PCR=PCR)


    @property
    def samples_seq(self):
        '''Sequencing samples that refer to this patient sample'''
        if self._sequenced_samples is None:

            #TODO: optimize this call
            from hivwholeseq.patients.filenames import get_mapped_to_initial_filename
            from hivwholeseq.sequencing.samples import load_samples_sequenced as lss
            samples_seq = lss()
            samples_seq = samples_seq.loc[samples_seq['patient sample'] == self.name]
            self._sequenced_samples = samples_seq

        return self._sequenced_samples.copy()


    @samples_seq.setter
    def samples_seq(self, value):
        self._sequenced_samples = value


    def get_fragmented_roi(self, roi, VERBOSE=0, **kwargs):
        '''Get a region of interest in fragment coordinates'''
        from hivwholeseq.patients.get_roi import get_fragmented_roi
        if isinstance(roi, basestring):
            roi = (roi, 0, '+oo')
        refseq = self.get_reference('genomewide', 'gb')
        return get_fragmented_roi(refseq, roi, VERBOSE=VERBOSE, **kwargs)


    def get_fragments_covered(self, roi, VERBOSE=0, include_coordinates=False):
        '''Get fragments covering a certain roi'''
        from .get_roi import get_fragments_covered
        return get_fragments_covered(self, roi, VERBOSE=VERBOSE,
                                     include_coordinates=include_coordinates)


    def get_reference_filename(self, fragment, format='fasta'):
        '''Get filename of the reference for mapping'''
        from hivwholeseq.patients.filenames import get_initial_reference_filename
        return get_initial_reference_filename(self.patient, fragment, format)


    def get_reference(self, region, format='fasta'):
        '''Get the reference for a genomic region'''
        from Bio import SeqIO
        fragments = ['F'+str(i) for i in xrange(1, 7)] + ['genomewide']

        if region in fragments:
            fragment = region
        else:
            (fragment, start, end) = self.get_fragmented_roi((region, 0, '+oo'),
                                                             include_genomewide=True)

        refseq = SeqIO.read(self.get_reference_filename(fragment, format=format),
                            format)

        if format in ('gb', 'genbank'):
            from hivwholeseq.utils.sequence import correct_genbank_features_load
            correct_genbank_features_load(refseq)

        if region not in fragments:
            refseq = refseq[start: end]

        return refseq


    def get_mapped_filenames(self, fragment, PCR=1):
        '''Get filename(s) of mapped and filtered reads'''
        samples_seq = self.samples_seq
        fns = [get_mapped_to_initial_filename(self.patient, self.name, samplename,
                                              PCR=PCR)
               for samplename, sample in samples_seq.iterrows()]
        return fns


    def get_mapped_filtered_filename(self, fragment, PCR=1, **kwargs):
        '''Get filename(s) of mapped and filtered reads'''
        from hivwholeseq.patients.filenames import get_mapped_filtered_filename
        return get_mapped_filtered_filename(self.patient, self.name, fragment,
                                            PCR=PCR, **kwargs)


    def get_number_reads(self, fragment, PCR=1, **kwargs):
        '''Get the number of reads (not pairs)'''
        from ..utils.mapping import get_number_reads
        return get_number_reads(self.get_mapped_filtered_filename(fragment,
                                                                  PCR=PCR,
                                                                  **kwargs))


    def get_allele_counts_filename(self, fragment, PCR=1, qual_min=30, type='nuc'):
        '''Get the filename of the allele counts'''
        from hivwholeseq.patients.filenames import get_allele_counts_filename
        return get_allele_counts_filename(self.patient, self.name, fragment,
                                          PCR=PCR, qual_min=qual_min, type=type)


    def get_allele_cocounts_filename(self, fragment, PCR=1, qual_min=30,
                                     compressed=True):
        '''Get the filename of the allele counts'''
        from hivwholeseq.patients.filenames import get_allele_cocounts_filename
        return get_allele_cocounts_filename(self.patient, self.name, fragment,
                                            PCR=PCR, qual_min=qual_min,
                                            compressed=compressed)


    def get_consensus_filename(self, fragment, PCR=1):
        '''Get the filename of the consensus of this sample'''
        from hivwholeseq.patients.filenames import get_consensus_filename
        return get_consensus_filename(self.patient, self.name, fragment, PCR=PCR)


    def get_consensus(self, fragment, PCR=1):
        '''Get consensu for this sample'''
        from Bio import SeqIO
        return SeqIO.read(self.get_consensus_filename(fragment, PCR=PCR), 'fasta')


    # TODO: implement this with depth checks
    #def get_allele_frequencies(self, region, PCR=1, qual_min=30, merge_read_types=True):
    #    pass


    def get_allele_counts(self, region, PCR=1, qual_min=30, merge_read_types=True):
        '''Get the allele counts'''
        import numpy as np

        # Fall back on genomewide counts if no single fragment is enough
        (fragment, start, end) = self.get_fragmented_roi((region, 0, '+oo'),
                                                         include_genomewide=True)

        ac = np.load(self.get_allele_counts_filename(fragment, PCR=PCR, qual_min=qual_min))
        ac = ac[:, :, start: end]

        if merge_read_types:
            ac = ac.sum(axis=0)
        return ac


    def get_allele_counts_aa(self, protein, PCR=1, qual_min=30):
        '''Get the amino acid allele counts'''
        import numpy as np
        ac = np.load(self.get_allele_counts_filename(protein, PCR=PCR,
                                                     qual_min=qual_min,
                                                     type='aa'))
        return ac


    def get_allele_cocounts(self, fragment, PCR=1, qual_min=30):
        '''Get the allele cocounts'''
        import numpy as np
        acc = np.load(self.get_allele_cocounts_filename(fragment, PCR=PCR,
                                                        qual_min=qual_min))['cocounts']
        return acc


    def get_coverage(self, region, PCR=1, qual_min=30, merge_read_types=True):
        '''Get the coverage'''
        ac = self.get_allele_counts(region, PCR=PCR, qual_min=qual_min,
                                    merge_read_types=merge_read_types)
        cov = ac.sum(axis=-2)
        return cov


    def get_local_haplotypes(self,
                             fragment, start, end,
                             VERBOSE=0,
                             maxreads=-1,
                             filters=None,
                             PCR=1):
        '''Get local haplotypes'''
        from hivwholeseq.patients.get_local_haplotypes import get_local_haplotypes
        bamfilename = self.get_mapped_filtered_filename(fragment, PCR=PCR)
        haplo = get_local_haplotypes(bamfilename,
                                     start, end,
                                     VERBOSE=VERBOSE,
                                     maxreads=maxreads,
                                     label=self.name)

        if filters is not None:
            if 'noN' in filters:
                hnames = [hname for hname in haplo.iterkeys() if 'N' in hname]
                for hname in hnames:
                    del haplo[hname]

            if 'nosingletons' in filters:
                hnames = [hname for hname, c in haplo.iteritems() if c <= 1]
                for hname in hnames:
                    del haplo[hname]

            if any('mincount=' in ft for ft in filters):
                for ft in filters:
                    if 'mincount=' in ft:
                        break
                cmin = int(ft[len('mincount='):])
                hnames = [hname for hname, c in haplo.iteritems() if c < cmin]
                for hname in hnames:
                    del haplo[hname]
            
            if any('freqmin=' in ft for ft in filters):
                for ft in filters:
                    if 'freqmin=' in ft:
                        break
                fmin = float(ft[len('minfreq='):])
                csum = sum(haplo.itervalues())
                hnames = [hname for hname, c in haplo.iteritems() if c < fmin * csum]
                for hname in hnames:
                    del haplo[hname]

        return haplo



# Functions
def itersample(samples):
    for samplename, sample in samples.iterrows():
        yield (samplename, SamplePat(sample))


def load_samples_sequenced(patients=None, include_empty=False):
    '''Load patient samples sequenced from general table'''
    sample_table = pd.read_excel(table_filename, 'Samples timeline sequenced',
                                 index_col=0)

    # Reindex DataFrame
    sample_table.index = pd.Index(map(str, sample_table.index))
    sample_table.loc[:, 'patient'] = map(str, sample_table.loc[:, 'patient'])

    # Note: this refers to the TOTAL # of templates, i.e. the factor 2x for
    # the two parallel RT-PCR reactions
    sample_table['n templates'] = sample_table['viral load'] * 0.4 / 12 * 2

    if not include_empty:
        ind = [i for i, sample in sample_table.iterrows()
               if (sample[['F1', 'F2', 'F3', 'F4', 'F5', 'F6']] != 'miss').any()]
        sample_table = sample_table.loc[ind]

    if patients is not None:
        sample_table = sample_table.loc[sample_table.loc[:, 'patient'].isin(patients)]

    return sample_table


def load_sample_sequenced(samplename):
    '''Load one patient sample'''
    return SamplePat(load_samples_sequenced().loc[samplename])

