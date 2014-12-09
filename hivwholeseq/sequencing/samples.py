# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/10/13
content:    Description module for HIV samples from patients.

            Note: it is possible to read the Excel XLSX table with xlrd. It also
            works ok, but it's painful and error-prone code, and it assumes that
            table is up-to-date. Better write down the machine-relevant info
            manually down here (until otherwise specified).
'''
# Modules
import numpy as np
import pandas as pd

from hivwholeseq.sequencing.filenames import table_filename



# Globals
_sequencing_runs = None



# Classes
class SampleSeq(pd.Series):
    '''A sequenced sample (if something has been sequenced twice, they are separate)'''

    _sequencing_run = None
    _sample_pat = None


    def __init__(self, *args, **kwargs):
        '''Initialize a sequenced sample'''
        super(SampleSeq, self).__init__(*args, **kwargs)

        from hivwholeseq.sequencing.filenames import get_seqrun_foldername
        from hivwholeseq.sequencing.adapter_info import foldername_adapter
        seq_run = self.loc['seq run']
        adaID = self.loc['adapter']
        self['folder'] = str(get_seqrun_foldername(seq_run)+foldername_adapter(adaID))
        self['seqrun_folder'] = str(get_seqrun_foldername(seq_run))


    @property
    def _constructor(self):
        return SampleSeq


    @property
    def sequencing_run(self):
        '''The sequencing run object of this sample'''
        if self._sequencing_run is None:
            self._sequencing_run = load_sequencing_run(self.loc['seq run'])
        return self._sequencing_run


    @property
    def sample_pat(self):
        '''Patient sample of this sequencing sample'''
        if self._sample_pat is None:
            from hivwholeseq.patients.patients import load_samples_sequenced as lssp
            from hivwholeseq.patients.patients import SamplePat
            self._sample_pat = SamplePat(lssp(include_wrong=True).loc[self['patient sample']])
        return self._sample_pat


    @property
    def patientname(self):
        '''Name of the patient this sample belongs to'''
        return self.sample_pat.patient


    @property
    def regions_complete(self):
        '''Get the complete regions, e.g. F5ao'''
        if self.PCR == 1:
            PCR_suffix = 'o'
        elif self.PCR == 2:
            PCR_suffix = 'i'
        else:
            PCR_suffix = ''
        
        regions = self.regions.split(' ')
        if len(regions):
            return ['F'+fr+PCR_suffix for fr in regions]
        else:
            return []


    @property
    def regions_generic(self):
        '''Get the complete regions, e.g. F5ao'''        
        regions = self.regions.split(' ')
        if len(regions):
            return ['F'+fr[0] for fr in regions]
        else:
            return []


    def convert_region(self, fragment):
        '''Get the complete region of a generic or vice versa'''
        if fragment in self.regions_generic:
            return self.regions_complete[self.regions_generic.index(fragment)]
        elif fragment in self.regions_complete:
            return self.regions_generic[self.regions_complete.index(fragment)]
        else:
            raise IndexError('Region not in complete nor in generic regions!')


    def get_read_filenames(self, **kwargs):
        '''Get the filenames of the demultiplexed reads'''
        from hivwholeseq.sequencing.filenames import get_read_filenames as gfn
        return gfn(self.folder, **kwargs)


    def get_reference_premap_filename(self, **kwargs):
        '''Get the filename of the premapping refernce'''
        from hivwholeseq.sequencing.filenames import get_reference_premap_filename as gfn
        return gfn(self.seqrun_folder, self.adapter, **kwargs)


    def get_premapped_filename(self, **kwargs):
        '''Get the filename of the readed premapped to reference'''
        from hivwholeseq.sequencing.filenames import get_premapped_filename as gfn
        return gfn(self.folder, **kwargs)


    def get_divided_filename(self, fragment, **kwargs):
        '''Get the filename of the divided and trimmed reads'''
        from hivwholeseq.sequencing.filenames import get_divided_filename as gfn
        return gfn(self.folder, adaID=None, fragment=fragment, **kwargs)


    def get_consensus_filename(self, fragment, **kwargs):
        '''Get the filename of the consensus'''
        from hivwholeseq.sequencing.filenames import get_consensus_filename as gfn
        from hivwholeseq.sequencing.filenames import get_merged_consensus_filename as gfn2
        if fragment != 'genomewide':
            return gfn(self.folder, adaID=None, fragment=fragment, **kwargs)
        else:
            return gfn2(self.folder, adaID=None, **kwargs)


    def get_mapped_filename(self, fragment, **kwargs):
        '''Get the filename of the mapped reads'''
        from hivwholeseq.sequencing.filenames import get_mapped_filename as gfn
        return gfn(self.folder, adaID=None, fragment=fragment, **kwargs)


    def get_mapped_to_initial_filename(self, fragment, **kwargs):
        '''Get the filename of the reads mapped to patient initial reference'''
        from hivwholeseq.patients.filenames import get_mapped_to_initial_filename as gfn
        pname = self.patientname
        samplename_pat = self['patient sample']
        PCR = int(self.PCR)
        return gfn(pname, samplename_pat, self.name, fragment, PCR=PCR, **kwargs)


    def get_mapped_filtered_filename(self, fragment, **kwargs):
        '''Get the filename of the reads mapped to patient initial reference'''
        from hivwholeseq.patients.filenames import get_mapped_filtered_filename as gfn
        pname = self.patientname
        samplename_pat = self['patient sample']
        PCR = int(self.PCR)
        return gfn(pname, samplename_pat, fragment, PCR=PCR, **kwargs)


    def get_consensus(self, fragment):
        '''Get consensus sequence for mapping'''
        from Bio import SeqIO
        return SeqIO.read(self.get_consensus_filename(fragment), 'fasta')


    def get_allele_counts_filename(self, fragment):
        '''Get the filename with the allele counts'''
        from .filenames import get_allele_counts_filename
        return get_allele_counts_filename(self.seqrun_folder, self.adapter, fragment)


    def get_allele_counts(self, fragment, merge_read_types=False):
        '''Get the allele counts'''
        import numpy as np
        counts = np.load(self.get_allele_counts_filename(fragment))
        if merge_read_types:
            counts = counts.sum(axis=0)
        return counts


    def get_insert_size_distribution(self, fragment, VERBOSE=0, **kwargs):
        '''Get insert size distribution for a fragment'''
        from .check_insert_distribution import get_insert_size_distribution as gisd
        return gisd(self.seqrun_folder, self.adapter, fragment, VERBOSE=VERBOSE,
                    **kwargs)



class SamplesSeq(pd.DataFrame):
    '''Table of sequenced samples'''

    @property
    def _constructor(self):
        return SamplesSeq


    def filter_seq_run(self, run):
        '''Get only the samples from one sequencing run'''
        return self.loc[self.run == run]


    def filter_patient(self, patient_id, exclude_reps=True):
        '''Get only the samples from a specific patient'''
        samples = self.loc[self.patient == patient_id]
        if exclude_reps:
            ind = [(len(sn) < 2) or (sn[-2:] not in ('-2', '-3')) for sn in samples.name.tolist()]
            samples = samples[ind]
        return samples


class SequencingRun(pd.Series):
    '''Sequencing run'''

    def __init__(self, *args, **kwargs):
        '''Initialize a sequencing run'''
        super(SequencingRun, self).__init__(*args, **kwargs)

        from hivwholeseq.sequencing.filenames import get_seqrun_foldername
        self['folder'] = str(get_seqrun_foldername(self.name))

        self['samples'] = load_samples_sequenced(seq_runs=[self.name])


    @property
    def _constructor(self):
        return SequencingRun


    def itersamples(self):
        '''Generator for samples in this run, each with extended attributes'''
        for samplename, sample in self.samples.iterrows():
            yield SampleSeq(sample)


    @property
    def adapters(self):
        '''The adapters used in the sequencing run'''
        return self.samples.adapter


    def discard_nondivided_samples(self):
        '''Discard samples that have no divided reads (e.g. SA, random hexamers)'''
        import os
        from hivwholeseq.sequencing.filenames import get_divided_filename
        ind = []
        for sample in self.itersamples():
            frag = sample.regions_complete[0]
            div = os.path.isfile(get_divided_filename(self.folder, sample.adapter, frag))
            ind.append(div)
        self.samples = self.samples.loc[ind]



# Functions
def load_samples_sequenced(seq_runs=None):
    '''Load samples sequenced from general table'''
    sample_table = pd.read_excel(table_filename, 'Samples sequenced',
                                 index_col=0)
    sample_table.index = pd.Index(map(str, sample_table.index))
    sample_table.loc[:, 'patient sample'] = map(str, sample_table.loc[:, 'patient sample'])
    sample_table.loc[:, 'regions'] = map(str, sample_table.loc[:, 'regions'])
    sample_table.loc[sample_table.loc[:, 'regions'] == 'nan', 'regions'] = ''

    if seq_runs is not None:
        sample_table = sample_table.loc[sample_table.loc[:, 'seq run'].isin(seq_runs)]

    return SamplesSeq(sample_table)


def load_sample_sequenced(samplename):
    '''Load a sequenced sample from the general table'''
    return SampleSeq(load_samples_sequenced().loc[samplename])


def load_sequencing_runs(seq_runs=None):
    '''Load sequencing runs from general table'''
    global _sequencing_runs
    if _sequencing_runs is None:
        seq_runs_in = pd.read_excel(table_filename, 'Sequencing runs',
                                 index_col=0)
        _sequencing_runs = seq_runs_in

    if seq_runs is not None:
        seq_runs = _sequencing_runs.loc[_sequencing_runs.index.isin(seq_runs)]
        return seq_runs
    else:
        return _sequencing_runs


def load_sequencing_run(seq_run):
    '''Load a single sequencing run (with extended attributes)'''
    seq_runs = load_sequencing_runs()
    seq_run_obj = SequencingRun(seq_runs.loc[seq_run])
    return seq_run_obj

