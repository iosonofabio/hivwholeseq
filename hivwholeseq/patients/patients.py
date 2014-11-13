# vim: fdm=marker
'''
author:     Fabio Zanini
date:       23/10/13
content:    Description module for HIV patients.
'''
# Modules
import numpy as np
import pandas as pd

from hivwholeseq.sequencing.filenames import table_filename
from hivwholeseq.patients.samples import * # FIXME: lots of scripts import from here still



# Classes
class Patient(pd.Series):
    '''HIV patient'''

    def __init__(self, *args, **kwargs):
        '''Initialize a patient with all his samples'''
        super(Patient, self).__init__(*args, **kwargs)
        from hivwholeseq.patients.samples import load_samples_sequenced
        samples = load_samples_sequenced(patients=[self.name])
        self.samples = samples


    @property
    def _constructor(self):
        return Patient


    @property
    def folder(self):
        '''The folder with the data on this patient'''
        from hivwholeseq.patients.filenames import get_foldername
        return str(get_foldername(self.name))


    def discard_nonsequenced_samples(self):
        '''Discard all samples that have not been sequenced yet'''
        from hivwholeseq.sequencing.samples import load_samples_sequenced as lss
        samples_sequenced = lss()
        samples_sequenced_set = set(samples_sequenced.loc[:, 'patient sample']) - set(['nan'])
        samples = self.samples.loc[self.samples.index.isin(samples_sequenced_set)]

        ## Add info on sequencing
        ## FIXME: why is this here?!
        ## FIXME: this is buggy is so many ways... pandas is nto great at this
        #samples_seq_col = []
        #for samplename in samples.index:
        #    ind = samples_sequenced.loc[:, 'patient sample'] == samplename
        #    samples_seq_col.append(samples_sequenced.loc[ind])

        #samples.loc[:, 'samples seq'] = samples_seq_col

        self.samples = samples


    @property
    def dates(self):
        '''Get the dates of sampling'''
        return self.samples.date

    
    @property
    def times(self, unit='day'):
        '''Get the times from transmission'''
        return convert_date_deltas_to_float(self.dates - self.transmission_date, unit=unit)


    @property
    def viral_load(self):
        '''Get the time course of the viral load [molecules/ml of serum]'''
        return self.samples['viral load']


    @property
    def cell_count(self):
        '''Get the time course of the CD4+ cell count'''
        return self.samples['CD4+ count']


    @property
    def n_templates(self, n_reactions=6):
        '''Get the time course of the number of templates to PCR, limiting depth'''
        n = self.viral_load.copy()
        # We take 400 ul of serum
        n *= 0.4
        # We typically have 6 reactions with that total volume (plus the F4 dilution
        # series, but each of those uses only 0.1x template which is very little)
        n /= (n_reactions + 0.2)
        return n


    @property
    def initial_sample(self):
        '''The initial sample used as a mapping reference'''
        from .samples import SamplePat
        return SamplePat(self.samples.iloc[0])


    def itersamples(self):
        '''Generator for samples in this patient, each with extended attributes'''
        from hivwholeseq.patients.samples import SamplePat
        for samplename, sample in self.samples.iterrows():
            yield SamplePat(sample)


    def get_reference_filename(self, fragment, format='fasta'):
        '''Get filename of the reference for mapping'''
        from hivwholeseq.patients.filenames import get_initial_reference_filename
        return get_initial_reference_filename(self.name, fragment, format)


    def get_reference(self, fragment, format='fasta'):
        '''Get the reference for a fragment'''
        from Bio import SeqIO
        refseq = SeqIO.read(self.get_reference_filename(fragment, format=format), format)
        if format in ('gb', 'genbank'):
            from hivwholeseq.sequence_utils import correct_genbank_features_load
            correct_genbank_features_load(refseq)
        return refseq


    def get_consensi_alignment_filename(self, fragment):
        '''Get the multiple sequence alignment of all consensi'''
        from hivwholeseq.patients.filenames import get_consensi_alignment_filename
        return get_consensi_alignment_filename(self.name, fragment)


    def get_consensi_tree_filename(self, fragment):
        '''Get the filename of the consensi of the patient'''
        from hivwholeseq.patients.filenames import get_consensi_tree_filename
        return get_consensi_tree_filename(self.name, fragment)


    def get_initial_allele_counts(self, fragment):
        '''Get allele counts from the initial time point'''
        import os
        from hivwholeseq.patients.samples import SamplePat
        for i in xrange(len(self.samples)):
            sample = SamplePat(self.samples.iloc[i])
            if os.path.isfile(sample.get_allele_counts_filename(fragment)):
                return sample.get_allele_counts(fragment)


    def get_initial_allele_frequencies(self, fragment, cov_min=1):
        '''Get the allele frequencies from the initial time point'''
        counts = self.get_initial_allele_counts(fragment)
        cov = counts.sum(axis=0)
        af = np.ma.masked_where(np.tile(cov < cov_min, (counts.shape[0], 1)), counts)
        af.harden_mask()
        af = 1.0 * af / af.sum(axis=0)
        return af


    def get_coverage_trajectories(self, fragment, use_PCR1=1):
        '''Get coverage as a function of time'''
        (act, ind) = self.get_allele_count_trajectories(fragment, use_PCR1=use_PCR1)
        return (act.sum(axis=1), ind)


    def get_allele_frequency_trajectories(self, fragment, use_PCR1=1, cov_min=1,
                                          depth_min=None, **kwargs):
        '''Get the allele frequency trajectories from files
        
        Args:
          cov_min (int): minimal coverage accepted, anything lower are masked.
          depth_min (float): minimal depth, both by sequencing and template numbers.
            Time points with less templates are excluded, and positions are masked.
            For convenience depth is defined > 1, e.g. 100 takes frequencies down
            to 1%.
          **kwargs: passed down to the get_allele_count_trajectories method.
        '''
        (act, ind) = self.get_allele_count_trajectories(fragment,
                                                        use_PCR1=use_PCR1,
                                                        **kwargs)

        if depth_min is not None:
            indd = np.array(self.n_templates[ind] >= depth_min)
            act = act[indd]
            ind = ind[indd]
            cov_min = max(cov_min, depth_min)

        covt = act.sum(axis=1)
        mask = np.zeros_like(act, bool)
        mask.swapaxes(0, 1)[:] = covt < cov_min

        # NOTE: the hard mask is necessary to avoid unmasking part of the alphabet
        # at a certain site: the mask is site-wise, not allele-wise
        aft = np.ma.array((1.0 * act.swapaxes(0, 1) / covt).swapaxes(0, 1),
                          mask=mask,
                          hard_mask=True,
                          fill_value=0)

        aft[(aft < 1e-4)] = 0
        # NOTE: we'd need to renormalize, but it's a small effect

        return (aft, ind)


    def get_allele_count_trajectories(self, fragment, use_PCR1=1, **kwargs):
        '''Get the allele count trajectories from files
        
        Args:
          **kwargs: passed down to the function (VERBOSE, etc.).

        Note: the genomewide counts are currently saved to file.
        '''
        if fragment == 'genomewide':
            from hivwholeseq.patients.filenames import \
                    get_allele_count_trajectories_filename as get_fn
            fn = get_fn(self.name, fragment)
            npz = np.load(fn)
            act = npz['act']
            ind = npz['ind']

        else:
            from hivwholeseq.patients.one_site_statistics import \
                    get_allele_count_trajectories
            from operator import itemgetter
            (sns, act) = get_allele_count_trajectories(self.name, self.samples.index,
                                                       fragment,
                                                       use_PCR1=use_PCR1, **kwargs)
            ind = np.array([i for i, (_, sample) in enumerate(self.samples.iterrows())
                            if sample.name in map(itemgetter(0), sns)], int)

        return (act, ind)


    def get_mapped_filtered_filename(self, samplename, fragment, PCR=1):
        '''Get filename(s) of mapped and filtered reads for a sample'''
        from hivwholeseq.patients.filenames import get_mapped_filtered_filename
        return get_mapped_filtered_filename(self.patient, samplename, fragment, PCR=PCR)


    def get_divergence(self, fragment, **kwargs):
        '''Get divergence of a fragment
        
        Args:
          **kwargs: passed to the allele frequency trajectories.
        '''
        from hivwholeseq.patients.get_divergence_diversity import get_divergence
        aft, ind = self.get_allele_frequency_trajectories(fragment, **kwargs)
        return (get_divergence(aft), ind)


    def get_diversity(self, fragment, **kwargs):
        '''Get diversity of a fragment
        
        Args:
          **kwargs: passed to the allele frequency trajectories.
        '''
        from hivwholeseq.patients.get_divergence_diversity import get_diversity
        aft, ind = self.get_allele_frequency_trajectories(fragment, **kwargs)
        return (get_diversity(aft), ind)


    @property
    def transmission_date(self):
        '''The most likely time of transmission'''
        return self['last negative date'] + \
                (self['first positive date'] - self['last negative date']) / 2


    def get_map_coordinates_reference(self, fragment, refname='HXB2', roi=None):
        '''Get the map of coordinate to some external reference
        
        Returns:
          mapco (2D int array): the first column are the positions in the reference,
            the second column the position in the patient initial reference. 
        '''
        from hivwholeseq.patients.filenames import get_coordinate_map_filename
        fn = get_coordinate_map_filename(self.name, fragment, refname=refname)
        mapco = np.loadtxt(fn, dtype=int)
        if roi is None:
            return mapco
        else:
            ind = (mapco[:, 1] >= roi[0]) & (mapco[:, 1] < roi[1])
            return mapco[ind]



# Functions
def load_patients():
    '''Load patients from general table'''
    patients = pd.read_excel(table_filename, 'Patients', index_col=1)
    patients.index = pd.Index(map(str, patients.index))
    return patients


def load_patient(pname):
    '''Get the patient from the sequences ones'''
    patients = load_patients()
    patient = Patient(patients.loc[pname])
    return patient


def filter_patients_n_times(patients, n_times=3):
    '''Find what patients have at least n_times time points sequenced'''
    ind = np.zeros(len(patients), bool)
    for i, (pname, patient) in enumerate(patients.iterrows()):
        patient = Patient(patient)
        patient.discard_nonsequenced_samples()
        if len(patient.times) >= n_times:
            ind[i] = True

    return ind


def convert_date_deltas_to_float(deltas, unit='day'):
    '''Convert pandas date deltas into float'''
    nanoseconds_per_unit = {'day': 3600e9 * 24,
                            'month': 3600e9 * 24 * 365.25 / 12,
                            'year': 3600e9 * 24 * 365.25,
                           }
    return np.array(deltas, float) / nanoseconds_per_unit[unit]
