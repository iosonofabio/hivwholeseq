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
    def n_templates(self):
        '''Get the time course of the number of templates to PCR, limiting depth'''
        from hivwholeseq.patients.get_template_number import get_template_number
        n = [get_template_number(dilstr) for dilstr in self.samples.dilutions]
        n = np.ma.masked_invalid(n)
        return n


    @property
    def n_templates_viral_load(self):
        '''Get the number of templates, estimated from the viral load'''
        n = self.viral_load.copy()
        # We take 400 ul of serum
        n *= 0.4
        # We typically have 6 reactions with that total volume (plus the F4 dilution
        # series, but each of those uses only 0.1x template which is very little)
        n /= 6.1
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


    def get_fragmented_roi(self, roi, VERBOSE=0, **kwargs):
        '''Get a region of interest in fragment coordinates'''
        from hivwholeseq.patients.get_roi import get_fragmented_roi
        if isinstance(roi, basestring):
            roi = (roi, 0, '+oo')
        return get_fragmented_roi(self, roi, VERBOSE=VERBOSE, **kwargs)


    def get_reference_filename(self, fragment, format='fasta'):
        '''Get filename of the reference for mapping'''
        from hivwholeseq.patients.filenames import get_initial_reference_filename
        return get_initial_reference_filename(self.name, fragment, format)


    def get_reference(self, region, format='fasta'):
        '''Get the reference for a genomic region'''
        from Bio import SeqIO
        if region == 'genomewide':
            fragment = region
        else:
            (fragment, start, end) = self.get_fragmented_roi((region, 0, '+oo'), include_genomewide=True)
        refseq = SeqIO.read(self.get_reference_filename(fragment, format=format), format)
        if format in ('gb', 'genbank'):
            from hivwholeseq.sequence_utils import correct_genbank_features_load
            correct_genbank_features_load(refseq)
        if region != 'genomewide':
            refseq = refseq[start: end]
        return refseq


    def get_consensi_alignment_filename(self, region, format='fasta'):
        '''Get the filename of the multiple sequence alignment of all consensi'''
        from hivwholeseq.patients.filenames import get_consensi_alignment_filename
        return get_consensi_alignment_filename(self.name, region, format=format)


    def get_consensi_alignment(self, region, format='fasta'):
        '''Get the multiple sequence alignment of all consensi'''
        from Bio import AlignIO
        return AlignIO.read(self.get_consensi_alignment_filename(region,
                                                                 format=format),
                            format)


    def get_consensi_tree_filename(self, region, format='newick'):
        '''Get the filename of the consensi tree of the patient'''
        from hivwholeseq.patients.filenames import get_consensi_tree_filename
        return get_consensi_tree_filename(self.name, region, format=format)


    def get_consensi_tree(self, region, format='newick'):
        '''Get consensi tree from the patient'''
        import os.path

        if format == 'json':
            fn = self.get_consensi_tree_filename(region, format='json')
            if os.path.isfile(fn):
                from ..generic_utils import read_json
                from ..tree_utils import tree_from_json
                return tree_from_json(read_json(fn))

        fn = self.get_consensi_tree_filename(region, format='newick')
        if os.path.isfile(fn):
            from Bio import Phylo
            return Phylo.read(fn, 'newick')


    def get_local_tree_filename(self, region, format='json'):
        '''Get the filename of the consensi tree of the patient'''
        from hivwholeseq.patients.filenames import get_local_tree_filename
        return get_local_tree_filename(self.name, region, format=format)


    def get_local_tree(self, region):
        '''Get consensi tree from the patient'''
        import os.path

        fn = self.get_local_tree_filename(region, format='json')
        if os.path.isfile(fn):
            from ..generic_utils import read_json
            from ..tree_utils import tree_from_json
            return tree_from_json(read_json(fn))


    @staticmethod
    def get_initial_consensus_noinsertions(aft, VERBOSE=0):
        '''Make initial consensus from allele frequencies, keep coordinates and masked
        
        Args:
          aft (np.ma.ndarray): 3d masked array with the allele frequency trajectories

        Returns:
          np.ndarray: initial consensus, augmented with later time points at masked
          positions, with Ns if never covered
        '''
        af0 = aft[0]
        # Fill the masked positions with N...
        cons_ind = af0.argmax(axis=0)
        cons_ind[af0[0].mask] = 5
    
        # ...then look in later time points
        if aft.shape[0] == 1:
            return cons_ind
        for af_later in aft[1:]:
            cons_ind_later = af_later.argmax(axis=0)
            cons_ind_later[af_later[0].mask] = 5
            ind_Ns = (cons_ind == 5) & (cons_ind_later != 5)
            cons_ind[ind_Ns] = cons_ind_later[ind_Ns]
        return cons_ind


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


    def get_coverage_trajectories(self, region, use_PCR1=1):
        '''Get coverage as a function of time'''
        (act, ind) = self.get_allele_count_trajectories(region, use_PCR1=use_PCR1)
        return (act.sum(axis=1), ind)


    def get_allele_frequency_trajectories(self, region, use_PCR1=1, cov_min=1,
                                          depth_min=None, **kwargs):
        '''Get the allele frequency trajectories from files
        
        Args:
          region (str): region to study, a fragment or a genomic feature (e.g. V3)
          cov_min (int): minimal coverage accepted, anything lower are masked.
          depth_min (float): minimal depth, both by sequencing and template numbers.
            Time points with less templates are excluded, and positions are masked.
            For convenience depth is defined > 1, e.g. 100 takes frequencies down
            to 1%.
          **kwargs: passed down to the get_allele_count_trajectories method.
        '''
        (act, ind) = self.get_allele_count_trajectories(region,
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


    def get_allele_count_trajectories(self, region, use_PCR1=2, safe=False, **kwargs):
        '''Get the allele count trajectories from files
        
        Args:
          region (str): region to study, a fragment or a genomic feature (e.g. V3)
          **kwargs: passed down to the function (VERBOSE, etc.).

        Note: the genomewide counts are currently saved to file.
        '''
        from operator import itemgetter
        from .one_site_statistics import get_allele_count_trajectories

        # Fall back on genomewide counts if no single fragment is enough
        (fragment, start, end) = self.get_fragmented_roi((region, 0, '+oo'),
                                                         include_genomewide=True)
        (sns, act) = get_allele_count_trajectories(self.name, self.samples.index,
                                                   fragment,
                                                   use_PCR1=use_PCR1, **kwargs)
        # Select genomic region
        act = act[:, :, start: end]

        # Select time points
        ind = np.array([i for i, (_, sample) in enumerate(self.samples.iterrows())
                        if sample.name in map(itemgetter(0), sns)], int)

        # If safe, take only samples tagged with 'OK'
        if safe:
            ind_safe = np.zeros(len(ind), bool)
            for ii, i in enumerate(ind):
                sample = self.samples.iloc[i]
                from .get_roi import get_fragments_covered
                frags = get_fragments_covered(self, (fragment, start, end))
                ind_safe[ii] = all(getattr(sample, fr).upper() == 'OK'
                                   for fr in frags)

            act = act[ind_safe]
            ind = ind[ind_safe]


        return (act, ind)


    def get_mapped_filtered_filename(self, samplename, fragment, PCR=1):
        '''Get filename(s) of mapped and filtered reads for a sample'''
        from hivwholeseq.patients.filenames import get_mapped_filtered_filename
        return get_mapped_filtered_filename(self.patient, samplename, fragment, PCR=PCR)


    def get_divergence(self, region, **kwargs):
        '''Get genetic divergence of a region
        
        Args:
          **kwargs: passed to the allele frequency trajectories.
        '''
        from hivwholeseq.patients.get_divergence_diversity import get_divergence
        aft, ind = self.get_allele_frequency_trajectories(region, **kwargs)
        return (get_divergence(aft), ind)


    def get_diversity(self, region, **kwargs):
        '''Get geneticdiversity of a region
        
        Args:
          **kwargs: passed to the allele frequency trajectories.
        '''
        from hivwholeseq.patients.get_divergence_diversity import get_diversity
        aft, ind = self.get_allele_frequency_trajectories(region, **kwargs)
        return (get_diversity(aft), ind)


    @property
    def transmission_date(self):
        '''The most likely time of transmission'''
        return self['last negative date'] + \
                (self['first positive date'] - self['last negative date']) / 2


    def get_map_coordinates_reference_filename(self, fragment, refname='HXB2'):
        '''Get the filename of the coordinate map to an external reference'''
        from hivwholeseq.patients.filenames import get_coordinate_map_filename
        return get_coordinate_map_filename(self.name, fragment, refname=refname)


    def get_map_coordinates_reference(self, roi, refname='HXB2'):
        '''Get the map of coordinate to some external reference
        
        Returns:
          mapco (2D int array): the first column are the positions in the reference,
            the second column the position in the patient initial reference. 
        '''
        if isinstance(roi, basestring):
            region = roi
        else:
            region = roi[0]

        fn = self.get_map_coordinates_reference_filename(region, refname=refname)
        mapco = np.loadtxt(fn, dtype=int)

        if isinstance(roi, basestring):
            return mapco
        else:
            start = roi[1]
            ind = (mapco[:, 1] >= start)

            if roi[2] != '+oo':
                end = roi[2]
                ind &= (mapco[:, 1] < end)

            # The patient coordinates are referred to the roi itself!
            mapco[:, 1] -= start
            return mapco[ind]


    def get_local_haplotype_trajectories(self, fragment, start, end, VERBOSE=0,
                                         **kwargs):
        '''Get trajectories of local haplotypes'''
        if fragment not in ['F'+str(i) for i in xrange(1, 7)]:
            (fragment, start, end) = self.get_fragmented_roi((fragment, start, end),
                                                             VERBOSE=VERBOSE)
        ind = []
        haplos = []
        for i, sample in enumerate(self.itersamples()):
            try:
                haplo = sample.get_local_haplotypes(fragment, start, end,
                                                    VERBOSE=VERBOSE,
                                                    **kwargs)
            except IOError:
                continue

            haplos.append(haplo)
            ind.append(i)

        return (haplos, ind)


    def get_local_haplotype_count_trajectories(self, fragment, start, end, VERBOSE=0,
                                               **kwargs):
        '''Get trajectories of local haplotypes counts'''
        (haplos, ind) = self.get_local_haplotype_trajectories(fragment,
                                                              start, end,
                                                              VERBOSE=VERBOSE,
                                                              **kwargs)
        # Make trajectories of counts
        seqs_set = set()
        for haplo in haplos:
            seqs_set |= set(haplo.keys())
        seqs_set = list(seqs_set)
        hct = np.zeros((len(seqs_set), len(haplos)), int)
        for i, haplo in enumerate(haplos):
            for seq, count in haplo.iteritems():
                hct[seqs_set.index(seq), i] = count

        seqs_set = np.array(seqs_set, 'S'+str(np.max(map(len, seqs_set))))
        return (hct.T, ind, seqs_set)

    
    def get_region_count_trajectories(self, region, VERBOSE=0, **kwargs):
        '''Get trajectories of local haplotypes in a precompiled region'''
        from .get_tree_local import get_region_count_trajectories
        return get_region_count_trajectories(self, region, VERBOSE=VERBOSE)



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
