# vim: fdm=marker
'''
author:     Fabio Zanini
date:       08/11/14
content:    Store annotated initial references and other sequences for the
            website.
'''
# Modules
import os
import sys
from Bio import SeqIO
from datetime import date

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.patients import load_patients, Patient
from hivwholeseq.website.filenames import (get_patient_reference_filename,
                                           get_coordinate_map_filename)
from hivwholeseq.reference import load_custom_reference



# Functions
def delete_features(ref):
    '''Delete useless features'''
    fea_useless = ['gp120_noVloops', 'chunk']
    ref.features = [fea for fea in ref.features
                    if (fea.id not in fea_useless) and (fea.type not in fea_useless)]


def add_annotations_fragments(ref):
    '''Add annotations of fragments to reference'''
    from Bio.SeqFeature import SeqFeature, FeatureLocation

    if 'HXB2' not in ref.name:
        raise NotImplementedError

    from hivwholeseq.data.primers import (
        primers_coordinates_HXB2_outer as primers)
    for fr in ['F1', 'F2', 'F3B', 'F4', 'F5a', 'F6']:
        coo = primers[fr]   
        start = coo[0][0]
        end = coo[1][1]
        fea = SeqFeature(FeatureLocation(start, end, strand=+1),
                         id=fr[:2],
                         type="fragment")
        ref.features.append(fea)


def save_reference(ref, filename):
    '''Save reference to file for website'''
    from hivwholeseq.utils.sequence import correct_genbank_features_save

    fmt = filename.split('.')[-1]

    # NOTE: the genbank Biopython writer is buggy, or the Genbank format is
    # buggy, anyway we need to trick feature types into notes.
    if fmt == 'gb':
        correct_genbank_features_save(ref)

    SeqIO.write(ref, filename, fmt)


def save_mapco(mapco, filename, header):
    '''Save coordinate map to file for website'''
    import numpy as np

    dirname = os.path.dirname(filename)
    mkdirs(dirname)
    np.savetxt(filename, mapco, fmt='%d',
               delimiter='\t',
               header=header)
    



# Script
if __name__ == '__main__':

    refname = 'HXB2'

    patients = load_patients()
    for pname, patient in patients.iterrows():
        patient = Patient(patient)
        print patient.code, patient.name

        # Get reference
        ref = patient.get_reference('genomewide', format='gb')

        # Mask patient name into code
        ref.id = 'reference_'+patient.code
        ref.name = 'reference_'+patient.code
        ref.description = 'reference for patient '+patient.code+', genomewide'

        ref.annotations['organism'] = 'HIV-1'
        ref.annotations['date'] = date.today().strftime('%d-%b-%Y').upper()
        ref.annotations['accessions'] = [patient.code+'_genomewide']

        # Delete useless annotations
        delete_features(ref)

        # Save to file GENBANK
        fn_out = get_patient_reference_filename(patient.code, format='gb')
        save_reference(ref, fn_out)

        # Save to file FASTA
        fn_out = get_patient_reference_filename(patient.code, format='fasta')
        save_reference(ref, fn_out)

        # Get coordinate maps
        mapco = patient.get_map_coordinates_reference('genomewide', refname=refname)
        fn_out = get_coordinate_map_filename(patient.code, refname=refname)
        header = refname+'\t'+patient.code
        save_mapco(mapco, fn_out, header)

    # Store HXB2
    print refname
    ref = load_custom_reference(refname, format='gb')

    # Delete useless annotations
    delete_features(ref)

    # Add fragment annotations
    add_annotations_fragments(ref)

    # Save to file GENBANK
    fn_out = get_patient_reference_filename(refname, format='gb')
    save_reference(ref, fn_out)

    # Save to file FASTA
    fn_out = get_patient_reference_filename(refname, format='fasta')
    save_reference(ref, fn_out)
