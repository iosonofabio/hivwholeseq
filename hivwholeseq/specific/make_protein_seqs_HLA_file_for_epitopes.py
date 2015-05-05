# vim: fdm=indent
'''
author:     Fabio Zanini
date:       05/05/15
content:    Make files with the protein sequences and HLA alleles for patient
            for the online epitope predictor:

                http://tools.immuneepitope.org/mhci/
'''
# Modules
import os
import sys

from hivwholeseq.patients.patients import load_patients, iterpatient


# Functions
def write_proteins_file(patient):
    from Bio.SeqRecord import SeqRecord
    from Bio import SeqIO

    regions = ['gag', 'pol', 'vif', 'vpu', 'vpr', 'gp120', 'gp41', 'nef']
    seqs = []
    for region in regions:
        prot = patient.get_reference(region).seq.translate()
        prot = prot[:str(prot).find('*')]
        seqs.append(SeqRecord(prot,
                              id=patient.code+'_'+region,
                              name=patient.code+'_'+region,
                              description=''))

    SeqIO.write(seqs, '/home/fabio/Desktop/seqs_'+patient.code+'.fasta', 'fasta')


def write_allele_file(patient):
    hlastr = '\n'.join(['HLA-'+x[:-3]+',9' for x in patient.get_hla_type()])
    with open('/home/fabio/Desktop/allele_file_'+patient.code+'.csv', 'w') as f:
        f.write(hlastr)



# Script
if __name__ == '__main__':

    patients = load_patients()
    for pname, patient in iterpatient(patients):
        print patient.code
        write_proteins_file(patient)
        write_allele_file(patient)
