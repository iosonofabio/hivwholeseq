# vim: fdm=marker
'''
author:     Fabio Zanini
date:       17/06/14
content:    Check whether we need a new initial reference. Reasons might be a new
            sample sequenced that comes before all current ones, or the current
            reference has some genes not properly assembled.
'''
# Modules
import sys
import os
import numpy as np
import argparse
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import ambiguous_dna
from Bio import SeqIO

from hivwholeseq.samples import SampleSeq
from hivwholeseq.patients.patients import load_patients, load_patient, Patient
from hivwholeseq.patients.filenames import get_sample_foldername
from hivwholeseq.genome_info import locate_gene
from hivwholeseq.primer_info import fragments_genes


# Functions
def check_similarity_initial_sample(refseq, sample_seq, fragment, VERBOSE=0, maxdiff=10):
    '''Check whether the reference looks similar to the initial sample'''
    consseq = SeqIO.read(sample_seq.get_consensus_filename(fragment), 'fasta')

    from seqanpy import align_global
    (score, ali1, ali2) = align_global(str(refseq.seq), str(consseq.seq), band=50)
    alim = np.zeros((2, len(ali1)), 'S1')
    alim[0] = np.fromstring(ali1, 'S1')
    alim[1] = np.fromstring(ali2, 'S1')
    n_diff = (alim[0] != alim[1]).sum()
    if VERBOSE >= 2:
        print fragment+': difference between ref and initial consensus:', n_diff

    if n_diff > maxdiff:
        print 'ERROR: '+fragment+', reference is not similar to initial consensus ('+\
                str(sample_init_seq.name)+', '+\
                str(n_diff)+' differences)'
        return False
    elif VERBOSE >=3:
        print 'OK: reference is similar to initial consensus ('+\
                str(sample_init_seq.name)+', '+sample_seq['seq run']+' '+sample_seq.adapter+', '+\
                str(n_diff)+' differences)'

    return True


def check_has_complete_codons(gene, genename, VERBOSE=0):
    '''Check that the length is multiple of 3'''
    if len(gene) % 3:
        print 'ERROR: '+genename+' has a length not 3 * X'
        return False
    elif VERBOSE >= 3:
        print 'OK: '+genename+' has a length 3 * X'
    return True


def check_start_aminoacid(prot, genename, VERBOSE=0):
    '''Check whether the protein starts with an M'''
    # Pol starts with an F instead of an M
    from collections import defaultdict
    start_aa = defaultdict(lambda: 'M')
    start_aa['pol'] = 'F'
    if prot[0] != start_aa[genename]:
        print 'ERROR: '+genename+' does not start with an '+start_aa[genename]+'!'
        return False
    elif VERBOSE >= 3:
        print 'OK: '+genename+' starts with an '+start_aa[genename]
    return True


def check_has_end(prot_ref, genename, VERBOSE=0):
    '''Check whether it has a stop codon'''
    if prot_ref[-1] != '*':
        if VERBOSE >= 1:
            print 'ERROR: '+genename+' does not end!'
        return False
    elif VERBOSE >= 3:
        print 'OK: '+genename+' ends with a *'
    return True


def check_has_premature_stops(prot_ref, genename, VERBOSE=0):
    '''Check for premature stop codons'''
    protm = np.array(prot_ref)
    ind_star = (protm == '*').nonzero()[0]
    if (len(ind_star) != 1) or (ind_star[0] != len(protm) - 1):
        if VERBOSE >= 1:
            print 'ERROR: '+genename+' has premature stop codons!'
        return False
    elif VERBOSE >= 3:
        print 'OK: '+genename+' has no premature stop codons'
    return True


def check_has_similar_length(len_ref, len_HXB2, genename, VERBOSE=0, maxdiff=30):
    '''Does the gene have similar length like the HXB2?'''
    if len_ref < len_HXB2 - maxdiff:
        print 'ERROR: '+genename+' too short! (ref '+str(len_ref)+', HXB2 '+str(len_HXB2)+')'
        return False
    elif len_ref > len_HXB2 + maxdiff:
        print 'ERROR: '+genename+' too long! (ref '+str(len_ref)+', HXB2 '+str(len_HXB2)+')'
        return False
    elif VERBOSE >= 3:
        print 'OK: '+genename+' has the right length (ref '+str(len_ref)+', HXB2 '+str(len_HXB2)+')'
    return True


def check_has_premature_stops_noend(prot_ref, genename, VERBOSE=0):
    '''Check for premature stop codons'''
    protm = np.array(prot_ref)
    ind_star = (protm == '*').nonzero()[0]
    if len(ind_star):
        if VERBOSE >= 1:
            print 'ERROR: '+genename+' has premature stop codons!'
        return False
    elif VERBOSE >= 3:
        print 'OK: '+genename+' has no premature stop codons'
    return True


def get_fragment_length_HXB2(frag_spec):
    '''Get the length of a fragment in HXB2'''
    from hivwholeseq.primer_info import primers_coordinates_HXB2
    pr_coord_HXB2 = primers_coordinates_HXB2[frag_spec]
    len_HXB2 = pr_coord_HXB2[1][0] - pr_coord_HXB2[0][1]
    return len_HXB2


def get_gene_HXB2(genename):
    '''Get a gene or exon in HXB2'''
    from operator import attrgetter
    from hivwholeseq.reference import load_custom_reference
    HXB2 = load_custom_reference('HXB2', format='gb')
    if genename not in ('tat1', 'tat2', 'rev1', 'rev2'):
        gene_coord = HXB2.features[map(attrgetter('id'), HXB2.features).index(genename)]
        gene_HXB2 = gene_coord.extract(HXB2)
        return gene_HXB2

    else:
        exon_n = int(genename[-1])
        genename = genename[:-1]
        gene_coord = HXB2.features[map(attrgetter('id'), HXB2.features).index(genename)]
        exon_coord = gene_coord.location.parts[exon_n - 1]
        exon_HXB2 = exon_coord.extract(HXB2)
        return exon_HXB2


def check_length_fragment(refseq, frag_spec, VERBOSE=0, tolerance=50):
    '''Check the length of the fragment compared to HXB2'''
    fragment = frag_spec[:2]
    len_HXB2 = get_fragment_length_HXB2(frag_spec)
    if len(refseq) < len_HXB2 - tolerance:
        print 'ERROR: '+fragment+' too short! ('+str(len(refseq))+' vs '+str(len_HXB2)+' in HXB2)'
        return False
    elif len(refseq) > len_HXB2 + tolerance:
        print 'ERROR: '+fragment+' too long! ('+str(len(refseq))+' vs '+str(len_HXB2)+' in HXB2)'
        return False
    elif VERBOSE >= 3:
        print 'OK: '+fragment+' has approximately the right length ('+str(len(refseq))+' vs '+str(len_HXB2)+' in HXB2)'

    return True


def check_F1(refseq, spec, VERBOSE=0):
    '''Check fragment F1: gag, pol'''
    check = check_length_fragment(refseq, 'F1'+spec, VERBOSE=VERBOSE, tolerance=50)
    if not check:
        return False

    # Check gag (should be complete)
    genename = 'gag'
    (start, end, start_found, end_found) = locate_gene(refseq, genename, VERBOSE=VERBOSE)
    if (not start_found) or (not end_found):
        print 'ERROR: '+genename+' not found in F1!'
        return False
    elif VERBOSE >= 3:
        print 'OK: start and end of '+genename+' found'
    
    gene_HXB2 = get_gene_HXB2(genename)
    check = check_has_similar_length(end - start, len(gene_HXB2), genename, VERBOSE=VERBOSE, maxdiff=30)
    if not check:
        return False

    geneseq = refseq[start: end]
    gene = geneseq.seq
    check = check_has_complete_codons(gene, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    prot = gene.translate()
    check = check_start_aminoacid(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    check = check_has_end(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    check = check_has_premature_stops(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    # Check pol (there should be the start)
    genename = 'pol'
    (start, end, start_found, end_found) = locate_gene(refseq, genename, VERBOSE=VERBOSE)
    if (not start_found):
        print 'ERROR: start of '+genename+' not found in F1!'
        return False
    elif VERBOSE >= 3:
        print 'OK: start of '+genename+' found'

    geneseq = refseq[start:]
    geneseq = geneseq[:len(geneseq) - len(geneseq) % 3]
    gene = geneseq.seq
    prot = gene.translate()
    check = check_start_aminoacid(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    check = check_has_premature_stops_noend(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    return True


def check_F2(refseq, spec, VERBOSE=0):
    '''Check fragment F2: gag, pol'''
    check = check_length_fragment(refseq, 'F2'+spec, VERBOSE=VERBOSE, tolerance=50)
    if not check:
        return False

    # Check gag (there should be end)
    genename = 'gag'
    (start, end, start_found, end_found) = locate_gene(refseq, genename, VERBOSE=VERBOSE)
    if (not end_found):
        print 'ERROR: end of '+genename+' not found in F2!'
        return False
    elif VERBOSE >= 3:
        print 'OK: end of '+genename+' found'

    geneseq = refseq[:end]
    geneseq = geneseq[len(geneseq) % 3:]
    gene = geneseq.seq
    prot = gene.translate()
    check = check_has_end(prot, 'gag', VERBOSE=VERBOSE)
    if not check:
        return False

    check = check_has_premature_stops(prot, 'gag', VERBOSE=VERBOSE)
    if not check:
        return False

    # Check pol (there should be the start)
    genename = 'pol'
    (start, end, start_found, end_found) = locate_gene(refseq, genename, VERBOSE=VERBOSE)
    if (not start_found):
        print 'ERROR: start of '+genename+' not found in F2!'
        return False
    elif VERBOSE >= 3:
        print 'OK: start of '+genename+' found'

    geneseq = refseq[start:]
    geneseq = geneseq[:len(geneseq) - len(geneseq) % 3]
    gene = geneseq.seq
    prot = gene.translate()
    check = check_start_aminoacid(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    check = check_has_premature_stops_noend(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    return True


def check_F3(refseq, spec, VERBOSE=0):
    '''Check fragment F3: end of pol'''
    check = check_length_fragment(refseq, 'F3'+spec, VERBOSE=VERBOSE, tolerance=50)
    if not check:
        return False

    # Check pol: this depends on the spec: for F3bo there should be the end,
    # anything else has only the middle (it's all pol!)
    genename = 'pol'
    if spec == 'bo':
        (start, end, start_found, end_found) = locate_gene(refseq, genename, VERBOSE=VERBOSE)
        if (not end_found):
            print 'ERROR: end of '+genename+' not found in F3!'
            return False
        elif VERBOSE >= 3:
            print 'OK: end of '+genename+' found'

        geneseq = refseq[:end]
        geneseq = geneseq[len(geneseq) % 3:]
        gene = geneseq.seq
        prot = gene.translate()
        check = check_has_end(prot, genename, VERBOSE=VERBOSE)
        if not check:
            return False

        check = check_has_premature_stops(prot, genename, VERBOSE=VERBOSE)
        if not check:
            return False

    else:
        # Try all 3 reading frames
        for offset in xrange(3):
            geneseq = refseq[offset:]
            geneseq = geneseq[: len(geneseq) - (len(geneseq) % 3)]
            gene = geneseq.seq
            prot = gene.translate()

            check = check_has_premature_stops_noend(prot, genename, VERBOSE=0)
            if check:
                if VERBOSE >= 3:
                    print 'OK: '+genename+' has no premature stop codons'
                break
        else:
            if VERBOSE >= 1:
                print 'ERROR: '+genename+' has premature stop codons in all reading frames!'
            return False

    return True


def check_F4(refseq, spec, VERBOSE=0):
    '''Check fragment F4: vif, vpr, vpu, tat1, rev1, env'''
    check = check_length_fragment(refseq, 'F4'+spec, VERBOSE=VERBOSE, tolerance=50)
    if not check:
        return False

    # Check env (there should be the start)
    genename = 'env'
    (start, end, start_found, end_found) = locate_gene(refseq, genename, VERBOSE=VERBOSE)
    if (not start_found):
        print 'ERROR: start of '+genename+' not found in F4!'
        return False
    elif VERBOSE >= 3:
        print 'OK: start of '+genename+' found'

    geneseq = refseq[start:]
    geneseq = geneseq[:len(geneseq) - len(geneseq) % 3]
    gene = geneseq.seq
    prot = gene.translate()
    check = check_start_aminoacid(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    check = check_has_premature_stops_noend(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    # Check vif (should be complete)
    genename = 'vif'
    (start, end, start_found, end_found) = locate_gene(refseq, genename, VERBOSE=VERBOSE)
    if (not start_found) or (not end_found):
        print 'ERROR: '+genename+' not found in F4!'
        return False
    elif VERBOSE >= 3:
        print 'OK: start and end of '+genename+' found'
    
    gene_HXB2 = get_gene_HXB2(genename)
    check = check_has_similar_length(end - start, len(gene_HXB2), genename, VERBOSE=VERBOSE, maxdiff=15)
    if not check:
        return False

    geneseq = refseq[start: end]
    gene = geneseq.seq
    check = check_has_complete_codons(gene, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    prot = gene.translate()
    check = check_start_aminoacid(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    check = check_has_end(prot, genename, VERBOSE=0)
    if check:
        if VERBOSE >= 3:
            print 'OK: '+genename+' ends with a *'
    else:
        # Vif tends to be a bit longer than in HXB2
        for nc in xrange(1, 4):
            gene_ext = refseq[start: end + 3 * nc].seq
            prot_ext = gene_ext.translate()
            check = check_has_end(prot_ext, genename, VERBOSE=0)
            if check:
                gene = gene_ext
                prot = prot_ext
                if VERBOSE:
                    print 'WARNING: '+genename+' actually ends '+str(nc)+' codons downstream'
                break
        else:
            print 'ERROR: '+genename+' does not end, not even slightly downstream'
            return False

    check = check_has_premature_stops(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    # Check vpu (should be complete)
    genename = 'vpu'
    (start, end, start_found, end_found) = locate_gene(refseq, genename, VERBOSE=VERBOSE)
    if (not start_found) or (not end_found):
        print 'ERROR: '+genename+' not found in F4!'
        return False
    elif VERBOSE >= 3:
        print 'OK: start and end of '+genename+' found'
    
    gene_HXB2 = get_gene_HXB2(genename)
    check = check_has_similar_length(end - start, len(gene_HXB2), genename, VERBOSE=VERBOSE, maxdiff=15)
    if not check:
        return False

    geneseq = refseq[start: end]
    gene = geneseq.seq
    check = check_has_complete_codons(gene, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    prot = gene.translate()
    check = check_start_aminoacid(prot, genename, VERBOSE=VERBOSE)
    if not check:
        print 'ERROR IN VPU STARTING CODON, CONTINUING!'
        #return False

    check = check_has_end(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    check = check_has_premature_stops(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    # Check vpr (should be complete)
    genename = 'vpr'
    (start, end, start_found, end_found) = locate_gene(refseq, genename, VERBOSE=VERBOSE)
    if (not start_found) or (not end_found):
        print 'ERROR: '+genename+' not found in F4!'
        return False
    elif VERBOSE >= 3:
        print 'OK: start and end of '+genename+' found'
    
    gene_HXB2 = get_gene_HXB2(genename)
    check = check_has_similar_length(end - start, len(gene_HXB2), genename, VERBOSE=VERBOSE, maxdiff=15)
    if not check:
        return False

    geneseq = refseq[start: end]
    gene = geneseq.seq
    check = check_has_complete_codons(gene, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    prot = gene.translate()
    check = check_start_aminoacid(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    check = check_has_end(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    check = check_has_premature_stops(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    # Check tat1 (first exon of tat, should be complete)
    genename = 'tat1'
    (start, end, start_found, end_found) = locate_gene(refseq, genename, VERBOSE=VERBOSE)
    if (not start_found) or (not end_found):
        print 'ERROR: '+genename+' not found in F4!'
        return False
    elif VERBOSE >= 3:
        print 'OK: start and end of '+genename+' found'
    
    gene_HXB2 = get_gene_HXB2(genename)
    check = check_has_similar_length(end - start, len(gene_HXB2), genename, VERBOSE=VERBOSE, maxdiff=15)
    if not check:
        return False

    geneseq = refseq[start: end]
    geneseq = geneseq[:len(geneseq) - len(geneseq) % 3]
    gene = geneseq.seq
    prot = gene.translate()
    check = check_start_aminoacid(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    check = check_has_premature_stops_noend(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    # Check rev1 (first exon of rev, should be complete)
    genename = 'rev1'
    (start, end, start_found, end_found) = locate_gene(refseq, genename, VERBOSE=VERBOSE)
    if (not start_found) or (not end_found):
        print 'ERROR: '+genename+' not found in F4!'
        return False
    elif VERBOSE >= 3:
        print 'OK: start and end of '+genename+' found'
    
    gene_HXB2 = get_gene_HXB2(genename)
    check = check_has_similar_length(end - start, len(gene_HXB2), genename, VERBOSE=VERBOSE, maxdiff=15)
    if not check:
        return False

    geneseq = refseq[start: end]
    geneseq = geneseq[:len(geneseq) - len(geneseq) % 3]
    gene = geneseq.seq
    prot = gene.translate()
    check = check_start_aminoacid(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    check = check_has_premature_stops_noend(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    return True


def check_F5(refseq, spec, VERBOSE=0):
    '''Check fragment F5: env'''
    if spec == 'a+bo':
        spec_inner = 'bo'
    else:
        spec_inner = spec

    check = check_length_fragment(refseq, 'F5'+spec_inner, VERBOSE=VERBOSE, tolerance=70)
    if not check:
        return False

    # Check env (there should be the start)
    genename = 'env'
    (start, end, start_found, end_found) = locate_gene(refseq, genename, VERBOSE=VERBOSE)
    if (not start_found):
        print 'ERROR: start of '+genename+' not found in F5!'
        return False
    elif VERBOSE >= 3:
        print 'OK: start of '+genename+' found'

    geneseq = refseq[start:]
    geneseq = geneseq[:len(geneseq) - len(geneseq) % 3]
    gene = geneseq.seq
    prot = gene.translate()
    check = check_start_aminoacid(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    check = check_has_premature_stops_noend(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    return True


def check_F6(refseq, spec, VERBOSE=0):
    '''Check fragment F6: end of env, tat2, rev2'''
    check = check_length_fragment(refseq, 'F6'+spec, VERBOSE=VERBOSE, tolerance=50)
    if not check:
        return False

    # Check env (there should be end)
    genename = 'env'
    (start, end, start_found, end_found) = locate_gene(refseq, 'env', VERBOSE=VERBOSE)
    if (not end_found):
        print 'ERROR: end of '+genename+' not found in F6!'
        return False
    elif VERBOSE >= 3:
        print 'OK: end of '+genename+' found'

    geneseq = refseq[:end]
    geneseq = geneseq[len(geneseq) % 3:]
    gene = geneseq.seq
    prot = gene.translate()
    check = check_has_end(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    check = check_has_premature_stops(prot, genename, VERBOSE=VERBOSE)
    if not check:
        print prot
        return False

    # Check tat2 (second exon of tat, should be complete)
    genename = 'tat2'
    (start, end, start_found, end_found) = locate_gene(refseq, genename, VERBOSE=VERBOSE)
    if (not start_found) or (not end_found):
        print 'ERROR: '+genename+' not found in F6!'
        return False
    elif VERBOSE >= 3:
        print 'OK: start and end of '+genename+' found'
    
    gene_HXB2 = get_gene_HXB2(genename)
    check = check_has_similar_length(end - start, len(gene_HXB2), genename, VERBOSE=VERBOSE, maxdiff=15)
    if not check:
        return False

    geneseq = refseq[start: end]
    geneseq = geneseq[len(geneseq) % 3:]
    gene = geneseq.seq
    prot = gene.translate()
    check = check_has_end(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    check = check_has_premature_stops(prot, genename, VERBOSE=VERBOSE)
    if not check:
        print 'ERROR IN TAT2 PREMATURE STOPS, CONTINUING!'

    # Check rev2 (second exon of rev, should be complete)
    genename = 'rev2'
    (start, end, start_found, end_found) = locate_gene(refseq, genename, VERBOSE=VERBOSE)
    if (not start_found) or (not end_found):
        print 'ERROR: '+genename+' not found in F6!'
        return False
    elif VERBOSE >= 3:
        print 'OK: start and end of '+genename+' found'

    # NOTE: rev2 overlaps with env gp41 and can have insertions or deletions
    gene_HXB2 = get_gene_HXB2(genename)
    check = check_has_similar_length(end - start, len(gene_HXB2), genename, VERBOSE=VERBOSE, maxdiff=45)
    if not check:
        return False

    geneseq = refseq[start: end]
    geneseq = geneseq[len(geneseq) % 3:]
    gene = geneseq.seq
    prot = gene.translate()
    check = check_has_end(prot, genename, VERBOSE=VERBOSE)
    if not check:
        # rev2 can end a bit early
        end_new = prot.rfind('*')
        if end_new == -1:
            return False
        if len(prot) - 1 - end_new < 20:
            print 'REV2 ENDS '+str(len(prot) - end_new - 1)+' AMINO ACIDS TOO EARLY!'
            prot = prot[:end_new + 1]

    check = check_has_premature_stops(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    return True


def check_genomewide(refseq, VERBOSE=0):
    '''Check the integrity of all genes in the genomewide consensus'''
    # Check single-exon genes
    length_tolerance = {'gag': 30, 'pol': 30, 'env': 70, 'vpr': 15, 'vpu': 15}
    for genename, tol in length_tolerance.iteritems():
        (start, end, start_found, end_found) = locate_gene(refseq, genename, VERBOSE=VERBOSE)
        if (not start_found) or (not end_found):
            print 'ERROR: '+genename+' not found in genomewide!'
            return False
        elif VERBOSE >= 3:
            print 'OK: start and end of '+genename+' found'
        
        gene_HXB2 = get_gene_HXB2(genename)
        check = check_has_similar_length(end - start, len(gene_HXB2), genename, VERBOSE=VERBOSE, maxdiff=tol)
        if not check:
            return False

        geneseq = refseq[start: end]
        gene = geneseq.seq
        check = check_has_complete_codons(gene, genename, VERBOSE=VERBOSE)
        if not check:
            return False

        prot = gene.translate()
        check = check_start_aminoacid(prot, genename, VERBOSE=VERBOSE)
        if (not check):
            if genename != 'vpu':
                return False
            else:
                print 'ERROR IN VPU STARTING CODON, CONTINUING!'

        check = check_has_end(prot, genename, VERBOSE=VERBOSE)
        if not check:
            return False

        check = check_has_premature_stops(prot, genename, VERBOSE=VERBOSE)
        if not check:
            return False

    # Vif is special because it can be longer than in HXB2
    genename = 'vif'
    (start, end, start_found, end_found) = locate_gene(refseq, genename, VERBOSE=VERBOSE)
    if (not start_found) or (not end_found):
        print 'ERROR: '+genename+' not found in genomewide!'
        return False
    elif VERBOSE >= 3:
        print 'OK: start and end of '+genename+' found'
    
    gene_HXB2 = get_gene_HXB2(genename)
    check = check_has_similar_length(end - start, len(gene_HXB2), genename, VERBOSE=VERBOSE, maxdiff=15)
    if not check:
        return False

    geneseq = refseq[start: end]
    gene = geneseq.seq
    check = check_has_complete_codons(gene, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    prot = gene.translate()
    check = check_start_aminoacid(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    check = check_has_end(prot, genename, VERBOSE=0)
    if not check:
        # Vif tends to be a bit longer than in HXB2
        for nc in xrange(1, 4):
            gene_ext = refseq[start: end + 3 * nc].seq
            prot_ext = gene_ext.translate()
            check = check_has_end(prot_ext, genename, VERBOSE=0)
            if check:
                gene = gene_ext
                prot = prot_ext
                if VERBOSE:
                    print 'WARNING: '+genename+' actually ends '+str(nc)+' codons downstream'
                break
        else:
            print 'ERROR: '+genename+' does not end, not even slightly downstream'
            return False

    check = check_has_premature_stops(prot, genename, VERBOSE=VERBOSE)
    if not check:
        return False

    # Check 2-exon genes
    for genename_whole in ('tat', 'rev'):
        genename = genename_whole+'1'
        (start, end, start_found, end_found) = locate_gene(refseq, genename, VERBOSE=VERBOSE)
        if (not start_found) or (not end_found):
            print 'ERROR: '+genename+' not found in genomewide!'
            return False
        elif VERBOSE >= 3:
            print 'OK: start and end of '+genename+' found'
        
        gene_HXB2 = get_gene_HXB2(genename)
        check = check_has_similar_length(end - start, len(gene_HXB2), genename, VERBOSE=VERBOSE, maxdiff=15)
        if not check:
            return False

        geneseq = refseq[start: end]
        geneseq = geneseq[:len(geneseq) - len(geneseq) % 3]
        gene = geneseq.seq
        prot = gene.translate()
        check = check_start_aminoacid(prot, genename, VERBOSE=VERBOSE)
        if not check:
            return False

        start_exon1 = start
        end_exon1 = end

        genename = genename_whole+'2'
        (start, end, start_found, end_found) = locate_gene(refseq[end_exon1 + 2000:], genename, VERBOSE=VERBOSE)
        if (not start_found) or (not end_found):
            print 'ERROR: '+genename+' not found in genomewide!'
            return False
        elif VERBOSE >= 3:
            print 'OK: start and end of '+genename+' found'

        start += end_exon1 + 2000
        end += end_exon1 + 2000

        # NOTE: rev2 overlaps with env gp41 and can have insertions or deletions
        if genename == 'rev2':
            tol = 45
        else:
            tol = 15
        gene_HXB2 = get_gene_HXB2(genename)
        check = check_has_similar_length(end - start, len(gene_HXB2), genename, VERBOSE=VERBOSE, maxdiff=tol)
        if not check:
            return False

        geneseq = refseq[start: end]
        geneseq = geneseq[len(geneseq) % 3:]
        gene = geneseq.seq
        prot = gene.translate()
        check = check_has_end(prot, genename, VERBOSE=VERBOSE)
        if not check:
            return False

        check = check_has_premature_stops(prot, genename, VERBOSE=VERBOSE)
        if not check:
            return False

        start_exon2 = start
        end_exon2 = end

        genename = genename_whole
        gene_HXB2 = get_gene_HXB2(genename)

        from Bio.SeqFeature import FeatureLocation
        gene_loc = FeatureLocation(start_exon1, end_exon1, strand=+1) + \
                   FeatureLocation(start_exon2, end_exon2, strand=+1)
        geneseq = gene_loc.extract(refseq)
        gene = geneseq.seq

        check = check_has_complete_codons(gene, genename, VERBOSE=VERBOSE)
        if not check:
            return False

        prot = gene.translate()
        check = check_start_aminoacid(prot, genename, VERBOSE=VERBOSE)
        if not check:
            return False

        check = check_has_end(prot, genename, VERBOSE=VERBOSE)
        if not check:
            return False

        check = check_has_premature_stops(prot, genename, VERBOSE=VERBOSE)
        if not check:
            return False

    return True


def check_genes(refseq, frag_spec, VERBOSE=0):
    '''Check whether a gene is present and intact'''
    if frag_spec == 'genomewide':
        return check_genomewide(refseq, VERBOSE=VERBOSE)

    fragment = frag_spec[:2]
    spec = frag_spec[2:]

    if fragment == 'F1':
        return check_F1(refseq, spec, VERBOSE=VERBOSE)
    elif fragment == 'F2':
        return check_F2(refseq, spec, VERBOSE=VERBOSE)
    elif fragment == 'F3':
        return check_F3(refseq, spec, VERBOSE=VERBOSE)
    elif fragment == 'F4':
        return check_F4(refseq, spec, VERBOSE=VERBOSE)
    elif fragment == 'F5':
        return check_F5(refseq, spec, VERBOSE=VERBOSE)
    elif fragment == 'F6':
        return check_F6(refseq, spec, VERBOSE=VERBOSE)
    else:
        raise ValueError('Fragment '+fragment+' not implemented')



# Script
if __name__ == '__main__':


    # Parse input args
    parser = argparse.ArgumentParser(description='Check patient samples')
    parser.add_argument('--verbose', type=int, default=0,
                        help='Verbosity level [0-3]')
    parser.add_argument('--patients', nargs='+',
                        help='Patient to analyze')
    parser.add_argument('--fragments', nargs='+',
                        help='Fragment to map (e.g. F1 F6)')
    parser.add_argument('--force', action='store_true',
                        help='Ignore a single bad fragment and move to the next')

    args = parser.parse_args()
    VERBOSE = args.verbose
    pnames = args.patients
    use_force = args.force
    fragments = args.fragments

    if not fragments:
        fragments = ['F'+str(i) for i in xrange(1, 7)]
    if VERBOSE >= 3:
        print 'fragments', fragments

    patients = load_patients()
    if pnames is not None:
        patients = patients.loc[patients.index.isin(pnames)]

    for pname, patient in patients.iterrows():
        patient = Patient(patient)

        if VERBOSE:
            print 'Patient:', patient.name

        patient.discard_nonsequenced_samples()

        # Check whether patient has at least three time points, else ignore
        if patient.samples.shape[0] < 3:
            print 'WARNING: patient has less than three samples sequenced. Skipping.'
            continue
    
        sample_init = patient.initial_sample

        for fragment in fragments:
            print fragment

            if (patient.name in ('15241', '15319')) and (fragment in ('F4', 'F5', 'F6')):
                sample_init_seq = SampleSeq(sample_init['samples seq'].iloc[1])
            else:
                sample_init_seq = SampleSeq(sample_init['samples seq'].iloc[0])
            frag_specs = sample_init_seq.regions_complete

            if VERBOSE:
                print 'Initial sample:', sample_init_seq.name, sample_init_seq['seq run'], sample_init_seq.adapter


            # Check whether a reference exists at all
            ref_fn = patient.get_reference_filename(fragment)
            if not os.path.isfile(ref_fn):
                print 'ERROR: reference for fragment', fragment, 'not found!'
                continue
            elif VERBOSE >= 3:
                print 'OK: reference file found'
    
            refseq = SeqIO.read(ref_fn, 'fasta')

            # Check whether the consensus from the first sample is similar to
            # the reference. If not, it's not going to work
            check = check_similarity_initial_sample(refseq, sample_init_seq, fragment,
                                                    VERBOSE=VERBOSE)
            if not check:
                if not use_force:
                    sys.exit()

            # Check whether genes are fine
            frag_spec = [fr for fr in frag_specs if fragment in fr][0]
            check = check_genes(refseq, frag_spec, VERBOSE=VERBOSE)
            if not check:
                print 'ERROR in', fragment
                if use_force:
                    continue
                else:
                    sys.exit()

        # Check genomewide if present
        ref_fn = patient.get_reference_filename('genomewide')
        if not os.path.isfile(ref_fn):
            if VERBOSE >= 2:
                print 'WARNING: genomewide reference not found'
            continue

        refseq = SeqIO.read(ref_fn, 'fasta')
        #n_diff = check_similarity_initial_sample(refseq, sample_init_seq, 'genomewide',
        #                                        VERBOSE=VERBOSE)
        #if n_diff > 10:
        #    print 'ERROR: genomewide reference is not similar to initial consensus ('+\
        #            str(n_diff)+' differences)'
        #    continue
        
        check = check_genes(refseq, 'genomewide', VERBOSE=VERBOSE)
        if not check:
            print 'ERROR in genomewide'
            sys.exit()

