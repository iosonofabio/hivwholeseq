# vim: fdm=marker
'''
author:     Fabio Zanini
date:       06/11/14
content:    Store phylogenetic trees in a suitable format for the website.

            NOTE: the trees have been JSONed and anonymyzed already.
'''
# Modules
import os
import shutil
import sys
import argparse

from hivwholeseq.utils.generic import mkdirs
from hivwholeseq.patients.filenames import get_consensi_tree_filename as gfn_in
from hivwholeseq.website.filenames import get_consensi_tree_filename as gfn_out
from hivwholeseq.reference import load_custom_reference
from hivwholeseq.website import _regions



# Globals



# Script
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Align consensi and make trees',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)    
    parser.add_argument('--patients', nargs='+',
                        default=['p1', 'p2', 'p3', 'p5', 'p6', 'p8', 'p9', 'p10', 'p11'],
                        help='Patients to analyze')
    parser.add_argument('--regions', nargs='*', default=_regions,
                        help='Regions to analyze (e.g. V3 F6)')

    args = parser.parse_args()
    pcodes = args.patients
    regions = args.regions

    # GLOABL CONSENSI TREES
    print 'All patients'
    for region in regions:
        print region,

        reg_tmp = '_'.join(pcodes)+'_'+region
        fn_in = gfn_in('all', reg_tmp, format='json')
        if not os.path.isfile(fn_in):
            print 'SKIP'
            continue

        # Write output
        fn_out = gfn_out('all', region, format='json')
        shutil.copy(fn_in, fn_out)
        print 'OK'


    # PATIENT CONSENSI TREES
    print 'Single patients'
    for region in regions:
        print region
        for pname in pcodes:
            fn_in = gfn_in(pname, region, format='json')
            if not os.path.isfile(fn_in):
                print 'SKIP'
                continue
    
            # Write output
            fn_out = gfn_out(pname, region, format='json')
            if os.path.isfile(fn_out):
                os.remove(fn_out)
            shutil.copy(fn_in, fn_out)
            print 'OK'

