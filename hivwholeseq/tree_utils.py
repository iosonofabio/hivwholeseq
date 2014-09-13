# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/09/14
content:    Support module with tree utility functions.
'''
# Functions
def build_tree_fasttree(filename_or_ali, rootname=None, VERBOSE=0):
    '''Build phylogenetic tree using FastTree
    
    Args:
      filename_or_ali: filename of a FASTA multiple sequence alignment, or a
                       Biopython alignment itself
          
    '''
    import os
    import subprocess as sp
    import StringIO
    from Bio import Phylo

    if isinstance(filename_or_ali, basestring):
        filename = filename_or_ali
    else:
        ali = filename_or_ali
        tmp_folder = '/ebio/ag-neher/share/tmp/'
        filename = tmp_folder+'tmp_fasttree_'+str(np.random.randint(1000000000))+'.fasta.gz'
        from Bio import AlignIO
        with gzip.open(filename, 'w') as outfile:
            AlignIO.write(ali, outfile, 'fasta')

    try:
        if VERBOSE >= 3:
            output = sp.check_output(['fasttree', '-nt', filename])
        else:
            output = sp.check_output(['fasttree', '-nt', filename], stderr=sp.STDOUT)
        tree_string = output.split('\n')[-2]

        tree = Phylo.read(StringIO.StringIO(tree_string), 'newick')
        tree.root.branch_length = 0.001

        if rootname is not None:
            if VERBOSE >= 2:
                print 'Reroot'
            for leaf in tree.get_terminals():
                if leaf.name == rootname:
                    root = leaf
                    break
            else:
                raise ValueError('Initial reference not found in tree')

            tree.root_with_outgroup(leaf)

    finally:
        if filename_or_ali != filename:
            os.remove(filename)

    return tree
