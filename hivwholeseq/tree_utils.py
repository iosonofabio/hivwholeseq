# vim: fdm=marker
'''
author:     Fabio Zanini
date:       11/09/14
content:    Support module with tree utility functions.
'''
# Functions
def build_tree_fasttree(filename, VERBOSE=0, rootname=None):
    '''Build phylogenetic tree using FastTree'''
    import subprocess as sp

    if VERBOSE >= 3:
        output = sp.check_output(['fasttree', '-nt', filename])
    else:
        output = sp.check_output(['fasttree', '-nt', filename], stderr=sp.STDOUT)
    tree_string = output.split('\n')[-2]

    import StringIO
    from Bio import Phylo
    tree = Phylo.read(StringIO.StringIO(tree_string), 'newick')

    if rootname is not None:
        for leaf in tree.get_terminals():
            if leaf.name == rootname:
                root = leaf
                break
        else:
            raise ValueError('Initial reference not found in tree')

        tree.root_with_outgroup(leaf)

    return tree

