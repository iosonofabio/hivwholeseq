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
    import numpy as np

    if isinstance(filename_or_ali, basestring):
        filename = filename_or_ali
    else:
        ali = filename_or_ali
        tmp_folder = os.getenv('HOME')+'/tmp/'
        filename = tmp_folder+'tmp_fasttree_'+str(np.random.randint(1000000000))+'.fasta'
        from Bio import AlignIO
        with open(filename, 'w') as outfile:
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


def find_parent(tree, node):
    '''Find the parent node of a tree node'''
    return tree.root.get_path(node)[-2]


def get_path_toroot(tree, node):
    '''Get the path to root'''
    return tree.root.get_path(node)[::-1]


def tree_to_json(node,
                 fields=('DSI', 'seq', 'muts',
                         'fmax', 'freq',
                         'readcount',
                         'VL', 'CD4',
                         'confidence'),
                ):
    '''Convert tree in nested dictionary (JSON)'''
    from numpy import isnan

    if fields is None:
        fieldsnode = set(node.__dict__.keys())
        fieldsnode -= set(['clades', 'width', '_color'])

    else:
        fieldsnode = ['name', 'branch_length'] + list(fields)

    json = {}
    for field in fieldsnode:
        if not hasattr(node, field):
            json[field] = "undefined"
            continue

        val = getattr(node, field)
        if val is None:
            # The root is missing a branch length (maybe a FastTree or Biopython bug)
            if field == 'branch_length': 
                json[field] = 0
                continue

            json[field] = "undefined"
            continue

        if isinstance(val, basestring):
            if val.lower() == 'nan':
                json[field] = "undefined"
                continue

            json[field] = val
            continue

        if isnan(val):
            json[field] = "undefined"
            continue

        json[field] = val

    # repeat for all children
    if len(node.clades):
        json["children"] = []
        for ch in node.clades:
            json["children"].append(tree_to_json(ch, fields=fields))

    return json


def tree_from_json(json):
    '''Convert JSON into a Biopython tree'''
    from Bio import Phylo

    def node_from_json(json, node):
        '''Biopython Clade from json (for recursive call)'''
        for attr,val in json.iteritems():
            if attr == 'children':
                for sub_json in val:
                    child = Phylo.BaseTree.Clade()
                    node.clades.append(child)
                    node_from_json(sub_json, child)
            else:
                try:
                    node.__setattr__(attr, float(val))
                except:
                    node.__setattr__(attr, val)

    tree = Phylo.BaseTree.Tree()
    node_from_json(json, tree.root)
    tree.root.branch_length=0.01
    return tree


def add_mutations_tree(tree, translate=False,
                       sequence_attrname='sequence',
                       mutation_attrname='muts'):
    '''Add mutations to a tree
    
    NOTE: the nodes must have sequences already
    '''
    import numpy as np
    from Bio.Seq import translate as tran
    from hivwholeseq.utils.tree import find_parent

    # NOTE: the root has no mutations by definition
    for node in tree.get_terminals() + tree.get_nonterminals():
        try:
            parent = find_parent(tree, node)
        except IndexError:
            continue

        pseq = getattr(parent, sequence_attrname)
        nseq = getattr(node, sequence_attrname)

        if translate:
            pseq = tran(pseq)
            nseq = tran(nseq)

        pseqm = np.fromstring(pseq, 'S1')
        nseqm = np.fromstring(nseq, 'S1')
        posm = (pseqm != nseqm).nonzero()[0]
        muts = [pseqm[p]+str(p+1)+nseqm[p] for p in posm]

        setattr(node, mutation_attrname, ', '.join(muts))
