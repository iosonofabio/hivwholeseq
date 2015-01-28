# vim: fdm=marker
'''
author:     Richard Neher
date:       13/01/15
content:    Jsonification and annotation of trees.
'''
# Modules
import os

from hivwholeseq.utils.tree import tree_to_json
from hivwholeseq.utils.generic import read_json, write_json



# Functions
def annotate_from_name(T):
    '''
    loop over all nodes and extract information encoded in their names
    mostly useful for terminal nodes, attempted for internal nodes as well
    '''
    for node in T.get_terminals():
        anno = parse_leaf_name(node.name)
        for k, val in anno.iteritems():
            node.__setattr__(k, val)

    for node in T.get_nonterminals():
        anno = parse_internal_name(node.name)
        for k, val in anno.iteritems():
            node.__setattr__(k, val)


def parse_leaf_name(name):
    anno = {}
    entries = name.split('_')
    anno['name'] = entries[0]

    # go over fields defined globally above and use the next entry in the name
    # as value for the field. try int, float, otherwise string
    for field in fields:
        if field in entries:
            try:
                anno[field] = int(entries[entries.index(field)+1])
            except ValueError:
                try:
                    anno[field] = float(entries[entries.index(field)+1])
                except ValueError:
                    anno[field] = entries[entries.index(field)+1]
            except IndexError:
                anno[field] = "undefined"

    # add fields that require special parsing
    try:
        anno['DSI'] = int(entries[1].rstrip('days'))
    except:
        anno['DSI'] = 'nan' 

    return anno


def parse_internal_name(name):
    if isinstance(name, basestring):
        return parse_leaf_name(name)
    else:
        return {}


def annotate_from_dict(T, anno):
    for node in T.get_nonterminals()+T.get_terminals():
        if node.name in anno:
            for k, val in anno.iteritems():
                node.__setattr__(k, val)


def main(params):
    from Bio import Phylo
    if isinstance(params.tree, basestring):
        if os.path.isfile(params.tree):
            T = Phylo.read(params.tree,'newick')
            
            annotate_from_name(T.root)
            if hasattr(params,"anno") and params.anno is not None:
                if os.path.isfile(params.anno):
                    import pickle
                    with open(params.anno,'r') as infile:
                        anno = pickle.load(infile)
                        annotate_with_dict(T.root, anno)

            tree_json = tree_to_json(T.root)
            return tree_json
    print "No good tree file found", params
    return None



# Script
if __name__=="__main__":
    import argparse,sys
    parser = argparse.ArgumentParser(description='Add annotation to existing trees')
    parser.add_argument('--tree', required=True,help='tree to annotate')
    parser.add_argument('--anno',help='extra annotations (expecting a pickled dict)')
    params = parser.parse_args()

    tree_json = main(params)
    write_json(tree_json, '.'.join(params.tree.split('.')[:-1])+'.json')

