import json,os

fields = ['DSI','seq','muts','fmax', 'freq','readcount', 'VL', 'CD4']

def read_json(file_name):
    '''
    read the json file with name file_name into a dict.
    '''
    try:
        with open(file_name, 'r') as infile:
            data = json.load(infile)
    except IOError:
        print("Cannot open "+file_name)
    return data

def write_json(data, file_name, indent=None):
    '''
    dump a dict as a json to file
    params:
    data       -- dictionary
    file_name  -- name of file to save dict
    indent     -- indentation of json file, default None, 0 = only line breaks
    '''
    try:
        with open(file_name, 'w') as outfile:
            json.dump(data, outfile, indent=indent)
    except IOError:
        print("Cannot open "+file_name)

def tree_to_json(node):
    '''
    recursively go through the tree and create nested dictionaries
    each corresponding to one node. 
    '''
    json = {'name':node.name}
    json = {'branch_length':node.branch_length}
    for field in fields:
        if hasattr(node, field):
            val = node.__getattribute__(field)
            if val is None:
                json[field] = "undefined"
            else:
                json[field] = val

    # repeat for all children
    if len(node.clades):
        json["children"] = []
        for ch in node.clades:
            json["children"].append(tree_to_json(ch))

    return json

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

if __name__=="__main__":
    import argparse,sys
    parser = argparse.ArgumentParser(description='Add annotation to existing trees')
    parser.add_argument('--tree', required=True,help='tree to annotate')
    parser.add_argument('--anno',help='extra annotations (expecting a pickled dict)')
    params = parser.parse_args()

    tree_json = main(params)
    write_json(tree_json, '.'.join(params.tree.split('.')[:-1])+'.json')

