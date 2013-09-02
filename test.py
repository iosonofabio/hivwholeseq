#def pair_generator(iterable):
#    it = iter(iterable)
#    while True:
#        try:
#            a = it.next()
#            b = it.next()
#            yield (a, b)
#        except StopIteration:
#            raise
##
##a = [1, 2, 3, 4]
##b = [1, 2, 3, 4, 5]
##
##print [aa for aa in pair_generator(a)]
##print [aa for aa in pair_generator(b)]
#
#def identity(x):
#    '''Test function for import only'''
#    return x
#
#
#print __file__
#
#
#bs = []
#for aa in a:
#    b = np.zeros((4, 6, len(aa[0][0])), int)
#    for i, aaa in enumerate(aa):
#        for ii, aaaa in enumerate(aaa):
#            b[i, ii] = aaaa
#    bs.append(b)

#import mapping.sequence_utils
#from mapping.sequence_utils.annotate_HXB2 import load_HXB2

#import Bio.SeqIO as SeqIO
#from mapping.adapter_info import load_adapter_table, foldername_adapter
#def consensus_file(data_folder, adaID, fragment):
#    return data_folder+foldername_adapter(adaID)+'consensus_'+fragment+'.fasta'
#data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/'
#adaIDs = load_adapter_table(data_folder)['ID']
#fragments = ['F'+str(i) for i in xrange(1, 7)]
#for adaID in adaIDs:
#    for fragment in fragments:
#        fn = consensus_file(data_folder, adaID, fragment)
#        seq = SeqIO.read(fn, 'fasta')
#        print seq.name
#        #seq.name = seq.id = 'adaID_'+'{:02}'.format(adaID)+'_'+seq.id+'_consensus'
#        #with open(fn, 'w') as f:
#        #    SeqIO.write(seq, f, 'fasta')
#
#

