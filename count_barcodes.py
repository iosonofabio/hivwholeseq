# vim: fdm=indent
'''
author:     Richard Neher, Fabio Zanini
date:       02/08/13
content:    Count barcode abundance. We expect a high abundance for the true ones,
            plus sequencing errors around them.
'''
# Modules
from collections import defaultdict
from Bio import SeqIO
import matplotlib.pylab as plt



# Globals
data_folder = '/ebio/ag-neher/share/data/MiSeq_HIV_Karolinska/run28_test_samples/'
datafile = data_folder+'lane1_NoIndex_L001_R2_001.fastq'



# Script
bc_counts = defaultdict(int)
rc = 0
max_reads=-1
with open(datafile, 'r') as infile:
    for read in SeqIO.parse(infile, 'fastq'):
        bc_counts[read.seq.tostring()] += 1
        rc+=1
        if rc==max_reads:
            break

print sorted(bc_counts.items(), key=lambda x:x[1], reverse=True)[:20]

plt.figure()
ax=plt.subplot(111)
plt.plot(range(1,len(bc_counts)+1), sorted(bc_counts.values(), reverse=True))
ax.set_yscale('log')
ax.set_xscale('log')
plt.xlabel('barcode rank')
plt.ylabel('abundance')


