from Bio import AlignIO
import numpy as np
import matplotlib.pyplot as plt

fname = 'alignment_all_with_NL4-3_aligned.fasta'
all_refs = ['B.FR.'] #, 'B.US.', 'gi|19']
ref = all_refs[0]

smoothing=200

with open(fname, 'r') as infile:
    aln = AlignIO.read(infile, 'fasta')

seqs = np.zeros((len(aln)-len(all_refs), aln.get_alignment_length()), dtype='|S1')
seq_ids = []
si=0
for seq in aln:
    if seq.id[:len(ref)]==ref:
        print 'ref:', seq.id
        ref_seq = np.array(seq.seq)
    elif seq.id[:len(ref)] not in all_refs:
        print 'query:', seq.id
        seqs[si,:] = np.array(seq.seq)
        seq_ids.append(seq.id)
        si+=1

plt.figure()
plt.title('total difference')
for si,seq in enumerate(seqs):
    plt.plot(np.convolve(np.ones(smoothing, float)/smoothing, seq!=ref_seq), label = seq_ids[si])

plt.xlabel('alignment coordinate')
plt.ylabel('smoothed difference')
plt.legend(loc=2)

plt.figure()
plt.title('snp difference')
for si,seq in enumerate(seqs):
    plt.plot(np.convolve(np.ones(smoothing, float)/smoothing, (seq!=ref_seq)*(seq!='-')*(ref_seq!='-')), label = seq_ids[si])

plt.xlabel('alignment coordinate')
plt.ylabel('smoothed difference')
plt.legend(loc=2)

plt.figure()
plt.title('indels difference')
for si,seq in enumerate(seqs):
    plt.plot(np.convolve(np.ones(smoothing, float)/smoothing, (seq!=ref_seq)*((seq=='-')+(ref_seq=='-')>0)), label = seq_ids[si])

plt.xlabel('alignment coordinate')
plt.ylabel('smoothed difference')
plt.legend(loc=2)

