cou_corr = {}
for (seq, count) in cou.most_common():
    if not len(cou_corr):
        cou_corr[seq] = count
        continue

    cous = np.array([np.fromstring(key, 'S1') for key in cou_corr.keys()], 'S1', ndmin=2)
    seqm = np.fromstring(seq, 'S1')
    dm = (seqm != cous).sum(axis=1)
    imin = np.argmin(dm)
    if dm[imin] <= 6:
        cou_corr[''.join(cous[imin])] += count

    else:
        cou_corr[seq] = count

cou_corr_filt = Counter({seq: count for (seq, count) in cou_corr.iteritems()
                 if count > 1})

seqs_filt = [SeqRecord(Seq(seq, alphabet=ambiguous_dna), id=str(count),
                       name=str(count), description='')
             for (seq, count) in cou_corr_filt.most_common()]
