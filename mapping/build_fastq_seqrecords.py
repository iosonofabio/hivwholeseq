def build_fastq_records(reads):
    '''Build FASTQ SeqRecords out of BAM reads
    
    Note: chunks copied from Bio.SeqIO.QualityIO.py
    '''

    # Precompute conversion table
    SANGER_SCORE_OFFSET = ord("!")
    q_mapping = dict()
    for letter in range(0, 255):
        q_mapping[chr(letter)] = letter - SANGER_SCORE_OFFSET
    
    seqs = []
    for read in reads:
        # Get the sequence first
        descr = read.qname
        id = read.qname
        name = id
        from Bio.Alphabet import IUPAC
        record = SeqRecord(Seq(read.seq, IUPAC.ambiguous_dna),
                           id=id, name=name, description=descr)
    
        # Get the qualities second
        qualities = [q_mapping[letter] for letter in read.qual]
        if qualities and (min(qualities) < 0 or max(qualities) > 93):
            raise ValueError("Invalid character in quality string")

        #For speed, will now use a dirty trick to speed up assigning the
        #qualities. We do this to bypass the length check imposed by the
        #per-letter-annotations restricted dict (as this has already been
        #checked by FastqGeneralIterator). This is equivalent to:
        #record.letter_annotations["phred_quality"] = qualities
        dict.__setitem__(record._per_letter_annotations,
                         "phred_quality", qualities)

        seqs.append(record)

    return seqs
