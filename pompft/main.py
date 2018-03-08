import os
import types
import re
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def _read_genome(f):
    '''
    Loads the genome file into memory.
    '''
    name, ext = os.path.splitext(f)
    # SeqIO doesn't recognize the .
    ext = ext.split('.')[-1]
    records = SeqIO.parse(open(f), ext)
    assert isinstance(records, types.GeneratorType)
    return records

def _kmers(seq):
    for	i in xrange(0, len(st)-(k-1)): yield st[i:i+k]

def _np_kmers(seq, start, stop, k=5):
    '''
    Returns an 2D np.array where row 1 is the start position for the kmer,
    and row 2 is the kmer.
    '''
    gen = _kmers(seq)
    return gen

def _panseq_positions(record_id):
    pattern = r'\(.[0-9]*\.\.[0-9]*\)'
    bracketed = re.search(pattern, record_id).group()
    # TODO: clean this up
    l = bracketed.split('..')
    start = l[0][1:]
    stop = l[-1][:-1]
    return start, stop

def _panseq_kmers(records):
    for record in records:
        print(record.id)
        seq = record.seq
        start, stop = _panseq_positions(record.id)
        np_a = _np_kmers(seq, start, stop)

def main():
    pass
