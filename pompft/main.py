import os
import types
import re
import numpy as np
import networkx as nx
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def _read_genome(f):
    '''
    Loads the genome file into memory.
    '''
    name, ext = os.path.splitext(f)
    # SeqIO doesn't recognize the .
    ext = ext.split('.')[-1]
    # TODO: stop pretending it is a fasta.
    records = SeqIO.parse(open(f), 'fasta')
    assert isinstance(records, types.GeneratorType)
    return records

def _kmers(seq, k=5):
    '''
    Returns a generator for k-mers with 1 character overlapping on either end.
    '''
    for	i in range(0, len(seq)-(k-1), k-1): yield seq[i:i+k]

def _seed_graph(seq, G):
    p = 0
    q = 1
    while p < len(seq)-1:
        if seq[p] != '-':
            # Continue.
            while seq[q] != '-' or q == len(seq)-1:
                q += 1
            # Found unknown, store the seq connecting them in an edge.
            # Using int of the start/stop as the nodes.
            G.add_edge(p, q-1, seq=seq[p:q-1])
        else:
            # Advance q, and then create a '-' edge.
            while seq[q] == '-' or q == len(seq)-1:
                q += 1
            G.add_edge(p, q-1, seq='-')
        # Advance p. This will trip else on the next loop.
        p = q
        q += 1

def graph(records):
    G = nx.MultiGraph()

def _np_kmers(seq, k=5):
    '''
    Returns an 2D np.array where row 1 is the start position for the kmer,
    and row 2 is the kmer.
    '''
    generator = _kmers(seq, k)
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
