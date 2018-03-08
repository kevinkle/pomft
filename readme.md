# What is this?

A method for storing and updating prokaryotic pangenomes using a [Li-Stephens model](https://www.ncbi.nlm.nih.gov/pubmed/14704198), [effectively](https://www.ncbi.nlm.nih.gov/pubmed/27769991) a hidden markov model.
The storage is handled using [lemongraph](https://github.com/NationalSecurityAgency/lemongraph), a lmdb-backed graph store, and prototyped with [networkx](https://github.com/networkx/networkx).
We use a k-mer approach as described in [Graphtyper](https://github.com/DecodeGenetics/graphtyper), namely 5-mers with one overlapping character on either end.
However, instead of using a directed acyclic graph, we use a HMM as it benefits updating the graph - with new or existing edges have changing weights which ties into search performance (& possibly other uses).

# Installation

## Jupyter Problems
https://stackoverflow.com/questions/37891550/jupyter-notebook-running-kernel-in-different-env

# Future considerations

Right now, we focus on the storing and updating a pre-computed pangenome.
It might be possible to extend this method to create the pangenome in the beginning too, as nodes will naturally cluster by the number / weight (occurrence) of the edges when the graph is updated.

lemongraph also provides historical views of the graph, and there might be a way to extract meaningful data out of the differences as well.
