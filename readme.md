# What is this?

Currently a half-empty repository.

# What is this suppose to be?

A method for storing and updating prokaryotic pangenomes using a quasi-[Li-Stephens model](https://www.ncbi.nlm.nih.gov/pubmed/14704198), [effectively](https://www.ncbi.nlm.nih.gov/pubmed/27769991) a hidden markov model.
The storage is handled using [lemongraph](https://github.com/NationalSecurityAgency/lemongraph), a lmdb-backed graph store, and prototyped with [networkx](https://github.com/networkx/networkx).
We incorporate elements from [graphtyper](https://github.com/DecodeGenetics/graphtyper)'s directed acyclic graph, along with emission probabilities from Li-Stephens.
(graphtyper also has an interesting k-mer approach for indexing, namely 5-mers with one overlapping character on either end.)
When updating the graph with different aligned haplotypes, our HMM can accommodate how common specific segments are in the edge weights which isn't possible with a DAG.
When sequence reads are realigned, a HMM may also be better for unindexed search (& possibly other uses), and would benefit from more data being fed through it.

Uses an aligned core genome, either `core_gene_alignment.aln` from Roary or `snp.phylip` from Panseq (I think???).

# Installation

## Jupyter Problems
https://stackoverflow.com/questions/37891550/jupyter-notebook-running-kernel-in-different-env

# Future considerations

Right now, we focus on the storing and updating a pre-computed pangenome.
It might be possible to extend this method to create the pangenome in the beginning too, as nodes will naturally cluster by the number / weight (occurrence) of the edges when the graph is updated.

lemongraph also provides historical views of the graph, and there might be a way to extract meaningful data out of the differences as well.

# Sample data

The .fasta files was annotated using Prokka and the pangenome generated with Roary following https://github.com/microgenomics/tutorials/blob/master/pangenome.md (except with Escherichia).
