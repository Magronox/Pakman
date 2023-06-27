# Pakman
Pakman: Parallel Genome Assembly Algorithm

graph_walk.jl is the final file. If you run that file then you can use all the functions. Here's a code for a simple test

```
slen = 1000
input = randstring("ACGT",slen)

DNA_seq = string_to_DNASeq(input)
k=5
kmer_list = read_kmer(DNA_seq, length(input),k)

G = graph_creator(kmer_list,['A','C','G','T'], 5)
I = IS(G,['A','C','G','T'],k)
G_new, ~ = compact_graph(G,k,5)
output = run_walk(G)
```
