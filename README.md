# Pakman
Serial implementation of Pakman: Parallel Genome Assembly Algorithm

read_data.jl can be used to run the algorithm on data. 

```
slen = 1000
input = randstring("ACGT",slen)

DNA_seq = string_to_DNASeq(input)
k=5
kmer_list = read_kmer(DNA_seq, length(input),k)

G = graph_creator(kmer_list,['A','C','G','T'], 5)
I = IS(G,['A','C','G','T'],k)
G_new, ~ = compact_graph(G,slen/2,5)
output = run_walk(G)
```
