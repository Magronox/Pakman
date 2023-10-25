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
output = run_walk(G, [])

number_of_compactions = 20
name = "ecoli_illumina_10x_part_.fasta"
output, G = read_and_walk(k,number_of_compactions,name)

string_ = randstring("ACGT",1000)
output,G = read_and_walk_string(32,20,string_)
```
