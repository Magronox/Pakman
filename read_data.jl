using FastaIO

FastaReader("ecoli_illumina_10x_part_.fasta") do fr
    for (desc, seq) in fr
        DNAS
        break
    end
end

file = "ecoli_illumina_10x_part_.fasta"

function read_and_walk(k :: Int64, number_of_compactions :: Int64, name :: String)
    kmer_list = []
    for (~,seq) in FastaReader(name)
        seq = string_to_DNASeq(input)
        push!(kmer_list, read_kmer(DNA_seq, length(input),k))
        
    end
    coverage = 5
    contig_list = []
    G = graph_creator(kmer_list,['A','C','G','T'], coverage)
    G_new,ls = compact_graph(G,k,Number)
    run_walk(G,ls)
end


## Test:
## k = 3
## number_of_compactions = 4
## name = "ecoli_illumina_10x_part_.fasta"
## read_and_walk(k,number_of_compactions,name)