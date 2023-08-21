using FastaIO
include("graph_walk.jl")

file = "ecoli_illumina_10x_part_.fasta"

function read_and_walk(k :: Int64, number_of_compactions :: Int64, name :: String)
    kmer_list =  DefaultDict{DNASeq, Int64}(0)
    
    FastaReader("ecoli_illumina_10x_part_.fasta") do fr
        
        for (~,seq) in fr
            input = string_to_DNASeq(seq)
            temp_kmers = read_kmer(input, input[1].len + (length(input)-1)*64,k)
            for(t,v) in temp_kmers
                kmer_list[t] += v
            end

            break
        end
    end
        
    
    coverage = 5
    contig_list = []
    G = graph_creator(kmer_list,['A','C','G','T'], coverage)

    print("Starting size ",length(G),"\n")
    G,ls = compact_graph!(G,number_of_compactions)
    run_walk(G,ls)
    contig_list, G
end


## Test:
## k = 3
## number_of_compactions = 4
## name = "ecoli_illumina_10x_part_.fasta"
## read_and_walk(k,number_of_compactions,name)