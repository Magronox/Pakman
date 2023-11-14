using FastaIO
include("graph_walk.jl")

file = "../data/random_output_500_100_10.fasta"

function read_and_walk(k :: Int64, number_of_compactions :: Int64, name :: String)
    kmer_list =  DefaultDict{DNASeq, Int64}(0)
    
    #FastaReader(name) do fr
    #    
    for (~,seq) in FastaReader(name)
        input = string_to_DNASeq(seq)
        temp_kmers = read_kmer(input, input[1].len + (length(input)-1)*64,k)
        for (t,v) in temp_kmers
            kmer_list[t] += v
        end
        #break
    end
    #end
        
    
    coverage = 10
    global contig_list = []
    G = graph_creator(kmer_list,['A','C','G','T'], coverage)

    print("Starting size ",length(G),"\n")
    G,ls = compact_graph!(G,number_of_compactions)
    
    #print("length()",length(ls))
    
    #append!(contig_list,ls)
    
    for i in ls
        if !(i in contig_list)
            append!(contig_list,i)
        end
    end
    output = run_walk(G,ls)
    for i in output
        if !(i in contig_list)
            append!(contig_list,i)
        end
    end
    deleteat!(contig_list, findall(x->x==Any[],contig_list))


    #print("\n new contig list\n",contig_list)
    unique(contig_list), G
end


function rw(k :: Int64, number_of_compactions :: Int64, name :: String)

    kmer_list =  DefaultDict{DNASeq, Int64}(0)
    
    #FastaReader(name) do fr
    #    
    for (~,seq) in FastaReader(name)
        input = string_to_DNASeq(seq)
        temp_kmers = read_kmer(input, input[1].len + (length(input)-1)*64,k)
        for (t,v) in temp_kmers
            kmer_list[t] += v
        end
        
        #break
    end
    
    #for (t,v) in kmer_list
    #    if v<2
    #        delete!(kmer_list,t)
    #    end
    #end   
    coverage = 10
    global contig_list = []
    G = graph_creator(kmer_list,['A','C','G','T'], coverage)

    print("Starting size ",length(G),"\n")
    G,ls = compact_graph!(G,number_of_compactions)
    run_walk(G,ls)
    #print("\n new contig list\n",contig_list)
    contig_list, G
end

function read_and_walk_string(k :: Int64, number_of_compactions :: Int64, seq :: String)
    kmer_list =  DefaultDict{DNASeq, Int64}(0)
    
    #FastaReader(name) do fr
    #    
    input = string_to_DNASeq(seq)
    temp_kmers = read_kmer(input, input[1].len + (length(input)-1)*64,k)
    for (t,v) in temp_kmers
        kmer_list[t] += v
    end
    
    
    #end
        
    
    coverage = 5
    global contig_list = []
    G = graph_creator(kmer_list,['A','C','G','T'], coverage)

    print("Starting size ",length(G),"\n")
    G,ls = compact_graph!(G,number_of_compactions)
    #print("length()",length(ls))
    output = run_walk(G,ls)
    append!(contig_list,ls)
    
    for i in ls
        if !(i in contig_list)
            append!(contig_list,i)
        end
    end
    deleteat!(contig_list, findall(x->x==Any[],contig_list))

    for i in output
        if !(i in contig_list)
            append!(contig_list,i)
        end
    end
    deleteat!(contig_list, findall(x->x==Any[],contig_list))

    #print("\n new contig list\n",contig_list)
    unique(contig_list), G
end






## Test:
## k = 3
## number_of_compactions = 4
## name = "ecoli_illumina_10x_part_.fasta"
## read_and_walk(k,number_of_compactions,name)

## Test 
## output, G = read_and_walk(32,28,"ecoli_illumina_10x_part_.fasta")
## for i in output
##    for j in i
##        print(DNASeq_to_string(j))
##    end
##    print("\n\n")
## end


### String function Test
## name = randstring("ACGT",1000)
## output,G = read_and_walk_string(32,20,name)
