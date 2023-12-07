using FastaIO
include("graph_walk.jl")
coverage = 100
name = "/data/random_output_500_100_10.fasta"

function read_and_walk(k :: Int64, number_of_compactions :: Int64, name :: String, analysis :: Bool = false)
    kmer_list =  DefaultDict{DNASeq, Int64}(0)
    for (~,seq) in FastaReader(name)
        input = string_to_DNASeq(seq)
        #print(seq,"\n")
        temp_kmers = read_kmer(input, input[1].len + (length(input)-1)*64,k)
        for (t,v) in temp_kmers
            kmer_list[t] += v
        end

    end

    
    global output = Set()
    G = graph_creator(kmer_list,['A','C','G','T'], coverage)

    print("Starting size ",length(G),"\n")
    G,ls = compact_graph!(G,number_of_compactions)
    output = run_walk(G,ls,analysis)
    
        
    contig = Vector{DNASeq}()
    no = -1000
    for i in output
        if ((length(i)-1)*64+i[1].len)>no
            contig = i
            no = (length(i)-1)*64+i[1].len
        end
    end  
    contig, G
end


function read_and_walk_string(k :: Int64, number_of_compactions :: Int64, seq :: String, analysis :: Bool = false)
    kmer_list =  DefaultDict{DNASeq, Int64}(0)
    
    #FastaReader(name) do fr
    #    
    input = string_to_DNASeq(seq)
    temp_kmers = read_kmer(input, input[1].len + (length(input)-1)*64,k)
    for (t,v) in temp_kmers
        kmer_list[t] += v
    end
    
    
    #end
    
    global output = Set()
    G = graph_creator(kmer_list,['A','C','G','T'], coverage)

    print("Starting size ",length(G),"\n")
    G,ls = compact_graph!(G,number_of_compactions)
    #print("length()",length(ls))
    output = run_walk(G,ls, analysis)
    contig = Vector{DNASeq}()
    no = -1000
    for i in output
        if ((length(i)-1)*64+i[1].len)>no
            contig = i
            no = (length(i)-1)*64+i[1].len
        end
    end
    output, G
end



function print_output(output)

    for i in output
        if typeof(i)==DNASeq
            #print("1\n")
            print(DNASeq_to_string(i))
        #elseif typeof(i)==Vector{DNASeq}
        #    print("2\n")
        #    for j in i
        #        print(DNASeq_to_string(j),"\n")
        #    end
        else
            print("\n")
            for j in i
                #for jj in j
                #    print(DNASeq_to_string(jj))
                #end
                print(DNASeq_to_string(j))
            end
            print("\n")
        end
    end

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
