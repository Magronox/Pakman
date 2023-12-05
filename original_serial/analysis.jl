include("read_data.jl")
using Graphs, GraphPlot
output, G = read_and_walk(30,10,pwd()*"/data/ecoli_illumina_10x_part_.fasta")



"""
function visual_graph(G :: DefaultDict{Vector{DNASeq}})
    VG = SimpleDiGraph(0)
    indices = Dict()
    counter = 1
    labels = []
    for (i,j) in G
        indices[i] = counter
        counter += 1
        add_vertex!(g)
    end
    for (i,j) in G
        for (ii,jj) in j.prefixes
            if jj == VTerminal
                continue
            end
            indices[] = 0
        end
    end

    for (i,j) in G
        for (ii,jj) in j.suffixes
            if jj == VTerminal
                continue
            end
            succ_node = succ_neigh(i,jj)
            if succ_node in keys(G)
                add_edge!(VG, indices[i])
            end
        end


    end
    

end"""



function show_graph(G :: DefaultDict{Vector{DNASeq}, macro_node})
    for (i,j) in G
        print("Node ", DNASeq_to_string(i[1]),"\nPrefixes ")
        for (ii,jj) in j.prefixes
            print(DNASeq_to_string(jj[1])," ",j.prefixes_terminal[ii]," ",j.prefix_counts[ii],"\t")
        end
        print("\nSuffixes ")
        for (ii,jj) in j.suffixes
            print(DNASeq_to_string(jj[1])," ",j.suffixes_terminal[ii]," ", j.suffix_counts[ii],"\t")
        end
        print("\n")
    end
end


function show_kmers(file :: String, k :: Int64)
    kmer_list =  DefaultDict{DNASeq, Int64}(0)
    for (~,seq) in FastaReader(file)
        input = string_to_DNASeq(seq)
        #print(seq,"\n")
        temp_kmers = read_kmer(input, input[1].len + (length(input)-1)*64,k)
        for (t,v) in temp_kmers
            kmer_list[t] += v
        end
    end
    for(i,j) in kmer_list
        print(seq_to_int64(i),"\t",j,"\n")

    end
end
function show_wiring(G ::  DefaultDict{Vector{DNASeq}, macro_node})
    for (i,j) in G  
        print("node\t",DNASeq_to_string(i[1]),"\n")
        print("wire info\t", j.wire_info,"prefix_begin_info\t",j.prefix_begin_info,"\n")
    end
end

function seq_to_int64(seq,s=Int64(0))
    v = 1
    #for i in view(arr,length(arr):-1:1)
    j = 1
    for i in 64:-1:64-seq.len+1
        bit = seq.bit2[i]
        s += v*bit
        v <<= 1
        bit = seq.bit1[i]
        s += v*bit
        v <<= 1

    end 
    s
end


function check_suf_pref(G :: DefaultDict{Vector{DNASeq}})
    
    for (i,j) in G
        for (ii,jj) in j.suffixes
            if jj == VTerminal
                continue
            end
            succ_node = succ_neigh(i[1],jj)
            if !(succ_node in keys(G))
                continue
            end
            a = find_succ_ext(G,i,succ_node)
            if a == Nothing
                print(G[succ_node].prefixes_terminal[ii])
                print("node    ",DNASeq_to_string(i[1]),"\n")
                print("succ_node    ", DNASeq_to_string(succ_node[1]),"\n")
                print("suffix    ",DNASeq_to_string(jj[end]),"\n")
                print("error\n")
                return false
            end
        end
        for (ii,jj) in j.prefixes
            if jj == VTerminal
                continue
            end
            pref_node = pred_neigh(i[1],jj)
            if !(pref_node in keys(G))
                continue
            end
            a = find_pred_ext(G, pref_node, i)
            if a == Nothing
                print(suffixes_terminal[ii],"\n")
                print("error\n")
                print(DNASeq_to_string(pref_node[1]),"\n")
                print(DNASeq_to_string(i[1]),"\n")
                print(DNASeq_to_string(jj[1]),"\n")
                return false
            end
        end
    end
    return true
end

