include("read_data.jl")
using Graphs, GraphPlot
output, G = read_and_walk(30,10,"ecoli_illumina_10x_part_.fasta")




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
