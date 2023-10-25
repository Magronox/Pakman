include("kmer_counting.jl")


Terminal = DNASeq(zeros(Bool, 64),zeros(Bool, 64),0)

VTerminal = [Terminal]

STerminal = macro_node()
STerminal.isTerminal = true
PTerminal = macro_node()
PTerminal.isTerminal = true

### to get the maximum values in dictionary
Base.isless(p::Pair, q::Pair) =
               isless(p.second,q.second)
           


## this function is used for wiring:


function Comp_rev(counts)
    a(i,j) = (last(counts[i]) > last(counts[j])) || ((last(counts[i]) == last(counts[j])) && (first(counts[i])> first(counts[j])))
    return  a
end


function graph_creator(kmer_list :: DefaultDict, C :: Int64)
    G = DefaultDict{Int64,edge}(0)
    edge_id = 0
    mn_list = Set{mn_label}()
    vc = 0
    empty!(STerminal.prefix_edge_ids)
    empty!(PTerminal.suffix_edge_ids)

    for (xkey, ~) in kmer_list
        x_prime_list = read_mn_from_kmer(xkey)
        
        for x_prime in x_prime_list

            if !([x_prime[1]] in mn_list)
                print(x_prime)
                x_prime_key,~ = x_prime
                u = macro_node(x_prime_key, Set{Int64}(), Set{Int64}())
                push!(mn_list, u.label)

                #u.label = [x_prime_key]
                #u.id = id
                pid = 1
                sid = 1


                for c in alphabet
                    
                    temp = kmerge(c,x_prime_key, "from left")
                    node = [DNASeq(vcat(MVector{64-k+1, Bool}(zeros(Bool, 64-k+1)),temp.bit1[end-k+1:end-1]),vcat(MVector{64-k+1, Bool}(zeros(Bool, 64-k+1)),temp.bit2[end-k+1:end-1]),k-1)]

                    if  node!= [x_prime_key] && temp in keys(kmer_list)
                        new_edge = edge()
                        new_edge.pred_label = neigh_label(c, u, "pred")
                        new_edge.succ_label = u
                        new_edge.pred_suffix = DNASeq([u.label.bit1[end]], [u.label.bit2[end]], 1)
                        new_edge.succ_prefix = string_to_DNASeq(string(c))

                        if !(new_edge in values(G))
                            new_edge.edge_id  = edge_id + 1
                            u.prefix_edge_ids = edge_id + 1
                            edge_id += 1
                            vc = ceil(Int64,kmer_list[temp]/C)
                            new_edge.counts = (kmer_list[temp],vc)
                            G[edge_id] = new_edge
                        end

                    end

                    temp = kmerge(c,x_prime_key, "from right")
                    node = [DNASeq(vcat(MVector{64-k+1}(zeros(Bool,64-k+1)),temp.bit1[end-k+2:end]),vcat(MVector{64-k+1}(zeros(Bool, 64-k+1)),temp.bit2[end-k+2:end]),k-1)]

                    if node != [x_prime_key] && temp in keys(kmer_list)
                        new_edge = edge()
                        new_edge.succ_label = neigh_label(c, u, "succ")
                        new_edge.pred_label = u
                        new_edge.pred_suffix = string_to_DNASeq(string(c))
                        new_edge.succ_prefix = DNASeq([u.label.bit1[1]], [u.label.bit2[1]], 1)
                        if !(new_edge in values(G))
                            new_edge.edge_id  = edge_id + 1
                            u.suffix_edge_ids = edge_id + 1
                            edge_id += 1
                            vc = ceil(Int64,kmer_list[temp]/C)
                            new_edge.counts = (kmer_list[temp],vc)
                            G[edge_id] = new_edge
                        end
                    end
                end
                wiring_prep(u :: macro_node, G :: DefaultDict{Int64, edge}, edge_id)
                setup_wiring!(u :: macro_node, G :: DefaultDict{Int64, edge})
                #G[u.label] = u
            end
            
            
            
            #u = macro_node()
        end
        
        
        

    end
#    pref = macro_node()
#    pref.prefix_terminal = true
#    pref.suffix_terminal = false
#    suff = macro_node()
#    suff.prefix_terminal = false
#    suff.suffix_terminal = true

    G, mn_list

end

function wiring_prep(u :: macro_node, G :: DefaultDict{Int64, edge}, edge_id :: Int64)
    pTerminal = edge()
    pTerminal.pred_label = PTerminal
    pTerminal.succ_label = u.label
    pTerminal.isPredTerminal = true
    pTerminal.counts = (-1,-1)
    push!(PTerminal.suffix_edge_ids, edge_id + 1)
    edge_id += 1
    sTerminal = edge()
    sTerminal.succ_label = STerminal
    sTerminal.pred_label = u.label
    sTerminal.isSuccTerminal = true
    sTerminal.counts = (-1,-1)
    push!(STerminal.prefix_edge_ids, edge_id + 1)
    edge_id += 1
    return edge_id

end


"""
function setup_wiring!(u :: macro_node,)
    sc, pc = 0,0
    null_sid, null_pid = -1,-1
    

    for id in u.suffixes
        if u.suffixes[i] != VTerminal
            sc += last(u.suffix_counts[i])
        end
       
    end
    for (i,~) in u.prefixes
        if u.prefixes[i] != VTerminal
            pc += last(u.prefix_counts[i])
        end
    end

    for (i,j) in u.suffixes
        if length(j) == 1 && j[1].len == 0
            null_sid = i
            if last(u.suffix_counts[i]) == -1
                 u.suffix_counts[i] = (1,max(pc-sc,0))
                 break
            end
        end
    end

    for (i,j) in u.prefixes
        if length(j) == 1 && j[1].len == 0
            null_pid = i
            if last(u.prefix_counts[i]) == -1
                u.prefix_counts[i] = (1, max(sc-pc,0))
                break
                @assert(null_pid == length(u.prefixes))

            end
        end
    end

    leftover = sc + last(u.suffix_counts[null_sid])

    last_largest_pid, prefix_begin_pos = -1,-1
    wire_idx = 0
    var_p, var_s = 0,0
    top_p, top_s = 0,0
    p_size = 0
    offset_in_suffix = zeros(Int64, length(u.suffixes))

    indices_s = collect(1:length(keys(u.suffixes)))
    indices_p = collect(1:length(keys(u.prefixes)))
    
    indices_s = sort(indices_s, lt = Comp_rev(u.suffix_counts))
    indices_p = sort(indices_p, lt = Comp_rev(u.prefix_counts))

    @assert(sc + last(u.suffix_counts[null_sid]) == pc + last(u.prefix_counts[null_pid]))
   
    while leftover > 0
 
        largest_sid = indices_s[top_s + 1]-1;
        largest_pid = indices_p[top_p + 1]-1;
        count = min(last(u.prefix_counts[largest_pid+1]) - var_p,last(u.suffix_counts[largest_sid+1]) - var_s)
        push!(u.wire_info[wire_idx + 1], (largest_sid+1, offset_in_suffix[largest_sid + 1],count))
        if last_largest_pid != largest_pid
            prefix_begin_pos = wire_idx 
            last_largest_pid = largest_pid
        end
        p_size += 1
        wire_idx += 1
        leftover -= count
        var_p += count
        var_s += count
        offset_in_suffix[largest_sid + 1] += count

        if var_p == last(u.prefix_counts[largest_pid + 1])
            var_p = 0
            top_p += 1
            u.prefix_begin_info[largest_pid + 1] = (prefix_begin_pos , p_size )
            p_size = 0
        end
       
        if var_s == last(u.suffix_counts[largest_sid + 1])
            var_s = 0
            top_s += 1
        end
  
    end

end
### test
### G = graph_creator(kmer_list,['A','C','G','T'], 5)

"""
