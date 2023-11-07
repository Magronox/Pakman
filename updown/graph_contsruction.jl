include("kmer_counting.jl")


Terminal = DNASeq(zeros(Bool, 64),zeros(Bool, 64),0)

VTerminal = [Terminal]

NTerminal = macro_node()
NTerminal.label.isTerminal = true

### to get the maximum values in dictionary
Base.isless(p::Pair, q::Pair) =
               isless(p.second,q.second)
           


## this function is used for wiring:


function Comp_rev(edge_)
    a(i,j) = (last(edge_[i].counts) > last(edge_[j].counts)) || ((last(edge_[i].counts) == last(edge_[j].counts)) && (first(edge_[i].counts)> first(edge_[j].counts)))
    return  a
end


function graph_creator(kmer_list :: DefaultDict, C :: Int64)
    G = DefaultDict{Int64,edge}(0)
    edge_id = 0
    mn_list = Set{macro_node}()
    vc = 0
    empty!(NTerminal.prefix_edge_ids)
    #empty!(PTerminal.suffix_edge_ids)
    push!(mn_list, NTerminal)
    #push!(mn_list, PTerminal)

    for (xkey, ~) in kmer_list
        x_prime_list = read_mn_from_kmer(xkey)
        
        for x_prime in x_prime_list

            if !([x_prime[1]] in labels(mn_list))
                #print(x_prime)
                x_prime_key,~ = x_prime
                u = macro_node(x_prime_key, Set{Int64}(), Set{Int64}())
                push!(mn_list, u)

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
                        new_edge.pred_suffix = [DNASeq([u.label.bit1[end]], [u.label.bit2[end]], 1)]
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
                        new_edge.succ_prefix = [DNASeq([u.label.bit1[1]], [u.label.bit2[1]], 1)]
    
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
                edge_id = wiring_prep(u :: macro_node, G :: DefaultDict{Int64, edge}, edge_id)
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
    pTerminal.edge_id = edge_id + 1
    pTerminal.pred_label = NTerminal
    pTerminal.succ_label = u.label
    pTerminal.isPredTerminal = true
    pTerminal.counts = (-1,-1)
    G[edge_id + 1] = pTerminal
    push!(NTerminal.suffix_edge_ids, edge_id + 1)
    edge_id += 1
    

    sTerminal = edge()
    sTerminal.edge_id = edge_id + 1
    sTerminal.succ_label = NTerminal
    sTerminal.pred_label = u.label
    sTerminal.isSuccTerminal = true
    sTerminal.counts = (-1,-1)
    G[edge_id+1] = sTerminal
    push!(NTerminal.prefix_edge_ids, edge_id + 1)
    edge_id += 1
    
    return edge_id

end



function setup_wiring!(u :: macro_node, G :: DefaultDict{Int64, edge})
    
    null_sid = 0
    null_pid = 0
    sc = 0
    pc = 0
    
    for id in u.suffix_edge_ids
        if !G[id].isSuccTerminal
            sc += last(G[id].counts)
        end
       
    end
    for id in u.prefix_edge_ids
        if !G[id].isPredTerminal
            pc += last(G[id].counts)
        end
    end

    for id in u.suffix_edge_ids
        if G[id].isSuccTerminal || G[id].isPredTerminal
            null_sid = id
            if last(G[id].counts) == -1
                G[id].counts = (1,max(pc-sc,0))
                break
            end
        end
    end

    for id in u.prefix_edge_ids
        if G[id].isPredTerminal || G[id].isSuccTerminal
            null_pid = id
            if last(G[id].counts) == -1
                G[id].counts = (1, max(sc-pc,0))
                break

            end
        end
    end
    print((null_sid in keys(G)),"\n")
    leftover = sc + last(G[null_sid].counts)

    last_largest_pid, prefix_begin_pos = -1,-1
    wire_idx = 0
    var_p, var_s = 0,0
    top_p, top_s = 0,0
    p_size = 0
    offset_in_suffix = zeros(Int64, length(u.suffix_edge_ids))

    indices_s = collect(u.suffix_edge_ids)
    indices_p = collect(u.prefix_edge_ids)
    
    indices_s = sort(indices_s, lt = Comp_rev(G))
    indices_p = sort(indices_p, lt = Comp_rev(G))
    print(indices_s,"\n", indices_p,"\n")
    @assert(false)

    @assert(sc + last(G[null_sid].counts) == pc + last(G[null_pid].counts))
    s_map = 1:length(u.suffix_edge_ids) .=> sort(u.suffix_edge_ids)
    p_map = 1:length(u.prefix_edge_ids) .=> sort(u.prefix_edge_ids)
    while leftover > 0
 
        largest_sid = indices_s[top_s + 1]-1;
        largest_pid = indices_p[top_p + 1]-1;
        count = min(last(G[p_map[largest_pid+1]].counts) - var_p,last(G[s_map[largest_sid+1]].counts) - var_s)
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

        if var_p == last(G[p_map[largest_pid + 1]].counts)
            var_p = 0
            top_p += 1
            u.prefix_begin_info[largest_pid + 1] = (prefix_begin_pos , p_size )
            p_size = 0
        end
       
        if var_s == last(G[s_map[largest_sid + 1]].counts)
            var_s = 0
            top_s += 1
        end
  
    end

end
### test
### slen = 100
### input = randstring("ACGT",slen)
### seq = string_to_DNASeq(input)
### kmer_list = read_kmer(seq, length(input),k)
### G = graph_creator(kmer_list, 5)

