
include("graph_construction.jl")
include("kmer_counting_v2.jl")

function compact_graph( G :: DefaultDict{DNASeq,macro_node}, k:: Int64, phi :: Int64)
    num_mn = length(G)
    contigs = []
    G1 = copy(G)
    while(num_mn > phi)
        i_s = IS(G1,['A','C','G','T'],k)
        if isempty(i_s)
            print(i_s)
            break
        end
        pcontig_list, transfer_nodeInfo = iterate_pack(G1,I,k)
        serialize_transfer(G1,transfer_nodeInfo)
        num_mn = length(G1)
        push!(contigs,pcontig_list)
        
    end
    return G1, pcontig_list
end

@inline function succ_neigh(kmer:: DNASeq , lmer :: Vector{DNASeq})
    ## kmer is the key, lmer is the prefix
    #k::Int64, l :: Int64
    k = kmer.len + 1
    bit1 = lmer[1].bit1
    bit2 = lmer[1].bit2
    l = (length(lmer)-1)*64 + lmer[1].len
    if l >= (k-1)
        
        if l == 0 
            print("error")
            return Nothing
        end
        return kmer_seq(bit1[64-k+2:64],bit2[64-k+2:64],k-1)
    else
        rem = k-1-l
        if l == 0 
            print("error")
            return Nothing
        end
        return kmerge(kmer, lmer[1], rem, l)
    end
end

@inline function pred_neigh(kmer :: DNASeq , lmer :: Vector{DNASeq})
    ## kmer is the macro_node key(actually k-1 mer), lmer is the prefix
    k = kmer.len +1
    bit1 = lmer[1].bit1
    bit2 = lmer[1].bit2
    l = (length(lmer)-1)*64 + lmer[1].len
    if l>= (k-1)
        return kmer_seq(bit1[1:k-1],bit2[1:k-1],k-1)
    else
        rem = k-1-l
        return kmerge(lmer[1], DNASeq(kmer.bit1[64-kmer.len+1:64-kmer.len+rem],kmer.bit2[64-kmer.len+1:64-kmer.len+rem],rem), l, rem)
    end

end

function find_pred_ext(G::DefaultDict, kmer1 ::DNASeq, kmer2 :: DNASeq)
    ### find prefix of the successor(kmer2) that is connected to the suffix of kmer 1

    for (pid ,prefix) in G[kmer2].prefixes
        if length(prefix)==1 && prefix[1].len<(k-1)
            if prefix[1].bit1[end-prefix[1].len+1:end] == kmer1.bit1[end-kmer1.len+1:end - kmer1.len + prefix[1].len] && prefix[1].bit2[end-prefix[1].len+1:end] == kmer1.bit2[end-kmer1.len+1:end - kmer1.len + prefix[1].len]
                return pid, prefix
            end
        elseif prefix[end].bit1[end-(k-1)+1:end]==kmer1.bit1[end-(k-1):end] && prefix[end].bit2[end-(k-1)+1:end] == kmer1.bit2[end - (k-1):end]
            return pid,prefix
        end
    end
    
    print("nothing found")
end

function find_succ_ext(G :: DefaultDict, kmer1 ::DNASeq, kmer2 :: DNASeq)
    ### find suffix of the predecessor(kmer1) that is connected to the prefix of kmer2

    for (sid ,suffix) in G[kmer1].suffixes
        #print(kmer1)
        #print("sss",suffix[1].bit1[end-suffix[1].len+1:end],"\n")
        if length(suffix)==1 && suffix[1].len<(k-1)
            if suffix[1].bit1[end-suffix[1].len+1:end] == kmer2.bit1[end - suffix[1].len+1:end] && suffix[1].bit2[end-suffix[1].len+1:end] == kmer2.bit2[end - suffix[1].len+1: end]
                return sid, suffix
            end
        elseif suffix[end].bit1[1:k-1]==kmer2.bit1[end-(k-1):end] && suffix[end].bit2[1:k-1] == kmer2.bit2[end - (k-1):end]
            return sid,suffix
        end
    end
    
    print("nothing found")
end

@inline function symmetric_bool(var::Bool) 
    2*Int(var) -1
end

function kmerge(kmer_1::DNASeq, kmer_2::DNASeq, k::Int64, l :: Int64 )
    if l+k <=64
        bit11 = kmer_1.bit1
        bit12 = kmer_1.bit2
        bit21 = kmer_2.bit1
        bit22 = kmer_2.bit2
        return kmer_seq(vcat(bit11[end-k+1:end],bit21[end-l+1:end]), vcat(bit12[end-k+1:end],bit22[end-l+1:end]),k+l)
    else
        bit11 = kmer_1.bit1
        bit12 = kmer_1.bit2
        bit21 = kmer_2.bit1
        bit22 = kmer_2.bit2
        seq1 = kmer_seq(vcat(bit11[end-l+1:end],bit21),vcat(bit12[end-l+1:end],bit22),64)
        seq2 = kmer_seq(bit11[1:64-l],bit12[1:64-l],64-l)
        return [seq2, seq1]
    end


end

function kmerge(kmer_1::Union{Vector{DNASeq},DNASeq}, kmer_2::Union{Vector{DNASeq},DNASeq})
    if typeof(kmer_1) == Vector{DNASeq} 
        k = (length(kmer_1)-1)*64 + kmer_1[1].len
    else
        k = kmer_1.len
    end
    if typeof(kmer_2) == Vector{DNASeq}
        l = (length(kmer_2)-1)*64 + kmer_2[1].len
    else
        l = kmer_2.len
    end
    
    return kmerge(kmer_1,kmer_2,k,l)
end

function kmerge(kmer_1::Vector{DNASeq}, kmer_2::DNASeq, k::Int64, l :: Int64 )
    k1 = k%64
    res = []
    range = k>=(64-l) ? (l+1:64) : 64-k+1:64
    bit11 = kmer_1[end].bit1[range]
    bit12 = kmer_1[end].bit2[range]
    bit21 = kmer_2.bit1
    bit22 = kmer_2.bit2
   
    if k == k1
        pushfirst!(res, kmer_seq(vcat(bit11,bit21), vcat(bit12,bit22),length(range)+l))
    else
        pushfirst!(res, kmer_seq(vcat(bit11,bit21), vcat(bit12,bit22),64))
    end

   for i in length(kmer_1):-1:2

        bit11 = kmer_1[i-1].bit1[l+1:end]
        bit12 = kmer_1[i-1].bit2[l+1:end]
        bit21 = kmer_1[i].bit1[1:l]
        bit22 = kmer_1[i].bit2[1:l]
        pushfirst!(res, kmer_seq(vcat(bit11,bit21), vcat(bit12,bit22),64))
    end
    if l+k1 > 64
        bit1 = kmer_1[1].bit1[1:l]
        bit2 = kmer_1[1].bit2[1:l]
        pushfirst!(res,kmer_seq(bit1,bit2,l))
    end
    return res
end

function kmerge(kmer_1::DNASeq, kmer_2:: Vector{DNASeq}, k::Int64, l :: Int64 )
    l1 = l%64
    res = l1[2:end]
    if l1+k<=64
        bit11 = kmer_1.bit1[end-k+1:end] 
        bit12 = kmer_1.bit2[end-k+1:end]
        bit21 = kmer_2[1].bit1[end-l1+1:end]
        bit22 = kmer_2[1].bit2[end-l1+1:end]
        pushfirst1(res, kmer_seq(vcat(bit11,bit21),vcat(bit12,bit22),l1+k))
    end
    res
end

function kmerge(kmer_1::Vector{DNASeq}, kmer_2::Vector{DNASeq}, k::Int64, l :: Int64 )
    #l1 = l÷ 64
    l1_m = l%64
    #k1 = k÷64
    k1_m = k%64
    if l1_m == 0
        return vcat(kmer_1,kmer_2)
    end

    res = kmer_2[2:end]
    r = (l1_m+k1_m)>64

    bit11 = kmer_1[end].bit1[end-k1_m+1:end]
    bit12 = kmer_1[end].bit2[end-k1_m+1:end]
    bit21 = kmer_2[1].bit1[end-l1_m+1:end]
    bit22 = kmer_2[1].bit2[end-l1_m+1:end]
    
    if l+k<64
         pushfirst!(res,kmer_seq(vcat(bit11,bit21),vcat(bit12,bit22),l+k))
    
    else
        pushfirst!(res,kmer_seq(vcat(bit11,bit21),vcat(bit12,bit22),64))
    
    end

    
    
    lk1 = length(kmer_1)
    
    for i in 1:lk1-1
        ind2 = lk1 + 1 - i
        ind1 = lk1  - i 
        bit11 = kmer_1[ind1].bit1[l1_m+1:end]
        bit12 = kmer_1[ind1].bit2[l1_m+1:end]
        bit21 = kmer_1[ind2].bit1[1:l1_m]
        bit22 = kmer_1[ind2].bit2[1:l1_m]
        pushfirst!(res,kmer_seq(vcat(bit11,bit21),vcat(bit12,bit22),64))
    end
    if r
        ind1 = 1
        bit1 = kmer_1[ind1].bit1[end - k1_m + 1: l1_m]
        bit2 = kmer_1[ind1].bit2[end - k1_m + 1: l1_m]
        pushfirst!(res,kmer_seq(bit1,bit2, l1_m + k1_m - 64))
    end

    return res

end


function IS(G::DefaultDict, Alphabet:: Vector{Char}, k :: Int64)
    IS_ = Set()

    for node in keys(G)
        max_kmer = node
        for (pid, p_kmer) in G[node].prefixes
            if p_kmer[1].len == 0 ## Terminal
                continue
            end
            
            pred_node = pred_neigh(node, p_kmer)
            if pred_node>node
                max_kmer = pred_node
                break
            end
        end
        for (sid, s_kmer) in G[node].suffixes
            if s_kmer[1].len == 0 ## Terminal
                continue
            end
            succ_node = succ_neigh(node, s_kmer)
            if succ_node > node
                max_kmer = succ_node
                break
            end
        end
        if max_kmer == node
            push!(IS_, node) 
        end
    end
    return IS_

end


function iterate_pack(G :: DefaultDict, IS_ :: Set, k:: Int64)
    transfer_nodeInfo = Set{Tuple}([])
    pcontig_list = Vector{}()
    #pred_node = DNASeq()
    for node in IS_

        for (i,iset) in G[node].wire_info.prefix_info
            pid = i
            if G[node].prefix_terminal_id != pid #&& length(G[node].prefixes)>0
                pred_node = pred_neigh(node, G[node].prefixes[pid] )
                ppid,pred_ext = find_succ_ext(G,pred_node, node)
                #G[node].prefixes[pid]
                
                
            end
            
            for (j, offset, visit_count) in iset
                
                sid = j
                if G[node].prefix_terminal_id == pid && G[node].suffix_terminal_id == sid
                    contig = kmerge(kmerge(G[node].prefixes[pid],node),G[node].suffixes[pid])
                    push!(pcontig_list,contig)
                else
                    succ_node = succ_neigh(node, G[node].suffixes[sid])
                    ssid,succ_ext = find_pred_ext(G,node, succ_node)
                    #G[node].suffixes[sid]
                    
                    
                    if G[node].prefix_terminal_id != pid 
                        if G[node].suffix_terminal_id == sid
                            new_ext = pred_ext
                        else
                          new_ext = kmerge(pred_ext,succ_ext)
                        end
                        push!(transfer_nodeInfo,(pred_node,pred_ext,new_ext,visit_count,1))
                    end

                    if G[node].suffix_terminal_id != sid 
                        if G[node].prefix_terminal_id == pid
                            new_ext = succ_ext
                        else
                            new_ext = kmerge(pred_ext,succ_ext)
                        end
                        
                        push!(transfer_nodeInfo,(succ_node,succ_ext,new_ext,visit_count,0))

                        
                    end

                end
            end
            #else
            #    succ_node = succ_neigh(node, G[node].suffixes)
            #end
            
        end

    end
    
    return pcontig_list, transfer_nodeInfo
end


function serialize_transfer(G :: DefaultDict, transfer_nodeInfo :: Set)
   
    rewirelist = []
    for i in transfer_nodeInfo
        (node,ext,new_ext,visit_count,node_type) = i

        if node_type == 1
            check = 0
            for (i,p_ext) in G[node].suffixes
                if p_ext == ext 
                    check = 1
                    G[node].suffixes[i] = new_ext
                    G[node].suffix_counts[i] = (first(G[node].suffix_counts[i]),visit_count)
                end
            end
            if check == 0

                print("did not work\n")
            end
        else
            check = 0
            for (i,s_ext) in G[node].prefixes
               
                
                if s_ext == ext 
                    check = 1
                    G[node].prefixes[i] = new_ext
                    #print(G[node].suffix_counts[i])
                    G[node].prefix_counts[i] = (first(G[node].prefix_counts[i]), visit_count)
                end
            end
            if check == 0
\
                print("did not work\n")
                
            end
            push!(rewirelist,node)
        end


    end


end


## Test
## I = IS(G,['A','C','G','T'],k)
## compact_graph(G,3,2)