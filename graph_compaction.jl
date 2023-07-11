
include("graph_construction.jl")
include("kmer_counting_v2.jl")

function compact_graph!( G :: DefaultDict{Vector{DNASeq},macro_node}, k:: Int64, phi :: Int64)
    num_mn = length(G)
    contigs = []
    pcontig_list= []

    while(num_mn > phi)
        i_s = IS(G)
        
        if isempty(i_s)
            print("error\n")
            break
        end

        pcontig_list, transfer_nodeInfo = iterate_pack(G,i_s,k)

        
        serialize_transfer!(G,transfer_nodeInfo)


        
        for i in i_s
            
            delete!(G,i)
        end
        num_mn = length(G)
        push!(contigs,pcontig_list)

        print("loop\n",num_mn,"\n")
        
    end
    return G, pcontig_list
end

@inline function succ_neigh(kmer:: DNASeq , lmer :: Vector{DNASeq})
    ## kmer is the key, lmer is the suffix
    #k::Int64, l :: Int64
    k = kmer.len + 1
    bit1 = lmer[end].bit1
    bit2 = lmer[end].bit2

    l = 0
    for i in lmer
        l += i.len
    end
    if l >= (k-1)
        
        if l == 0 
            print("error")
            return Nothing
        end
        #print("chera\n")
        #print("inja?",DNASeq_to_string(lmer[1]),"\n\n",kmer_seq(bit1[64-k+2:64],bit2[64-k+2:64],k-1),"\n\n\n")
        if kmer_seq(bit1[64-k+2:64],bit2[64-k+2:64],k-1)!= kmer
            if !([kmer_seq(bit1[64-k+2:64],bit2[64-k+2:64],k-1)] in keys(G))
                return Nothing
            end
            return [kmer_seq(bit1[64-k+2:64],bit2[64-k+2:64],k-1)]
        else return Nothing
        end
    else
        rem = k-1-l
        if l == 0 
            print("error")
            return Nothing
        end
        #print("aziatnakon",DNASeq_to_string(kmer),"\n",DNASeq_to_string(lmer[1]),"\n",rem,"  ", l ,"\n\n",kmerge(kmer, lmer[1], rem, l),"\n\n\n")
        if kmerge(kmer,lmer[1],rem,l)!= kmer
            if !([kmerge(kmer, lmer[1], rem, l)] in keys(G))
                return Nothing
            end
            return [kmerge(kmer, lmer[1], rem, l)]
        else return Nothing
        end
    end
end

@inline function pred_neigh(kmer :: DNASeq , lmer :: Vector{DNASeq})
    ## kmer is the macro_node key(actually k-1 mer), lmer is the prefix
    k = kmer.len +1
    bit1 = lmer[1].bit1
    bit2 = lmer[1].bit2
    l = 0
    for i in lmer
        l+= i.len
    end


    ll = lmer[1].len
    
    
    if l>= (k-1)
        if ll>=(k-1)
            if kmer_seq(bit1[64-ll+1:64-ll+ k-1],bit2[64-ll+1:64-ll+k-1],k-1) != kmer
                if !([kmer_seq(bit1[64-ll+1:64-ll+ k-1],bit2[64-ll+1:64-ll+k-1],k-1)] in keys(G))
                    return Nothing
                end
                return [kmer_seq(bit1[64-ll+1:64-ll+ k-1],bit2[64-ll+1:64-ll+k-1],k-1)]
            else
                return Nothing
            end
            
        else
            
            rem = k-1 - ll
            if  kmer_seq(vcat(bit1[64-ll+1:64],lmer[2].bit1[1:rem]),vcat(bit2[64-ll+1:64],lmer[2].bit2[1:rem]),k-1) != kmer
                if !([kmer_seq(vcat(bit1[64-ll+1:64],lmer[2].bit1[1:rem]),vcat(bit2[64-ll+1:64],lmer[2].bit2[1:rem]),k-1)] in keys(G))
                    return Nothing
                end
                return [kmer_seq(vcat(bit1[64-ll+1:64],lmer[2].bit1[1:rem]),vcat(bit2[64-ll+1:64],lmer[2].bit2[1:rem]),k-1)]
            else
                return Nothing
            end
        end
    else
        rem = k-1-l
        if kmerge(lmer[1], DNASeq(kmer.bit1[64-kmer.len+1:64-kmer.len+rem],kmer.bit2[64-kmer.len+1:64-kmer.len+rem],rem), l, rem) != kmer
            if !([kmerge(lmer[1], DNASeq(kmer.bit1[64-kmer.len+1:64-kmer.len+rem],kmer.bit2[64-kmer.len+1:64-kmer.len+rem],rem), l, rem)] in keys(G))
                return Nothing
            end
            
            return [kmerge(lmer[1], DNASeq(kmer.bit1[64-k+2:64-l],kmer.bit2[64-kmer.len+1:64-l],rem), l, rem)]
        else
            return Nothing
        end
    end

end

function find_succ_ext(G::DefaultDict, kmer1 ::DNASeq, kmer2 :: DNASeq)
    ### find prefix of the successor(kmer2) that is connected to the suffix of kmer 1
    k = kmer1.len +1
   
    for (pid ,prefix) in G[[kmer2]].prefixes
        if DNASeq_to_string(kmer1)=="CT"
        
            
        end
        if prefix[1] == Terminal
            continue
        elseif length(prefix)==1 && prefix[1].len<(k-1)
            
            if prefix[1].bit1[end-prefix[1].len+1:end] == kmer1.bit1[end-kmer1.len+1:end - kmer1.len + prefix[1].len] && prefix[1].bit2[end-prefix[1].len+1:end] == kmer1.bit2[end-kmer1.len+1:end - kmer1.len + prefix[1].len]
                return pid, prefix
            end
            
        elseif prefix[1].bit1[end-prefix[1].len+1:end-prefix[1].len+k-1]==kmer1.bit1[end-(k-1)+1:end] && prefix[1].bit2[end-prefix[1].len + 1:end-prefix[1].len+k-1] == kmer1.bit2[end - (k-1)+1:end]
            
            return pid,prefix
        end
    end
   
    print("nothing found")
    return Nothing
end

function find_pred_ext(G :: DefaultDict, kmer1 ::DNASeq, kmer2 :: DNASeq)
    ### find suffix of the predecessor(kmer1) that is connected to the prefix of kmer2
    k = kmer1.len +1
    for (sid ,suffix) in G[[kmer1]].suffixes
        #print(suffix[1].len,suffix[1],suffix[1]," ",k-1,"  ", length(suffix), DNASeq_to_string(suffix[1]),"\n")
        if length(suffix)==1 && suffix[1].len<(k-1)
            if suffix[1].bit1[end-suffix[1].len+1:end] == kmer2.bit1[end - suffix[1].len+1:end] && suffix[1].bit2[end-suffix[1].len+1:end] == kmer2.bit2[end - suffix[1].len+1: end]
                return sid, suffix
            end
        elseif suffix[end].bit1[end-k+2:end]==kmer2.bit1[end-(k-1)+1:end] && suffix[end].bit2[end-k+2:end] == kmer2.bit2[end - (k-1)+1:end]
            return sid,suffix
        end
    end

    
   
    print("nothing found")
    return Nothing
end

@inline function symmetric_bool(var::Bool) 
    2*Int(var) -1
end

function kmerge(kmer_1::DNASeq, kmer_2::DNASeq, k::Int64, l :: Int64 )
    if kmer_1 == Terminal 
        if kmer_2 == Terminal
            print("error")
        end
        return [kmer_1, Terminal]
    elseif kmer_2 == Terminal
        return [kmer_1,Terminal]
            
    end
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
    bit21 = kmer_2.bit1[64-l+1:end]
    bit22 = kmer_2.bit2[64-l+1:end]
   
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
    res = Vector{DNASeq}(res)
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
    if l1_m==0
        l1_m=64
        if l == 0
            return vcat(kmer_1,kmer_2)
        end
        
    end

    #k1 = k÷64
    k1_m = k%64
    if k1_m == 0 
        if k>0
            k1_m = 64
        elseif l==0
            print("error\n")
        else
            return vcat(kmer_1,kmer_2)
        end
        
    end



    res = kmer_2[2:end]
    r = (l1_m+k1_m)>64

    if (l+k1_m)<64

        bit11 = kmer_1[end].bit1[end-k1_m+1:end]
        bit12 = kmer_1[end].bit2[end-k1_m+1:end]
        bit21 = kmer_2[1].bit1[end-l1_m+1:end]
        bit22 = kmer_2[1].bit2[end-l1_m+1:end]
        
        
        pushfirst!(res,kmer_seq(vcat(bit11,bit21),vcat(bit12,bit22),l+k1_m))
    
    else

    bit11 = kmer_1[end].bit1[l1_m+1:end]
    bit12 = kmer_1[end].bit2[l1_m+1:end]
    bit21 = kmer_2[1].bit1[end-l1_m+1:end]
    bit22 = kmer_2[1].bit2[end-l1_m+1:end]
    
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


function IS(G::DefaultDict)#, Alphabet:: Vector{Char}, k :: Int64)
    IS_ = Set()

    for node in keys(G)
        max_kmer = node
        print("node",DNASeq_to_string(node[1]),"\n\n")
        for (pid, p_kmer) in G[node].prefixes
            if p_kmer[1].len == 0 ## Terminal
                continue
            end
            print(DNASeq_to_string(p_kmer[1]),"\n")
            
            pred_node = pred_neigh(node[1], p_kmer)
            if pred_node == Nothing
                continue
            end

            
            if pred_node>node
                max_kmer = pred_node
                
                break
            end
        end
        for (sid, s_kmer) in G[node].suffixes
            
            if s_kmer[end].len == 0 ## Terminal
                continue
            end
            succ_node = succ_neigh(node[1], s_kmer)
            if succ_node == Nothing
                continue
            end
            
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
    transfer_nodeInfo = DefaultDict{Vector{DNASeq}, Set{Tuple}}(Set{Tuple}())
    pcontig_list = Vector{}()
    
    
    for node in IS_


        for (i,iset) in G[node].wire_info.prefix_info
            pid = i
            if !(pid in G[node].prefix_terminal_id)  #&& length(G[node].prefixes)>0
                pred_node = pred_neigh(node[1], G[node].prefixes[pid] )
                if pred_node == Nothing
                    continue
                end
                
                print(find_pred_ext(G,pred_node[1], node[1]),"\n")
                ssid,pred_ext = find_pred_ext(G,pred_node[1], node[1])
 
            
                
                
            end
            
            for (j, offset, visit_count) in iset
                sid = j
                if (pid in G[node].prefix_terminal_id)  && (sid in G[node].suffix_terminal_id)
                    print("\n\n\n\npleeease\n\n\n")
                    
            
                    contig = kmerge(kmerge(G[node].prefixes[pid],node),G[node].suffixes[pid])
                    push!(pcontig_list,contig)
                else
                    if G[node].suffixes[sid] == VTerminal || G[node].suffixes[sid] == VTerminal
                        succ_node = Nothing
                        succ_ext = Nothing
                        
                        
                    else
                        succ_node = succ_neigh(node[1], G[node].suffixes[sid])

                        if succ_node == Nothing
                            continue
                        
                        end
                        ppid,succ_ext = find_succ_ext(G,node[1], succ_node[1])
                    end
                    
                    if !(pid in G[node].prefix_terminal_id) && pred_ext!= VTerminal ## I added the second condition but seemed necessary

                        new_ext = kmerge(pred_ext,G[node].suffixes[sid])
                        if length(new_ext)>1 && new_ext[end-1] == Terminal
                            new_ext = new_ext[1:end-1]
                        end
                        push!(transfer_nodeInfo[node] , (pred_node,pred_ext,new_ext,visit_count,1,ssid,sid))
                    end

                    if !(sid in G[node].suffix_terminal_id)
                        
                        #print("jafar",G[node].prefixes[pid][1],"\n\n\n")
                        new_ext = kmerge(G[node].prefixes[pid],succ_ext)
                        if length(new_ext)>1 && new_ext[2] == Terminal
                            new_ext = new_ext[2:end]
                           
                        end
                        
                        push!(transfer_nodeInfo[node] , (succ_node,succ_ext,new_ext,visit_count,0,ppid,pid))

                        
                    end

                end
            end
        end

    end
    
    return pcontig_list, transfer_nodeInfo
end


function serialize_transfer!(G :: DefaultDict, transfer_nodeInfo :: DefaultDict{Vector{DNASeq}, Set{Tuple}, Set{Tuple}})
   
    rewirelist = []
    #print("transfer",transfer_nodeInfo,"\n\n\n\n")

    for node in keys(transfer_nodeInfo)
       


        index_suf = 0
        index_pref = 0
        
        for i in transfer_nodeInfo[node]
            (nnode,ext,new_ext,visit_count,node_type,id,pid) = i
            
            if node_type == 1
                
                
                if index_suf == 0 && !(new_ext in values(G[nnode].suffixes))
                    if new_ext[end] == Terminal
                        push!(G[nnode].suffix_terminal_id,id)
                    end
                #    print("prev_suff",G[node].suffixes[id],"\n")
                    G[nnode].suffixes[id] = new_ext
                #    print("new_suff",DNASeq_to_string(new_ext[1]))
                    
                    G[nnode].suffix_counts[id] = (first(G[nnode].suffix_counts[id]),visit_count)
                    index_suf += 1
                elseif !(new_ext in values(G[nnode].suffixes))
                #    print("prev_suff",DNASeq_to_string(G[node].suffixes[id][1]),"\n\n")
                #    print("new_ext",DNASeq_to_string(new_ext[1]))
                    
                    G[nnode].suffixes[1 + maximum(keys(G[nnode].suffixes))] = new_ext
                    push!(G[nnode].wire_info.prefix_info[pid],(maximum(keys(G[nnode].suffixes)),0,visit_count))
                    G[nnode].suffix_counts[maximum(keys(G[nnode].suffixes))] = (first(G[nnode].suffix_counts[id]),visit_count)
                    if new_ext[end] == Terminal
                        push!(G[nnode].suffix_terminal_id,maximum(keys(G[nnode].suffixes)))
                    end
                end

            else
                #print("prev_pref",DNASeq_to_string(G[node].prefixes[id][1]),"\n")
                #print("new ext",DNASeq_to_string(new_ext[1]),"\n\n\n")
                

                if index_pref == 0 && !(new_ext in values(G[nnode].prefixes))
                    if new_ext[end] == Terminal
                        push!(G[nnode].prefix_terminal_id,id)
                    end
                    G[nnode].prefixes[id] = new_ext
                    G[nnode].prefix_counts[id] = (first(G[nnode].prefix_counts[id]), visit_count)
                    index_pref += 1

                elseif !(new_ext in values(G[nnode].prefixes))
                    
                    G[nnode].prefixes[1+maximum(keys(G[nnode].prefixes))] = new_ext
                    G[nnode].prefix_counts[maximum(keys(G[nnode].prefixes))] = (first(G[nnode].prefix_counts[id]), visit_count)
                    push!(G[nnode].wire_info.prefix_info[pid],(maximum(keys(G[nnode].suffixes)),1,visit_count))
                    if new_ext[end] == Terminal
                        push!(G[nnode].prefix_terminal_id,maximum(keys(G[nnode].suffixes)))
                    end
                end
                push!(rewirelist,node)
            end
    
    
        end
    

    end

end


## Test
## slen = 1000
## input = randstring("ACGT",slen)

## DNA_seq = string_to_DNASeq(input)
## k=3
## kmer_list = read_kmer(DNA_seq, length(input),k)

## G = graph_creator(kmer_list,['A','C','G','T'], 5)
## I = IS(G,['A','C','G','T'],k)
## G_new, ~ = compact_graph!(G,k,floor(Int64,length(G)*4/5))