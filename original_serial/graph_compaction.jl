
include("graph_construction.jl")
include("kmer_counting_v2.jl")

##coverage
C = 5
##========


function compact_graph!( G :: DefaultDict{Vector{DNASeq},macro_node}, compaction_times :: Int64)
    num_mn = length(G)
    #contigs = []
    pcontig_list= []
    temp = []
    #transfer_nodeInfo = []

    for t in 1:compaction_times
        #print("second\n",G[string_to_DNASeq("CAAATAGGGCGTGGTCCACACAATATCCGCC")].prefixes,"\n")
        #print("first\n",G[string_to_DNASeq("CGCAAATAGGGCGTGGTCCACACAATATCCG")].suffixes,"\n")
        #print("third\n", G[string_to_DNASeq("AAACTTCAACTCCCAAGCAACCAATTTATAT")].prefixes,"\n")
        #print("fourth\n",G[string_to_DNASeq("CTAAACTTCAACTCCCAAGCAACCAATTTAT")].suffixes,"\n")
        """print("fifth\n")
        print(G[string_to_DNASeq("ATCCGCCCTTTTTATTTAAGAAGCATAGAGG")].wire_info,G[string_to_DNASeq("ATCCGCCCTTTTTATTTAAGAAGCATAGAGG")].prefix_begin_info,"\n",G[string_to_DNASeq("ATCCGCCCTTTTTATTTAAGAAGCATAGAGG")].suffixes,"\n",G[string_to_DNASeq("ATCCGCCCTTTTTATTTAAGAAGCATAGAGG")].prefixes,"\n\n\n")
        for (j,i) in G[string_to_DNASeq("ATCCGCCCTTTTTATTTAAGAAGCATAGAGG")].suffixes
            if (succ_neigh(string_to_DNASeq("ATCCGCCCTTTTTATTTAAGAAGCATAGAGG")[1],i) in keys(G))
                print(G[succ_neigh(string_to_DNASeq("ATCCGCCCTTTTTATTTAAGAAGCATAGAGG")[1],i)].wire_info,"\n",i,"\n",DNASeq_to_string(succ_neigh(string_to_DNASeq("ATCCGCCCTTTTTATTTAAGAAGCATAGAGG")[1],i)[1]))
                print(G[succ_neigh(string_to_DNASeq("ATCCGCCCTTTTTATTTAAGAAGCATAGAGG")[1],i)].prefixes,"\n",G[succ_neigh(string_to_DNASeq("ATCCGCCCTTTTTATTTAAGAAGCATAGAGG")[1],i)].suffixes)
            else
                print("shefte\n",succ_neigh(string_to_DNASeq("ATCCGCCCTTTTTATTTAAGAAGCATAGAGG")[1],i),"\n",i,"\n")
            end
            print(succ_neigh(string_to_DNASeq("ATCCGCCCTTTTTATTTAAGAAGCATAGAGG")[1],i))
            print("kefte\n\n\n")
            #print(G[succ_neigh(string_to_DNASeq("ATCCGCCCTTTTTATTTAAGAAGCATAGAGG")[1],i)].wire_info,"\n")
        end
        """
        
        i_s = IS(G)
        if length(G)==length(i_s)
            break
        end
        ### sanity check
        
        if isempty(i_s)
            print("empty independent set error\n")
            break
        end

        temp, transfer_nodeInfo = iterate_pack(G,i_s)
        
        
        #for (i,j) in transfer_nodeInfo
        #    print(DNASeq_to_string(i[1]),"\n")
        #    print(j,"\n\n")
        #end

        push!(pcontig_list,temp)
        
        
        rewire_list = serialize_transfer!(G,transfer_nodeInfo)
        #print("rew",rewire_list,"\n")
        
        for i in rewire_list
            if i in keys(G)
            #print(i)
                G[i].wire_info = DefaultDict{Int64, Set{Tuple}}(Set())
                G[i].prefix_begin_info = DefaultDict{Int64,Tuple}((-1,-1))
            end
        end
        
        
        for i in i_s
            delete!(G,i)
        end

    
        num_mn = length(G)
        #append!(contigs,pcontig_list)
        
        for i in rewire_list
            if (i in keys(G))
                empty!(G[i].wire_info)
                empty!(G[i].prefix_begin_info)
                #initiate_wiring!(G[i])
                setup_wiring!(G[i])
            end
        end
        print("loop\nsize :",num_mn,"\n", "rewire list size : ",length(rewire_list),"\n") 
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
    if l == 0 || (lmer[end] == Terminal && lmer[end-1] != Terminal)
        return Nothing
    end

    if l >= (k-1)
        
        """
        if kmer_seq(bit1[64-k+2:64],bit2[64-k+2:64],k-1)!= kmer
            if !([kmer_seq(bit1[64-k+2:64],bit2[64-k+2:64],k-1)] in keys(G))
                return Nothing
            end
            return [kmer_seq(bit1[64-k+2:64],bit2[64-k+2:64],k-1)]
        else return Nothing
        end
        """
        return [kmer_seq(bit1[64-k+2:64],bit2[64-k+2:64],k-1)]
    else
        rem = k-1-l
        """
        if kmerge(kmer,lmer[1],rem,l)!= kmer
            if !([kmerge(kmer, lmer[1], rem, l)] in keys(G))
                return Nothing
            end
            return [kmerge(kmer, lmer[1], rem, l)]
        else return Nothing
        end
        """
        #if !([kmerge(kmer, lmer[1], rem, l)] in keys(G))
        #    print("error\n")
        #    return Nothing
        #end
        return kmerge(kmer, lmer[1], rem, l)
    end
end

@inline function pred_neigh(kmer :: DNASeq , lmer :: Vector{DNASeq})
    ## kmer is the macro_node key(actually k-1 mer), lmer is the prefix
    k = kmer.len +1
    bit1 = lmer[1].bit1
    bit2 = lmer[1].bit2
    l = 0
    for i in lmer
        l += i.len
    end


    ll = lmer[1].len
    if ll==0
        print("possible issue\n")
        return kmer
    end
    
    if l>= (k-1)
        if ll>=(k-1)
            """
            if kmer_seq(bit1[64-ll+1:64-ll+ k-1],bit2[64-ll+1:64-ll+k-1],k-1) != kmer
                if !([kmer_seq(bit1[64-ll+1:64-ll+ k-1],bit2[64-ll+1:64-ll+k-1],k-1)] in keys(G))
                    return Nothing
                end
                return [kmer_seq(bit1[64-ll+1:64-ll+ k-1],bit2[64-ll+1:64-ll+k-1],k-1)]
            else
                return Nothing
            end
            """
            return [kmer_seq(bit1[64-ll+1:64-ll+ k-1],bit2[64-ll+1:64-ll+k-1],k-1)]
        else
            
            rem = k-1 - ll
            """
            if  kmer_seq(vcat(bit1[64-ll+1:64],lmer[2].bit1[1:rem]),vcat(bit2[64-ll+1:64],lmer[2].bit2[1:rem]),k-1) != kmer
                if !([kmer_seq(vcat(bit1[64-ll+1:64],lmer[2].bit1[1:rem]),vcat(bit2[64-ll+1:64],lmer[2].bit2[1:rem]),k-1)] in keys(G))
                    return Nothing
                end
                return [kmer_seq(vcat(bit1[64-ll+1:64],lmer[2].bit1[1:rem]),vcat(bit2[64-ll+1:64],lmer[2].bit2[1:rem]),k-1)]
            else
                return Nothing
            end
            """
            if length(lmer) > 1
                return [kmer_seq(vcat(bit1[64-ll+1:64],lmer[2].bit1[1:rem]),vcat(bit2[64-ll+1:64],lmer[2].bit2[1:rem]),k-1)]
            else
                print("error\n")
            end
        end
    else
        rem = k-1-l
        """
        if kmerge(lmer[1], DNASeq(kmer.bit1[64-kmer.len+1:64-kmer.len+rem],kmer.bit2[64-kmer.len+1:64-kmer.len+rem],rem), l, rem) != kmer
            if !([kmerge(lmer[1], DNASeq(kmer.bit1[64-kmer.len+1:64-kmer.len+rem],kmer.bit2[64-kmer.len+1:64-kmer.len+rem],rem), l, rem)] in keys(G))
                return Nothing
            end
            
            return [kmerge(lmer[1], DNASeq(kmer.bit1[64-k+2:64-l],kmer.bit2[64-kmer.len+1:64-l],rem), l, rem)]
        else
            return Nothing
        end
        """
        #if !([kmerge(lmer[1], DNASeq(kmer.bit1[64-kmer.len+1:64-kmer.len+rem],kmer.bit2[64-kmer.len+1:64-kmer.len+rem],rem), l, rem)] in keys(G))
        #    return Nothing
        #end
        
        if length(lmer) > 1
            return [kmer_seq(vcat(bit1[64-ll+1:64],lmer[2].bit1[1:rem]),vcat(bit2[64-ll+1:64],lmer[2].bit2[1:rem]),k-1)]
        else
            return [kmer_seq(vcat(bit1[64-ll+1:64],kmer.bit1[64-(k-1)+1:64-(k-1)+rem]), vcat(bit2[64-ll+1:64], kmer.bit2[64-(k-1)+1:64-(k-1)+rem]), k-1)]
        end
    end

end

function find_succ_ext(G::DefaultDict, kmer1 :: Vector{DNASeq}, kmer2 :: Vector{DNASeq})
    ### find prefix of the successor(kmer2) that is connected to the suffix of kmer 1
    k = kmer1[1].len +1
    TerminalEnd = true
    idx_1 = 0
    idx_2 = 0
    kmer2_ = copy(kmer2)
    while(TerminalEnd)
        if kmer2[1+idx_1] == Terminal
            deleteat!(kmer2_,1)
            idx_1 += 1
        elseif kmer2[end-idx_2] != Terminal
            TerminalEnd = false
        end

        if kmer2[end-idx_2] == Terminal
            deleteat!(kmer2_,length(kmer2_))
            idx_2 += 1
        elseif kmer2[1+idx_1] != Terminal
            TerminalEnd = false
        end
        
    end
    
    for (pid ,prefix) in G[kmer2_].prefixes
    
        if prefix[1] == Terminal
            continue
        elseif length(prefix)==1
             
            if prefix[1].len<(k-1) 
                if prefix[1].bit1[end-prefix[1].len+1:end] == kmer1[1].bit1[end-kmer1[1].len+1:end - kmer1[1].len + prefix[1].len] && prefix[1].bit2[end-prefix[1].len+1:end] == kmer1[1].bit2[end-kmer1[1].len+1:end - kmer1[1].len + prefix[1].len]
                    return pid, prefix
                end
            elseif prefix[1].bit1[end-prefix[1].len+1:end-prefix[1].len+k-1]==kmer1[1].bit1[end-(k-1)+1:end] && prefix[1].bit2[end-prefix[1].len + 1:end-prefix[1].len+k-1] == kmer1[1].bit2[end - (k-1)+1:end]
                return pid, prefix
            end
        elseif prefix[1].len > (k-1)
            if prefix[1].bit1[end-prefix[1].len+1:end-prefix[1].len+k-1]==kmer1[1].bit1[end-(k-1)+1:end] && prefix[1].bit2[end-prefix[1].len + 1:end-prefix[1].len+k-1] == kmer1[1].bit2[end - (k-1)+1:end]
                return pid,prefix
            end
        elseif prefix[1].bit1[end-prefix[1].len+1:end]==kmer1[1].bit1[end-(k-1)+1:end-(k-1)+prefix[1].len] && prefix[1].bit2[end-prefix[1].len + 1:end] == kmer1[1].bit2[end - (k-1)+ 1:end-(k-1)+prefix[1].len] 
            if prefix[2].bit1[1:k-1-prefix[1].len] == kmer1[1].bit1[end-(k-1)+prefix[1].len+1:end]  && prefix[2].bit2[1:k-1-prefix[1].len] == kmer1[1].bit2[end-(k-1)+prefix[1].len+1:end]
                return pid,prefix
            end
        end
    end
    return Nothing
end

function find_pred_ext(G :: DefaultDict, kmer1 ::Vector{DNASeq}, kmer2 :: Vector{DNASeq})
    ### find suffix of the predecessor(kmer1) that is connected to the prefix of kmer2
    k = kmer1[1].len +1
    TerminalEnd = true
    idx_1 = 0
    idx_2 = 0
    kmer1_ = copy(kmer1)
    while(TerminalEnd)
        if kmer1[1+idx_1] == Terminal
            deleteat!(kmer1_,1)
            idx_1 += 1
        elseif kmer1[end-idx_2] != Terminal
            TerminalEnd = false
        end

        if kmer1[end-idx_2] == Terminal
            deleteat!(kmer1_,length(kmer1_))
            idx_2 += 1
        elseif kmer1[1+idx_1] != Terminal
            TerminalEnd = false
        end
        
    end
    
    for (sid ,suffix) in G[kmer1_].suffixes
        #len = 0
        #len += (length(suffix)-1)*64 
        #if suffix[end] != Terminal
        #    len += suffix[end].len
        #end
        
        if length(suffix)==1 && suffix[1].len<(k-1) && suffix[1] != VTerminal
            if suffix[1].bit1[end-suffix[1].len+1:end] == kmer2[1].bit1[end - suffix[1].len+1:end] && suffix[1].bit2[end-suffix[1].len+1:end] == kmer2[1].bit2[end - suffix[1].len+1: end]
                return sid, suffix
            end
        elseif suffix[end].bit1[end-k+2:end]==kmer2[1].bit1[end-(k-1)+1:end] && suffix[end].bit2[end-k+2:end] == kmer2[1].bit2[end - (k-1)+1:end]
            return sid,suffix
        end
    end

    print("Nothing found\n")
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
        return [kmer_2]#, Terminal]
    elseif kmer_2 == Terminal
        return [kmer_1]#,Terminal]
            
    end
    if l+k <=64
        bit11 = kmer_1.bit1
        bit12 = kmer_1.bit2
        bit21 = kmer_2.bit1
        bit22 = kmer_2.bit2
        return [kmer_seq(vcat(bit11[end-k+1:end],bit21[end-l+1:end]), vcat(bit12[end-k+1:end],bit22[end-l+1:end]),k+l)]
    else
        bit11 = kmer_1.bit1
        bit12 = kmer_1.bit2
        bit21 = kmer_2.bit1
        bit22 = kmer_2.bit2
        rem = l+k-64
        #print(l,k,"\n",bit12[l+1:end],"\n",bit22[end-l+1])
        seq1 = kmer_seq(vcat(bit11[l+1:end],bit21[end-l+1:end]),vcat(bit12[l+1:end],bit22[end-l+1:end]),64)
        seq2 = kmer_seq(bit11[64-k+1:l],bit12[64-k+1:l],rem)
        return [seq2, seq1]
    end
end

function kmerge(kmer_1::Union{Vector{DNASeq},DNASeq}, kmer_2::Union{Vector{DNASeq},DNASeq})
    """if typeof(kmer_1) == Vector{DNASeq} 
        k = (length(kmer_1)-1)*64 + kmer_1[1].len
        if kmer_1[end] == Terminal && kmer_1[1] != Terminal
            print("Middle Terminal Error\n")
        end
    else
        k = kmer_1.len
    end
    if typeof(kmer_2) == Vector{DNASeq}
        
        l = (length(kmer_2)-1)*64 + kmer_2[1].len
        if kmer_2[end] == Terminal && length(kmer_2)>1
            print("A Terminal Issue\n")
            l -= 64
        end
    else
        l = kmer_2.len
    end"""
    if typeof(kmer_2) == Vector{DNASeq}
        return vcat(kmerge_(kmer_1,kmer_2[1]),kmer_2[2:end])
    else return kmerge_(kmer_1,kmer_2)
    end
end

function kmerge_(kmer_1::Vector{DNASeq}, kmer_2::DNASeq)
    l1 = (length(kmer_1)-1)*64+kmer_1[1].len
    l2 = kmer_2.len
    if (l1+l2)<=64
        bit11 = kmer_1[1].bit1[64-l1+1:end]
        bit12 = kmer_1[1].bit2[64-l1+1:end]
        bit21 = kmer_2.bit1[64-l2+1:end]
        bit22 = kmer_2.bit2[64-l2+1:end]
        return kmer_seq(vcat(bit11,bit21),vcat(bit12,bit22),l1+l2)
    else
        output_ = Vector{DNASeq}()
        cut_idx = 100
        idx = length(kmer_1)
        push!(output_, kmer_seq(vcat(kmer_1[end].bit1[l2+1:end],kmer_2.bit1[64-l2+1:end]),vcat(kmer_1[end].bit2[l2+1:end],kmer_2.bit2[64-l2+1:end]),64))
        while(idx > 2)
            pushfirst!(output_,kmer_seq(vcat(kmer_1[idx-1].bit1[l2+1:end],kmer_1[idx].bit1[1:l2]),vcat(kmer_1[idx-1].bit2[l2+1:end],kmer_1[idx].bit2[1:l2]),64))
            idx-=1;
        end
        
        if (64-kmer_1[1].len+1)<(l2+1)
            if length(kmer_1)> 1 
                pushfirst!(output_,kmer_seq(vcat(kmer_1[1].bit1[l2+1:end],kmer_1[2].bit1[1:l2]),vcat(kmer_1[1].bit2[l2+1:end],kmer_1[2].bit2[1:l2]),64))
            end
        pushfirst!(output_,kmer_seq(kmer_1[1].bit1[64-kmer_1[1].len+1:l2],kmer_1[1].bit2[64-kmer_1[1].len+1:l2],l2-64+kmer_1[1].len))
        else
            pushfirst!(output_,kmer_seq(vcat(kmer_1[1].bit1[64-kmer_1[1].len+1:end],kmer_1[2].bit1[1:l2]),vcat(kmer_1[1].bit2[64-kmer_1[1].len+1:end],kmer_1[2].bit2[1:l2]),l2+kmer_1[1].len))
        end
        return output_
    end



end

"""
function kmerge_(kmer_1::Vector{DNASeq}, kmer_2::DNASeq  )
    k = (length(kmer_1)-1)*64+kmer_1[1].len
    l = kmer_2.len
    k1 = k%64
    print(" k ",k,"k1 ",k1, " l ",l)
    res = []
    range = k>=(64-l) ? (l+1:64) : 64-k+1:64
    bit11 = kmer_1[end].bit1[range]
    bit12 = kmer_1[end].bit2[range]
    bit21 = kmer_2.bit1[64-l+1:end]
    bit22 = kmer_2.bit2[64-l+1:end]

    if k == k1
        print("here\n")
        pushfirst!(res, kmer_seq(vcat(bit11,bit21), vcat(bit12,bit22),length(range)+l))
    else
        print("her2\n")
        pushfirst!(res, kmer_seq(vcat(bit11,bit21), vcat(bit12,bit22),64))
    end

   for i in length(kmer_1):-1:3
        print("sher3\n")
        bit11 = kmer_1[i-1].bit1[l+1:end]
        bit12 = kmer_1[i-1].bit2[l+1:end]
        bit21 = kmer_1[i].bit1[1:l]
        bit22 = kmer_1[i].bit2[1:l]
        pushfirst!(res, kmer_seq(vcat(bit11,bit21), vcat(bit12,bit22),64))
    end

    
    if l+k1 > 64
        print("fok\n")
        bit11 = kmer_1[1].bit1[l+1:end]
        bit12 = kmer_1[1].bit2[l+1:end]
        bit21 = kmer_1[2].bit1[1:l]
        bit22 = kmer_1[2].bit2[1:l]
        pushfirst!(res, kmer_seq(vcat(bit11,bit21), vcat(bit12,bit22),64))

        bit11 = kmer_1[1].bit1[end-k1+1:l]
        bit12 = kmer_1[1].bit2[end-k1+1:l]
        pushfirst!(res,kmer_seq(bit11,bit12,l+k1-64))
    else
        print("chock\n")
        bit11 = kmer_1[1].bit1[end-k1+1:end]
        bit12 = kmer_1[1].bit2[end-k1+1:end]
        bit21 = kmer_1[2].bit1[1:l]
        bit22 = kmer_1[2].bit2[1:l]
        pushfirst!(res,kmer_seq(vcat(bit11,bit21),vcat(bit12,bit22),l+k1))
    end
    res = Vector{DNASeq}(res)
    return res
end
"""
function kmerge(kmer_1::DNASeq, kmer_2:: Vector{DNASeq})
    return vcat(kmerge(kmer_1,kmer_2[1],kmer_1.len+kmer_2[1].len),kmer_2[2:end])
end



function kmerge(kmer_1 :: Vector{DNASeq}, kmer_2 :: Vector{DNASeq}, k :: Int64, l :: Int64)
    l1_m = l%64
    k1_m = k%64
    l1_m = (l1_m == 0) ? 64 : l1_m
    k1_m = (k1_m == 0) ? 64 : k1_m
    res = Vector{DNASeq}()
    #print(l1_m, " ",k1_m ," ",k," ",l,"\n")
    if length(kmer_1) == 1 
        if length(kmer_2) == 1
            return kmerge(kmer_1[1],kmer_2[1],k,l)
        elseif k == 0
            return kmer_2
        elseif l1_m+k1_m <= 64
            return vcat([kmer_seq(vcat(kmer_1[1].bit1[end-k+1:end],kmer_2[1].bit1[end-l1_m+1:end]), vcat(kmer_1[1].bit2[end-k+1:end],kmer_2[1].bit2[end-l1_m+1 : end]),l1_m+k1_m)],kmer_2[2:end])
        else
            tempi = vcat([kmer_seq(vcat(kmer_1[1].bit1[l1_m+1:end],kmer_2[1].bit1[end-l1_m+1:end]), vcat(kmer_1[1].bit2[l1_m+1:end],kmer_2[1].bit2[end-l1_m+1 : end]),64)],kmer_2[2:end])
            tempo = l1_m+k1_m-64
            return vcat(kmer_seq(kmer_1[1].bit1[end-k+1:end-k+tempo],kmer_1[1].bit2[end-k+1:end-k+tempo],tempo),tempi)
        end
    elseif length(kmer_2)==1
        return kmerge(kmer_1,kmer_2[1])
    else
        res = kmer_2[2:end]
        pushfirst!(res,kmer_seq(vcat(kmer_1[end].bit1[l1_m+1 : end],kmer_2[1].bit1[end-l1_m+1:end]),vcat(kmer_1[end].bit2[l1_m+1 : end],kmer_2[1].bit2[end-l1_m+1:end]),64))
        
        for i in 1:(length(kmer_1)-2)
            pushfirst!(res,kmer_seq(vcat(kmer_1[end-i].bit1[l1_m+1:end],kmer_1[end-i+1].bit1[1:l1_m]),vcat(kmer_1[end-i].bit2[l1_m+1:end],kmer_1[end-i+1].bit2[1:l1_m]),64))
            
        end

        if k1_m+l1_m <=64
            pushfirst!(res,kmer_seq(vcat(kmer_1[1].bit1[end-k1_m+1:end],kmer_1[2].bit1[1:l1_m]),vcat(kmer_1[1].bit2[end-k1_m+1:end],kmer_1[2].bit2[1:l1_m]),k1_m+l1_m))
            
        else
            pushfirst!(res,kmer_seq(vcat(kmer_1[1].bit1[l1_m+1:end],kmer_1[2].bit1[1:l1_m]),vcat(kmer_1[1].bit2[l1_m+1:end],kmer_1[2].bit2[1:l1_m]),64))
            
            pushfirst!(res, kmer_seq(kmer_1[1].bit1[end-k1_m+1:l1_m],kmer_1[1].bit2[end-k1_m+1:l1_m],l1_m+k1_m-64))
        end
    end

    @assert(k != 0)
    return res

end

function IS(G::DefaultDict)#, Alphabet:: Vector{Char}, k :: Int64)
    IS_ = Set()
    pred_node = Vector{DNASeq}()
    succ_node = Vector{DNASeq}()
    key_notin_idset = false

    for node in keys(G)
        max_kmer = node
        for (pid, p_kmer) in G[node].prefixes
            if !G[node].prefixes_terminal[pid] #&& p_kmer[1].len != 0
                pred_node = pred_neigh(node[1], p_kmer)
                if !(pred_node in keys(G))
                    continue
                end
                if pred_node > node
                    max_kmer = pred_node
                    key_notin_idset = true
                    break
                end
                
            end

        end


        if !key_notin_idset
            for (sid, s_kmer) in G[node].suffixes
            
                if !G[node].suffixes_terminal[sid] #s_kmer[end].len != 0 && 
                    succ_node = succ_neigh(node[1], s_kmer)
                    if !(succ_node in keys(G))
                        continue
                    end
                    if succ_node == Nothing || !(succ_node in keys(G))
                        print("None\n")
                        print("Node ",node,"\nsucc_node ",succ_node,"\n s_kmer ",s_kmer,"\n\n")
                        continue
                    end

                    if succ_node > node
                        max_kmer = succ_node
                        key_notin_idset = true
                        break
                    end
                end

            end
        end

        if !key_notin_idset
            @assert(max_kmer == node)
            push!(IS_, node) 
        else
            @assert(max_kmer[1].len == node[1].len)
        
        end
        key_notin_idset = false
    end
    return IS_
end


function iterate_pack(G :: DefaultDict, IS_ :: Set)
    #print("iter call\n\n")
    #print(G[string_to_DNASeq("CTAAACTTCAACTCCCAAGCAACCAATTTAT")].wire_info,"\n")
    #print(G[string_to_DNASeq("ATCCGCCCTTTTTATTTAAGAAGCATAGAGG")].wire_info,"\n")
    transfer_nodeInfo = DefaultDict{Vector{DNASeq}, Set{Tuple}}(Set(()))
    pcontig_list = Vector{}()
    self_loop = [false, false]
    itr_p = 1
    itr_s = 1
    ppid = 0
    ssid = 0
    succ_ext = Vector{DNASeq}()
    pred_ext = Vector{DNASeq}()
    succ_node = Vector{DNASeq}()
    pred_node = Vector{DNASeq}()
    succ_ext_in = false
    pred_ext_in = false
    
    for node in IS_
        succ_ext_in = false
        pred_ext_in = false
        #print("IS node",DNASeq_to_string(node[1]),"\n\n")
        itr_p = 1
        itr_s = 1
        for k in 1:length(G[node].prefix_begin_info)

            if last(G[node].prefix_begin_info[k])>0
                itr_p = k
                if G[node].prefixes[itr_p][1].len != 0 && !G[node].prefixes_terminal[itr_p]
                    pred_node = pred_neigh(node[1],G[node].prefixes[itr_p])
                    if !(pred_node == Nothing) && (pred_node in keys(G))
                        ssid,pred_ext = find_pred_ext(G,pred_node, node)
                        pred_ext_in = true
                    end
                    #print(G[node].wire_info,"\n",G[node].prefix_begin_info,"\n",k)
                    for (itr_s,~,~) in G[node].wire_info[last(G[node].prefix_begin_info[k])]

                        if ! (G[node].prefixes_terminal[itr_p] && G[node].suffixes_terminal[itr_s])
                            if G[node].suffixes[itr_s][1].len>0 && !G[node].suffixes_terminal[itr_s]
                                succ_node = succ_neigh(node[1],G[node].suffixes[itr_s])
                                
                                #print(DNASeq_to_string(succ_node[1]),"\n",length(keys(G)),"\n")
                                #print(DNASeq_to_string(node[1]),"\n")
                                if !(succ_node == Nothing) && (succ_node in keys(G))
                                    ppid,succ_ext = find_succ_ext(G,node, succ_node)
                                    succ_ext_in = true
                                end
                                if pred_node == node 
                                    self_loop[1] = true
                                end
                                if succ_node == node
                                    self_loop[2] = true
                                end
                            end
                        end
                    end
                end

            end

        end
    end

    for node in IS_
        succ_ext_in = false
        pred_ext_in = false
        for  k in 1:length(G[node].prefix_begin_info)
            if last(G[node].prefix_begin_info[k])>0
                itr_p = k
                if length(G[node].prefixes[itr_p]) >0  && G[node].prefixes[itr_p][1].len != 0 && !G[node].prefixes_terminal[itr_p]
                    pred_node = pred_neigh(node[1],G[node].prefixes[itr_p])
                    if !(pred_node == Nothing) && (pred_node in keys(G))
                        ssid, pred_ext = find_pred_ext(G,pred_node, node)
                        pred_ext_in = true
                    end
                end
                
                for (itr_s,~,count) in G[node].wire_info[last(G[node].prefix_begin_info[k])]
                
                    if G[node].prefixes_terminal[itr_p] && G[node].suffixes_terminal[itr_s]
                        print("pcontig found\n",kmerge(kmerge(G[node].prefixes[itr_p],node),G[node].suffixes[itr_s]),"\n\n")
                        push!(pcontig_list,kmerge(kmerge(G[node].prefixes[itr_p],node),G[node].suffixes[itr_s]))
                    else
                        if length(G[node].suffixes[itr_s]) > 0 && !G[node].suffixes_terminal[itr_s] 
                            succ_node = succ_neigh(node[1],G[node].suffixes[itr_s])
                            if !(succ_node == Nothing) && (succ_node in keys(G))
                                ppid,succ_ext = find_succ_ext(G,node, succ_node)
                                succ_ext_in = true
                            end
                        end
                               

                        if !G[node].prefixes_terminal[itr_p] && !(pred_node == Nothing)
                            if pred_node != node
                                
                                if G[node].suffixes_terminal[itr_s]
                                    new_pnode_type = true
                                else
                                    if succ_node == node
                                        new_pnode_type = true
                                    else
                                        new_pnode_type = false
                                    end
                                end
                                if pred_ext!=VTerminal && pred_ext_in
                                    #print(pred_ext,"\n",G[node].suffixes[itr_s],"\n")
                                    new_ext = kmerge(pred_ext,G[node].suffixes[itr_s])
                                    if length(new_ext)>1 && new_ext[end-1] == Terminal
                                        new_ext = new_ext[1:end-1]
                                    elseif length(new_ext)>1 && new_ext[1] == Terminal
                                        new_ext = new_ext[2:end]
                                    end
                                    #print("pred node", DNASeq_to_string(pred_node[1])," node ", DNASeq_to_string(node[1])," ",new_ext,"\n")
                                    push!(transfer_nodeInfo[copy(node)] , (copy(pred_node),copy(pred_ext),new_ext,min(first(G[node].prefix_counts[itr_p]),first(G[node].suffix_counts[itr_s])),count,1,ssid,itr_s,new_pnode_type,pred_ext_in))

                                end
                            end
                        end
                        if !G[node].suffixes_terminal[itr_s] && !(succ_node == Nothing)
                            
                            if succ_node != node && succ_ext_in
                                if G[node].prefixes_terminal[itr_p]
                                    new_snode_type = true
                                else
                                    if pred_node == node
                                        new_snode_type = true
                                    else
                                        new_snode_type = false
                                    end
                                end
                                if succ_ext == VTerminal
                                    continue
                                end
                                
                                new_ext = kmerge(G[node].prefixes[itr_p],succ_ext)
                    
                                if length(new_ext)>1 && new_ext[2] == Terminal
                                    new_ext = new_ext[2:end]
                                end

                                #if succ_node == string_to_DNASeq("AAACTTCAACTCCCAAGCAACCAATTTATAT")
                                #    print("iterate info\n")
                                #    print(new_ext,"\n",G[string_to_DNASeq("AAACTTCAACTCCCAAGCAACCAATTTATAT")].prefixes[ppid],"\n",G[node].prefixes[itr_p],"\n\n")
                                #end
                                #print("succ node", DNASeq_to_string(succ_node[1])," node ", DNASeq_to_string(node[1])," ",new_ext,"\n")
                                push!(transfer_nodeInfo[copy(node)] , (copy(succ_node),copy(succ_ext),new_ext,min(first(G[node].suffix_counts[itr_s]),first(G[node].prefix_counts[itr_p])),count,0,ppid,itr_p,new_snode_type,succ_ext_in))
                                
                            end
                        end
                    end

                end
            end
        end

    end
    
    return pcontig_list, transfer_nodeInfo
end
   


function serialize_transfer!(G :: DefaultDict, transfer_nodeInfo :: DefaultDict{Vector{DNASeq}, Set{Tuple}, Set{Union{}}})
   
    rewirelist = []
    #print("transfer",transfer_nodeInfo,"\n\n\n\n")
    for node in keys(transfer_nodeInfo)
        for i in transfer_nodeInfo[node]
            (n_node,ext,new_ext,visit_count,ccount,direction,id,pid,type,ext_in) = i
            #if n_node == string_to_DNASeq("AAACTTCAACTCCCAAGCAACCAATTTATAT")
            #    print("serialize\n",new_ext,"\n")
            #end
            #direction == 1 ? print("whether\n",ext == G[n_node].suffixes[id]) : print("whether\n",ext == G[n_node].prefixes[id],"\n")
            
            push!(rewirelist,n_node)
            
            if ! (n_node in keys(G)) 
                print("MN node key was not found ")
                print(n_node,"\n")
            elseif direction == 1 
                if  !(G[n_node].suffixes[id] == ext) #|| !ext_in#!(new_ext in keys(G[n_node].suffixes)) && (ext in keys(G[n_node].suffixes))
                ##adds to suffix of next MN
                    if n_node == DNASeq[DNASeq(Bool[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 1, 1], Bool[0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1], 31)]
                        print("\n\nserialize10",ext,"\n",G[n_node].suffixes,"\n\n")
                    end
                    #print(ext_in,"\n")
                    #print("h1\n",DNASeq_to_string(n_node[1]),"\n",ext,"\n\n")
                    
                    l = length(G[n_node].suffixes) + 1
                    G[n_node].suffixes[l] = copy(new_ext) 
                    G[n_node].suffix_counts[l] = (visit_count,ccount)
                    G[n_node].suffixes_terminal[l] = type
                else
                    if G[n_node].suffixes_terminal[id]
                        l = length(G[n_node].suffixes) + 1
                        G[n_node].suffixes[l] = copy(new_ext) 
                        G[n_node].suffix_counts[l] = (visit_count,ccount)
                        G[n_node].suffixes_terminal[l] = type
                    else
                        G[n_node].suffixes[id] = copy(new_ext)
                        G[n_node].suffix_counts[id] = (visit_count,ccount)
                        G[n_node].suffixes_terminal[id] = type
                
                    end
                end
            else
                if !(G[n_node].prefixes[id] == ext) #|| !ext_in #!(new_ext in keys(G[n_node].prefixes)) && (ext in keys(G[n_node].suffixes))
                    ##adds to suffix of next MN
                        #print(ext_in,"\n")
                        #print("h2\n",DNASeq_to_string(n_node[1]),"\n",ext,"\n")
                        #print(new_ext,"\n")
                        
                        l = length(G[n_node].prefixes) + 1
                        G[n_node].prefixes[l] = copy(new_ext)
                        G[n_node].prefix_counts[l] = (visit_count,ccount)
                        G[n_node].prefixes_terminal[l] = type
    
                    else
                        if G[n_node].prefixes_terminal[id]
                            l = length(G[n_node].prefixes) + 1
                            G[n_node].prefixes[l] = copy(new_ext) 
                            G[n_node].prefix_counts[l] = (visit_count,ccount)
                            G[node].prefixes_terminal[l] = type
                    
                        else
                            G[n_node].prefixes[id] = copy(new_ext)
                            G[n_node].prefix_counts[id] = (visit_count,ccount)
                            G[n_node].prefixes_terminal[id] = type
                    
                        end
                        
    
                    end
            end
        end

    end
    return rewirelist
end



function initiate_wiring!(u:: macro_node)
    u.prefixes[length(u.prefixes)+1] = VTerminal
    u.suffixes[length(u.suffixes)+1] = VTerminal
    u.prefixes_terminal[length(u.prefixes)] = true
    u.suffixes_terminal[length(u.suffixes)] = true
    u.prefix_counts[length(u.prefixes)] = (-1,-1)
    u.suffix_counts[length(u.suffixes)] = (-1,-1)
end
