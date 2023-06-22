###WIP

include("graph_construction.jl")
include("kmer_counting_v2.jl")

function compact_graph( G :: DefaultDict{DNASeq,macro_node}, k:: Int64, phi :: Int64)
    num_mn = length(G)
    len = 1
    while(num_mn > phi)
        new_G, num_mn = IS(G,k,len,['A','C','G','T'])
        len *=2
    end
end

@inline function succ_neigh( kmer:: Union{Vector{DNASeq}, DNASeq} , lmer :: Vector{DNASeq}, k::Int64, l :: Int64)
    ## kmer is the key, lmer is the prefix
    bit1 = lmer[1].bit1
    bit2 = lmer[1].bit2
    if l >= (k-1)
        return kmer_seq(bit1[64-k+2:64],bit2[64-k+2:64],k-1)
    else
        rem = k-1-l
        return kmerge(kmer, lmer[1], rem, l)
    end
end
@inline function pred_neigh(kmer :: Union{Vector{DNASeq}, DNASeq} , lmer :: Vector{DNASeq} , k::Int64, l:: Int64)
    ## kmer is the key, lmer is the prefix
    bit1 = lmer[1].bit1
    bit2 = lmer[1].bit2
    if l>= (k-1)
        return kmer_seq(bit1[1:k-1],bit2[1:k-1],k-1)
    else
        rem = k-1-l
        return kmerge(lmer[1], kmer, l, rem)
    end

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
        return kmer_seq(vcat(bit11[64-k+1:64],bit21[64-l+1:64]), vcat(bit12[64-k+1:64],bit22[64-l+1:64]),k+l)
    else
        bit11 = kmer_1.bit1
        bit12 = kmer_1.bit2
        bit21 = kmer_2.bit1
        bit22 = kmer_2.bit2
        seq1 = kmer_seq(vcat(bit11[64-l+1:end],bit21),vcat(bit12[64-l+1:end],bit22),64)
        seq2 = kmer_seq(bit11[1:64-l],bit12[1:64-l],64-l)
        return [seq1, seq2]
    end


end
function kmerge(kmer_1::Vector{DNASeq}, kmer_2::DNASeq, k::Int64, l :: Int64 )
    k1 = k%64
    res = []

    bit11 = kmer_1[end].bit1[l+1:end]
    bit12 = kmer_1[end].bit2[l+1:end]
    bit21 = kmer_2.bit1
    bit22 = kmer_2.bit2
    pushfirst!(res, kmer_seq(vcat(bit11,bit21), vcat(bit12,bit22),64))
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

    bit11 = kmer_2[1].bit1[end-l1_m+1:end]
    bit12 = kmer_2[1].bit2[end-l1_m+1:end]
    bit21 = kmer_1[end].bit1[l1_m+1:end]
    bit22 = kmer_1[end].bit2[l1_m+1:end]
    pushfirst!(res,kmer_seq(vcat(bit21,bit11),vcat(bit22,bit12),64))
    
    
    lk1 = length(kmer_1)
    lk2 = length(kmer_2)
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


function IS(G::DefaultDict, Alphabet:: Vector{Char}, k :: Int64, l :: Int64)
    IS_ = Set()

    for node in keys(G)
        max_kmer = node
        for (pid, p_kmer) in G[node].prefixes
            if p_kmer[1] == Terminal
                continue
            end
            
            pred_node = pred_neigh(node, p_kmer, k,l)
            if pred_node>node
                max_kmer = pred_node
                break
            end
        end
        for (sid, s_kmer) in G[node].suffixes
            if s_kmer[1] == Terminal
                continue
            end
            succ_node = succ_neigh(node, s_kmer, k , l)
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


function iterate_pack(G :: DefaultDict{Int64,macro_node}, IS_ :: Set{DNASeq}, k:: Int64, l:: Int64)
    transfer_nodeInfo = DefaultDict{DNASeq,Set}(Set([]))
    pcontig_list = Vector{} 
    for node in IS_

        for (i,iset) in G[node].wire_info.prefix_info
            pid = i
            if G[node].prefix_terminal_id != pid && length(G[n].prefixes)>0
                pred_node = pred_neigh(node, G[node].prefixes[pid] , k, l)
                pred_ext = G[node].prefixes[pid]
            end

            for (j, offset, visit_count) in iset
                sid = j
                
                if G[node].prefix_terminal_id == pid && G[node].suffix_terminal_id == sid


                end
            end
            #else
            #    succ_node = succ_neigh(node, G[node].suffixes)
            #end
            for s in set
            end
        end

    end
end