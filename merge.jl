## merge and operator handling
include("datatypes.jl")

### Adjusting Base operators


### labels and macro_node

function Base.isequal(seq1 :: T, seq2 :: T) where T <: mn_label
    if seq1.bit1 == seq2.bit1 && seq1.bit2 == seq2.bit2 && seq1.isTerminal == seq2.isTerminal
        return true
    else
        return false
    end
end

Base.:(==)(seq1 :: T, seq2 :: T) where T <: mn_label = isequal(seq1, seq2)

function Base.isequal(m1 :: T, m2 :: T) where T <: macro_node
    return isequal(m1.label, m2.label)
end
Base.:(==)(m1 :: T, m2 :: T) where T<: macro_node = isequal(m1, m2)

### sequences


function Base.isequal(seq1::T, seq2::T) where T <: DNASeq
    if seq1.bit1 == seq2.bit1 && seq1.bit2 == seq2.bit2 && seq1.len == seq2.len
        return true
    else 
        return false
    end
end

function Base.isequal(seq1::Vector{T}, seq2::T) where T <: DNASeq
    if length(seq1)==1
        if seq1[1].bit1 == seq2.bit1 && seq1[1].bit2 == seq2.bit2 && seq1[1].len == seq2.len
         return true
        else 
            return false
        end
    else return false
    end
end

#function Base.isequal(seq1 :: t, seq2 :: VecT)

function Base.isequal(seq1::Vector{T}, seq2::Vector{T}) where T <: DNASeq
    if length(seq2)!= length(seq1)
        return false
    else
        for i in 1:length(seq1)
            if !isequal(seq1[i],seq2[i])
                return false
            end
        end
        return true
    end
end

Base.:(==)(seq1 ::Union{Vector{DNASeq},DNASeq}, seq2 :: Union{Vector{DNASeq},DNASeq}) = isequal(seq1, seq2)


function Base.isless(seq1 :: T, seq2 :: T) where T<: mn_label
    if seq1.isTerminal 
        return !seq2.isTerminal
    elseif seq2.isTerminal
        return false
    else
        if seq1.bit1 == seq2.bit1
            return isless(seq1.bit2, seq2.bit2)

        elseif isless(seq1.bit1, seq2.bit1)
            return true
        else
            return false
        end
    end
end

function Base.isless(seq1 :: T, seq2 :: T) where T<: DNASeq
    return seq1.bit1 == seq2.bit1 ? seq1.bit2<seq2.bit2 : seq1.bit1 < seq2.bit1
    #return [(i%2==1)*seq1.bit1[((i+1)÷2)+(i%2==0)-(i==64)]+seq1.bit2[(i÷2)+(i%2==1)]*(i%2==0) for i in 1:128]<[(i%2==1)*seq2.bit1[((i+1)÷2)+(i%2==0)-(i==128)] + seq2.bit2[(i÷ 2)+(i%2==1)]*(i%2==0) for i in 1:128]
end

function Base.isless(seq1 :: Vector{T}, seq2 :: Vector{T}) where T<: DNASeq

    for j in 1:length(seq1)
        temp1 = seq1[j]#[(i%2==1)*seq1[j].bit1[((i+1)÷2)+(i%2==0)-(i==128)]+seq1[j].bit2[(i÷2)+(i%2==1)]*(i%2==0) for i in 1:128]
        temp2 = seq2[j]#[(i%2==1)*seq2[j].bit1[((i+1)÷2)+(i%2==0)-(i==128)] + seq2[j].bit2[(i÷ 2)+(i%2==1)]*(i%2==0) for i in 1:128]
        if temp1== temp2
            continue
        else 
            return isless(temp1,temp2)
        end
    end
    return false
end

function Base.isless(seq1::T, seq2::Vector{T}) where T<:DNASeq
    return Base.isless(seq1,seq2[1])
end

function Base.hash(seq1::T) where T <: DNASeq
    return Base.hash(bitarr_to_int(vcat(seq1.bit1,seq1.bit2)))
end

function Base.hash(seq1::T) where T <: mn_label
    return Base.hash(bitarr_to_int(vcat(seq1.bit1,seq1.bit2)))
end


function Base.hash(seq1::Vector{T}) where T <: DNASeq
    if length(seq1)<=1
        return Base.hash(bitarr_to_int(vcat(seq1[1].bit1,seq1[1].bit2)))
    else
        return Base.hash(bitarr_to_int(vcat(seq1[1].bit1,seq1[1].bit2,seq1[2].bit1,seq1[2].bit1,seq1[2].bit2)))
    end
end

### edges

function Base.isequal(a :: edge, b :: edge)
    return Base.isequal(a.pred_label,b.pred_label) && Base.isequal(a.succ_label, b.succ_label) && Base.isequal(a.pred_suffix, b.pred_suffix) && Base.isequal(a.succ_prefix, b.succ_prefix)
end

function Base.hash(a :: edge)
    return Base.hash(a.id)
end


####

@inline function symmetric_bool(var::Bool) 
    2*Int(var) -1
end


##########
##########
## Merging





@inline function kmerge(c :: Char , kmer::DNASeq, order :: String)
    k = kmer.len
    if order == "from left"
        bit1 = zeros(Bool, k+1)
        bit2 = zeros(Bool, k+1)
        bit1[2:k+1] = kmer.bit1[65-k:64]
        bit2[2:k+1] = kmer.bit2[65-k:64]
        bit1[1], bit2[1] = char_to_int(c)
        lmer = DNASeq(bit1,bit2,k+1)
        
    else 
        bit1 = zeros(Bool, k+1)
        bit2 = zeros(Bool, k+1)
        bit1[1:k] = kmer.bit1[65-k:64]
        bit2[1:k] = kmer.bit2[65-k:64]
        bit1[k+1], bit2[k+1] = char_to_int(c)
        lmer = DNASeq(bit1,bit2,k+1)
    end

    lmer
end

@inline function kmerge(c :: Char , kmer::mn_label, order :: String)
    kmer = DNASeq(kmer.bit1, kmer.bit2, 31)
    k = kmer.len
    if order == "from left"
        bit1 = zeros(Bool, k+1)
        bit2 = zeros(Bool, k+1)
        bit1[2:k+1] = kmer.bit1[65-k:64]
        bit2[2:k+1] = kmer.bit2[65-k:64]
        bit1[1], bit2[1] = char_to_int(c)
        lmer = DNASeq(bit1,bit2,k+1)
        
    else 
        bit1 = zeros(Bool, k+1)
        bit2 = zeros(Bool, k+1)
        bit1[1:k] = kmer.bit1[65-k:64]
        bit2[1:k] = kmer.bit2[65-k:64]
        bit1[k+1], bit2[k+1] = char_to_int(c)
        lmer = DNASeq(bit1,bit2,k+1)
    end

    lmer
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
        return [DNASeq(vcat(bit11[end-k+1:end],bit21[end-l+1:end]), vcat(bit12[end-k+1:end],bit22[end-l+1:end]),k+l)]
    else
        bit11 = kmer_1.bit1
        bit12 = kmer_1.bit2
        bit21 = kmer_2.bit1
        bit22 = kmer_2.bit2
        rem = l+k-64
        #print(l,k,"\n",bit12[l+1:end],"\n",bit22[end-l+1])
        seq1 = DNASeq(vcat(bit11[l+1:end],bit21[end-l+1:end]),vcat(bit12[l+1:end],bit22[end-l+1:end]),64)
        seq2 = DNASeq(bit11[64-k+1:l],bit12[64-k+1:l],rem)
        return [seq2, seq1]
    end
end

function kmerge(kmer_1::Union{Vector{DNASeq},DNASeq}, kmer_2::Union{Vector{DNASeq},DNASeq})
    if typeof(kmer_1) == Vector{DNASeq} 
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
        pushfirst!(res, DNASeq(vcat(bit11,bit21), vcat(bit12,bit22),length(range)+l))
    else
        pushfirst!(res, DNASeq(vcat(bit11,bit21), vcat(bit12,bit22),64))
    end

   for i in length(kmer_1):-1:3
        bit11 = kmer_1[i-1].bit1[l+1:end]
        bit12 = kmer_1[i-1].bit2[l+1:end]
        bit21 = kmer_1[i].bit1[1:l]
        bit22 = kmer_1[i].bit2[1:l]
        pushfirst!(res, DNASeq(vcat(bit11,bit21), vcat(bit12,bit22),64))
    end

    
    if l+k1 > 64
        bit11 = kmer_1[1].bit1[l+1:end]
        bit12 = kmer_1[1].bit2[l+1:end]
        bit21 = kmer_1[2].bit1[1:l]
        bit22 = kmer_1[2].bit2[1:l]
        pushfirst!(res, DNASeq(vcat(bit11,bit21), vcat(bit12,bit22),64))

        bit11 = kmer_1[1].bit1[end-k1+1:l]
        bit12 = kmer_1[1].bit2[end-k1+1:l]
        pushfirst!(res,DNASeq(bit11,bit12,l+k1-64))
    else
        bit11 = kmer_1[1].bit1[end-k1+1:end]
        bit12 = kmer_1[1].bit2[end-k1+1:end]
        bit21 = kmer_1[2].bit1[1:l]
        bit22 = kmer_1[2].bit2[1:l]
        pushfirst!(res,DNASeq(vcat(bit11,bit21),vcat(bit12,bit22),l+k1))
    end
    res = Vector{DNASeq}(res)
    return res
end

function kmerge(kmer_1::DNASeq, kmer_2:: Vector{DNASeq}, k::Int64, l :: Int64 )
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
            return vcat([DNASeq(vcat(kmer_1[1].bit1[end-k+1:end],kmer_2[1].bit1[end-l1_m+1:end]), vcat(kmer_1[1].bit2[end-k+1:end],kmer_2[1].bit2[end-l1_m+1 : end]),l1_m+k1_m)],kmer_2[2:end])
        else
            tempi = vcat([DNASeq(vcat(kmer_1[1].bit1[l1_m+1:end],kmer_2[1].bit1[end-l1_m+1:end]), vcat(kmer_1[1].bit2[l1_m+1:end],kmer_2[1].bit2[end-l1_m+1 : end]),64)],kmer_2[2:end])
            tempo = l1_m+k1_m-64
            return vcat(DNASeq(kmer_1[1].bit1[end-k+1:end-k+tempo],kmer_1[1].bit2[end-k+1:end-k+tempo],tempo),tempi)
        end
    elseif length(kmer_2)==1
        return kmerge(kmer_1,kmer_2[1])
    else
        res = kmer_2[2:end]
        pushfirst!(res,DNASeq(vcat(kmer_1[end].bit1[l1_m+1 : end],kmer_2[1].bit1[end-l1_m+1:end]),vcat(kmer_1[end].bit2[l1_m+1 : end],kmer_2[1].bit2[end-l1_m+1:end]),64))
        
        for i in 1:(length(kmer_1)-2)
            pushfirst!(res,DNASeq(vcat(kmer_1[end-i].bit1[l1_m+1:end],kmer_1[end-i+1].bit1[1:l1_m]),vcat(kmer_1[end-i].bit2[l1_m+1:end],kmer_1[end-i+1].bit2[1:l1_m]),64))
            
        end

        if k1_m+l1_m <=64
            pushfirst!(res,DNASeq(vcat(kmer_1[1].bit1[end-k1_m+1:end],kmer_1[2].bit1[1:l1_m]),vcat(kmer_1[1].bit2[end-k1_m+1:end],kmer_1[2].bit2[1:l1_m]),k1_m+l1_m))
            
        else
            pushfirst!(res,DNASeq(vcat(kmer_1[1].bit1[l1_m+1:end],kmer_1[2].bit1[1:l1_m]),vcat(kmer_1[1].bit2[l1_m+1:end],kmer_1[2].bit2[1:l1_m]),64))
            
            pushfirst!(res, DNASeq(kmer_1[1].bit1[end-k1_m+1:l1_m],kmer_1[1].bit2[end-k1_m+1:l1_m],l1_m+k1_m-64))
        end
    end

    @assert(k != 0)
    return res

end

