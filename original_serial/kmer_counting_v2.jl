using DataStructures, Random



struct DNASeq
    bit1 :: BitArray{1}
    bit2 :: BitArray{1}
    len :: Int64
end

#struct DNAVseq
#    bit1 :: Vector{BitArray{1}}
#    bit2 :: Vector{BitArray{1}}
#    len :: Int64
#end



import Base.isequal, Base.hash


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
function Base.isequal(seq1::T, seq2::Vector{T}) where T <: DNASeq
    if length(seq2)==1
        if seq1.bit1 == seq2[1].bit1 && seq1.bit2 == seq2[1].bit2 && seq1.len == seq2[1].len
         return true
        else 
            return false
        end
    else return false
    end
end
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


Base.:(==)(seq1 ::Union{Vector{DNASeq},DNASeq}, seq2 :: Union{Vector{DNASeq},DNASeq}) = isequal(seq1,seq2)


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
function Base.hash(seq1::Vector{T}) where T <: DNASeq
    if length(seq1)<=1
        return Base.hash(bitarr_to_int(vcat(seq1[1].bit1,seq1[1].bit2)))
    else
        return Base.hash(bitarr_to_int(vcat(seq1[1].bit1,seq1[1].bit2,seq1[2].bit1,seq1[2].bit1,seq1[2].bit2)))
    end
end

function bitarr_to_int(arr,s=Int128(0))
    v = 1
    for i in view(arr,length(arr):-1:1)
        s += v*i
        v <<= 1
    end 
    s
end


function char_to_int( c :: Union{Char, String})
    #@assert(sizeof(c)==1)
    if     only(c) == 'A' return false, false
    elseif only(c) == 'C' return false, true
    elseif only(c) == 'G' return true , true
    elseif only(c) == 'T' return true , false
    else return 42,42
    end
end

function DNASeq_to_int(DNAseq :: DNASeq)
    return bitarr_to_int(vcat(DNAseq.bit1,DNAseq.bit2))
end

function string_to_int( seq :: String)
    temp = char_to_int.(only.(split(seq,"")))
    output = hcat(first.(temp), last.(temp))
    output
end


function string_to_DNASeq(seq :: String)
   
    n = ceil(Int64, length(seq)/64)
    m = length(seq)%64
    if m == 0
        m = 64
    end
    
    bitseq = Vector{DNASeq}()
    for i in 1:n-1
        temp_1 = zeros(Bool,64)
        temp_2 = zeros(Bool,64)
        temp = string_to_int(seq[64*(i-1)+m+1:64*i+m])
        temp_1[65-length(temp[:,1]):64] = temp[:,1]
        temp_2[65-length(temp[:,2]):64] = temp[:,2]
        #print(temp_1)
        push!(bitseq,DNASeq(BitArray(vec(temp_1)),BitArray(vec(temp_2)),64))
    end
    #print(",",m,"\n")
    temp_1 = zeros(Bool,64)
    temp_2 = zeros(Bool,64)
    temp = string_to_int(seq[1:m])
    temp_1[65-length(temp[:,1]):64] = temp[:,1]
    temp_2[65-length(temp[:,2]):64] = temp[:,2]
    pushfirst!(bitseq,DNASeq(BitArray(vec(temp_1)),BitArray(vec(temp_2)),m))
    bitseq

end

 function DNASeq_to_string(seq:: DNASeq)
    k = seq.len
    bit1 = seq.bit1
    bit2 = seq.bit2
    str =""
    for i in 64-k+1:64
        if bit1[i]==0 && bit2[i] ==0
            str = str* "A"
        elseif bit1[i]==0 && bit2[i]==1
            str = str* "C"
        elseif bit1[i]==1 && bit2[i]==1
            str = str*"G"
        else
            str = str*"T"
        end
    end
    str
end


function kmer_seq( bit1 :: BitArray, bit2 :: BitArray, k :: Int64)
    # @assert(length(bit1)==length(bit2))
    # @assert(length(bit1)==k)
    # @assert(k<=64)
    Bit1 = BitArray(undef,64)
    Bit2 = BitArray(undef,64)
    Bit1[end-k+1 : end] = bit1
    Bit2[end-k+1 : end] = bit2

    output = DNASeq(Bit1, Bit2,k)
    output
end

function read_kmer(seq :: Vector, len :: Int64, k :: Int64)
    @assert(k<64)
    kmer_list = DefaultDict{DNASeq, Int64}(0)
    if length(seq) == 1
        if len<k
            print("Error\n")
        else
            for i in 64-len+1:64-k+1
                #print(i,"\n")
                bit1 = seq[1].bit1[i:i+k-1]
                bit2 = seq[1].bit2[i:i+k-1]
                kmer_list[kmer_seq(bit1,bit2,k)] += 1
            end
            return kmer_list
        end
    else
        len1 = seq[1].len
        for i in 64-len1+1:64
            ind2 = ((i+k-1)>64) ? 2 : 1
            if ind2 == 1
                bit1 = seq[1].bit1[i:i+k-1]
                bit2 = seq[1].bit2[i:i+k-1]
                kmer_list[kmer_seq(bit1,bit2,k)] += 1
            else
                bit1 = vcat(seq[1].bit1[i:end],seq[2].bit1[1:k+i-65])
                bit2 = vcat(seq[1].bit2[i:end],seq[2].bit2[1:k+i-65])
                kmer_list[kmer_seq(bit1,bit2,k)] += 1
            end
        end
        for idx in 2:length(seq)
            for i in 1:64
                ind2 = ((i+k-1)>64) ? idx+1 : idx
                if ind2 == idx
                    bit1 = seq[idx].bit1[i:i+k-1]
                    bit2 = seq[idx].bit2[i:i+k-1]
                    kmer_list[kmer_seq(bit1,bit2,k)] += 1
                elseif (idx == length(seq))
                    break
                else
                    bit1 = vcat(seq[idx].bit1[i:end],seq[ind2].bit1[1:k+i-65])
                    bit2 = vcat(seq[idx].bit2[i:end],seq[ind2].bit2[1:k+i-65])
                    kmer_list[kmer_seq(bit1,bit2,k)] += 1
                end
            end
        end
        return kmer_list
    end
    print("Read Error\n")
end


function read_lmer_from_kmer(kmer :: DNASeq, l :: Int64)
    k = kmer.len 
    @assert(l<k)
    ## l is length of the kmers we are searching for
    bit1 = kmer.bit1
    bit2 = kmer.bit2
    lmer_list = DefaultDict{DNASeq, Int64}(0)
   
   
    for i in 64-k+1:64-l+1
        lmer = kmer_seq(bit1[i:i+l-1],bit2[i:i+l-1],l)
        lmer_list[lmer] += 1
    end

    lmer_list
end



## Test:
## slen = 1000
## input = randstring("ACGT",slen)

## DNA_seq = string_to_DNASeq(input)
## k=3
## kmer_list = read_kmer(DNA_seq, length(input),k)


