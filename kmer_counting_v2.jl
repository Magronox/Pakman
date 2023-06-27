using DataStructures, Random



struct DNASeq
    bit1::BitArray{1}
    bit2::BitArray{1}
    len :: Int64
end

import Base.isequal, Base.hash


function Base.isequal(seq1::T, seq2::T) where T <: DNASeq
    if seq1.bit1 == seq2.bit1 && seq1.bit2 == seq2.bit2
        return true
    else 
        return false
    end
end

function Base.isless(seq1 :: T, seq2 :: T) where T<: DNASeq
    return [(i%2==1)*seq1.bit1[((i+1)÷2)+(i%2==0)-(i==128)]+seq1.bit2[(i÷2)+(i%2==1)]*(i%2==0) for i in 1:128]<[(i%2==1)*seq2.bit1[((i+1)÷2)+(i%2==0)-(i==128)] + seq2.bit2[(i÷ 2)+(i%2==1)]*(i%2==0) for i in 1:128]
end
function Base.hash(seq1::T) where T <: DNASeq
    return Base.hash(bitarr_to_int(vcat(seq1.bit1,seq1.bit2)))
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
    if     only(c) == 'A' return 0,0
    elseif only(c) == 'C' return 0,1
    elseif only(c) == 'G' return 1,0
    elseif only(c) == 'T' return 1,1
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
    bitseq = []
    for i in 1:n-1
        temp_1 = zeros(Int64,64)
        temp_2 = zeros(Int64,64)
        temp = string_to_int(seq[64*(i-1)+1:64*i])
        temp_1[65-length(temp[:,1]):64] = temp[:,1]
        temp_2[65-length(temp[:,2]):64] = temp[:,2]
        #print(temp_1)
        push!(bitseq,DNASeq(BitArray(vec(temp_1)),BitArray(vec(temp_2)),64))
    end
    temp_1 = zeros(Int64,64)
    temp_2 = zeros(Int64,64)
    temp = string_to_int(seq[64*(n-1)+1:end])
    temp_1[65-length(temp[:,1]):64] = temp[:,1]
    temp_2[65-length(temp[:,2]):64] = temp[:,2]
    push!(bitseq,DNASeq(BitArray(vec(temp_1)),BitArray(vec(temp_2)),n%64))
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
            str = str*"T"
        else
            str = str*"G"
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

function read_kmer(seq:: Vector, len::Int64, k :: Int64)
    ## k is the kmer length
    kmer_list = DefaultDict{DNASeq, Int64}(0)
    max_ind = length(seq)
    len_max_ind = len%64
    if (len>=k)
        for i in k:len

            ind1 = floor(Int64,(i-k+1)/64)
            if (i-k+1)%64 !=0 
                ind1 += 1
            end

            ind2 = floor(Int64, i/64)
            if (i%64) != 0
                ind2 += 1
            end
            
            if ind1==ind2
               
                offset = i%64
                if offset==0
                    offset = 64
                end
                
                
                if ind1==max_ind
                    kmer = kmer_seq(seq[ind1].bit1[64-len_max_ind+offset-k+1:64-len_max_ind+offset], seq[ind1].bit2[64-len_max_ind+offset-k+1:64-len_max_ind+offset],k)
                else
                    kmer = kmer_seq(seq[ind1].bit1[offset-k+1:offset], seq[ind1].bit2[offset-k+1:offset],k)
                end
                
                kmer_list[kmer] += 1

            else
               
                ind1 = floor(Int64,(i-k+1)/64)
                offset_1 = (i-k+1)%64
                if offset_1 == 0
                     offset_1 = 64
                else
                    ind1 += 1
                end
                ind2 = floor(Int64, i/64)
                offset_2 = i%64
                if offset_2 == 0
                    offset_2 = 64
                else 
                    ind2 += 1
                end

                chunk11 = seq[ind1].bit1[offset_1:end]
                chunk12 = seq[ind1].bit2[offset_1:end]
                chunk21 = seq[ind2].bit1[1:offset_2]
                chunk22 = seq[ind2].bit2[1:offset_2]
                #print(chunk11,"\n",chunk21,"\n",i," ",ind1,ind2)
                bit1 = vcat(chunk11,chunk21)
                bit2 = vcat(chunk12,chunk22)
                kmer = kmer_seq(bit1,bit2,k)
                kmer_list[kmer] += 1
               
            end
        end
    else
        print("error")
    end
    kmer_list
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


