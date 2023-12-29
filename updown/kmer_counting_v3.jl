using DataStructures, FastaIO

struct DNASeq
    bit1 :: BitArray{1}
    bit2 :: BitArray{1}
    len :: Int32
end

Base.@kwdef mutable struct mn_label
    bit1 = BitArray(undef, 31)
    bit2 = BitArray(undef, 31)
end

### to use macro_nodes in dict:


function Base.isequal(a :: mn_label, b :: mn_label)
    return Base.isequal(a.bit1,b.bit1) && Base.isequal(a.bit2, b.bit2)
end

function Base.isless(seq1 :: T, seq2 :: T) where T<: mn_label
    return seq1.bit1 == seq2.bit1 ? seq1.bit2<seq2.bit2 : seq1.bit1 < seq2.bit1
end


function Base.hash(seq :: T) where T <: mn_label
    return Base.hash(bitarr_to_int(vcat(seq.bit1,seq.bit2)))
end




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
    if length(seq2)== 0
        return false
    end
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

function bitarr_to_int(arr,s=Int64(0))
    v = 1
    for i in view(arr,length(arr):-1:1)
        s += v*i
        v <<= 1
    end 
    s
end

@inline function char_to_int( c :: Union{Char, String})
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
    temp = char_to_int.(Vector{Char}(seq))
    hcat(first.(temp), last.(temp))
end

function string_to_DNASeq(seq :: String)
    len = length(seq)
    bitseq = Vector{DNASeq}()
    idx = 1
    #print(floor(Int64,len/32),"\n")
    while(floor(Int64,len/32)>0)
        temp = string_to_int(seq[end-(32*idx)+1:end-(32*(idx-1))])
        pushfirst!(bitseq,DNASeq(BitArray(vec(temp[:,1])),BitArray(vec(temp[:,2])),32))
        idx += 1
        len -= 32
    end

    if len%32 != 0
        temp_1 = zeros(Bool,32)
        temp_2 = zeros(Bool,32)
        temp = string_to_int(seq[1:(len%32)])
        temp_1[32-length(temp[:,1])+1:32] = temp[:,1]
        temp_2[32-length(temp[:,2])+1:32] = temp[:,2]
        pushfirst!(bitseq,DNASeq(BitArray(vec(temp_1)),BitArray(vec(temp_2)),length(temp[:,1])))
    end
    bitseq

end

function DNASeq_to_string(seq:: DNASeq)
    k = seq.len
    bit1 = seq.bit1
    bit2 = seq.bit2
    str =""
    for i in 32-k+1:32
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


function kmer_seq( bit1 :: Union{BitArray,SubArray}, bit2 :: Union{BitArray,SubArray}, k :: Int64)
    Bit1 = BitArray(undef,32)
    Bit2 = BitArray(undef,32)
    Bit1[end-k+1 : end] = bit1
    Bit2[end-k+1 : end] = bit2
    output = DNASeq(Bit1, Bit2,k)
    output
end


@inline function read_kmer(seq :: Vector, k :: Int64)
    @assert(k<33)
    kmer_list = DefaultDict{DNASeq, Int64}(0)
    temp = copy(seq)
    while(length(temp)>2)
        for i in 1:33-k
            bit1 = @view temp[end].bit1[i:i+k-1]
            bit2 = @view temp[end].bit2[i:i+k-1]
            kmer_list[kmer_seq(bit1,bit2,k)] += 1
        end
        for i in 1:k-1
            kmer_list[kmer_seq(vcat(temp[end-1].bit1[end-k+i+1:end],temp[end].bit1[1:i]), vcat(temp[end-1].bit2[end-k+i+1:end],temp[end].bit2[1:i]),k)] += 1
        end
        pop!(temp)
    end

    if length(temp)==2
        for i in 1:33-k
            bit1 = @view temp[end].bit1[i:i+k-1]
            bit2 = @view temp[end].bit2[i:i+k-1]
            kmer_list[kmer_seq(bit1,bit2,k)] += 1
        end
        for i in k-1:-1:1
            if temp[1].len<k-i
                return kmer_list
            else
                kmer_list[kmer_seq(vcat(temp[end-1].bit1[end-k+i+1:end],temp[end].bit1[1:i]), vcat( temp[end-1].bit2[end-k+i+1:end],temp[end].bit2[1:i]),k)] += 1
            end
        end
        pop!(temp)
    end
    
    if length(temp)==1 && temp[1].len>=k
        for i in 32-temp[1].len+1:32-k+1
            bit1 = @view temp[1].bit1[i:i+k-1]
            bit2 = @view temp[1].bit2[i:i+k-1]
            kmer_list[kmer_seq(bit1,bit2,k)] += 1
        end
    end
    return kmer_list
    print("Read Error\n")
end


### to print kmers in the output:
function seq_to_int64(seq,s=Int64(0))
    v = 1
    #for i in view(arr,length(arr):-1:1)
    j = 1
    for i in 64:-1:64-2*seq.len+1
        if i%2 == 0
            bit = seq.bit2[Int64(i/2)]
            s += v*bit
            v <<= 1
        else
            bit = seq.bit1[Int64((i+1)/2)]
            s += v*bit
            v <<= 1
        end

    end 
    s
end

function read_from_kmer(kmer :: DNASeq)
    k = kmer.len
    return mn_label(kmer.bit1[1:31],kmer.bit2[1:31]), mn_label(kmer.bit1[2:32],kmer.bit2[2:32])

end

function save_output(file :: String, k :: Int64, wfile :: String = Nothing)
    kmer_list =  DefaultDict{DNASeq, Int64}(0)
    for (~,seq) in FastaReader(file)
        input = string_to_DNASeq(seq)
        #print(seq,"\n")
        temp_kmers = read_kmer(input,k)
        for (t,v) in temp_kmers
            kmer_list[t] += v
        end
    end
    print("Finished Reading kmers\n")
    output = ""
    for(i,j) in kmer_list
        #print(seq_to_int64(i),"\t",j,"\n")
        output = output*"$(seq_to_int64(i))\t$j\n"
    end
    print("Finished Creating Output String")
    if wfile!= Nothing
        open(wfile, "w") do file_
            write(file_, output)
        end
    else
        for(i,j) in kmer_list
            print(seq_to_int64(i),"\t",j,"\n")
            #output = output*"$(seq_to_int64(i))\t$j\n"
        end

    end
end


function save_kmer_dist(file :: String, k :: Int64, wfile :: String = Nothing)
    kmer_list =  DefaultDict{DNASeq, Int64}(0)
    idx = 0
    print("Started Reading")
    for (~,seq) in FastaReader(file)
        idx += 1
        input = string_to_DNASeq(seq)
        #print(seq,"\n")
        temp_kmers = read_kmer(input,k)
        for (t,v) in temp_kmers
            kmer_list[t] += v
        end
        if idx %100000 == 0
            print("reading\n")
        end
    end
    print("Finished Reading kmers\n")
    
    temp = sort(collect(values(kmer_list)))
    len = temp[end]
 
    arr = zeros(Int64, len)
    print(len,"\n")
    for i in temp
        arr[i] += 1
    end
    
    print("Finished Creating Output Array")
    if wfile!= Nothing
        open(wfile, "w") do file_
            for i in eachindex(arr)
                write(file_,"$i\t $(arr[i])\n")
            end
        end
    end

end
