
include("merge.jl")

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
    m = length(seq)%64
    if m == 0
        m = 64
    end
    
    bitseq = Vector{DNASeq}()
    for i in 1:n-1
        temp_1 = zeros(Int64,64)
        temp_2 = zeros(Int64,64)
        temp = string_to_int(seq[64*(i-1)+m+1:64*i+m])
        temp_1[65-length(temp[:,1]):64] = temp[:,1]
        temp_2[65-length(temp[:,2]):64] = temp[:,2]
        #print(temp_1)
        push!(bitseq,DNASeq(vec(temp_1),vec(temp_2),64))
    end
    #print(",",m,"\n")
    temp_1 = zeros(Int64,64)
    temp_2 = zeros(Int64,64)
    temp = string_to_int(seq[1:m])
    temp_1[65-length(temp[:,1]):64] = temp[:,1]
    temp_2[65-length(temp[:,2]):64] = temp[:,2]
    pushfirst!(bitseq,DNASeq(vec(temp_1),vec(temp_2),m))
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



function read_kmer(seq:: Vector, len::Int64, k :: Int64)
    ## k is the kmer length
    kmer_list = DefaultDict{DNASeq, Int64}(0)
    max_ind = length(seq)
    len_max_ind = len%64
    if len_max_ind == 0
        len_max_ind = 64
    end
    
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
                    #kmer = kmer_seq(seq[ind1].bit1[64-len_max_ind+offset-k+1:i], seq[ind1].bit2[64-len_max_ind+offset-k+1:64-len_max_ind+offset],k)
                
                    kmer = DNASeq(seq[ind1].bit1[64-len_max_ind+offset-k+1:64-len_max_ind+offset], seq[ind1].bit2[64-len_max_ind+offset-k+1:64-len_max_ind+offset],k)
                else
                    kmer = DNASeq(seq[ind1].bit1[offset-k+1:offset], seq[ind1].bit2[offset-k+1:offset],k)
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
                kmer = DNASeq(bit1,bit2,k)
                kmer_list[kmer] += 1

               
            end
        end
    else
        print("error")
    end
    kmer_list
end


@inline function read_mn_from_kmer(kmer :: DNASeq)
    mn_list = DefaultDict{mn_label, Int64}(0)
    for i in 64-k+1:64-k+2
        lmer = mn_label(kmer.bit1[i:i+k-2],kmer.bit2[i:i+k-2],false)
        mn_list[lmer] += 1
       
    end

    mn_list
end

@inline function neigh_label( c :: Char, u :: Union{mn_label, macro_node}, order :: String )
    if typeof(u) == macro_node
        u = u.label
    end
    if order == "pred"
        if c =='A'
            return mn_label(vcat(false, u.bit1[1:end-1]), vcat(false, u.bit1[1:end-1]), false)
        elseif c == 'C'
            return mn_label(vcat(false, u.bit1[1:end-1]), vcat(true, u.bit1[1:end-1]), false)
        elseif c == 'G' 
            return mn_label(vcat(true, u.bit1[1:end-1]), vcat(false, u.bit1[1:end-1]), false)
        else
            return mn_label(vcat(true, u.bit1[1:end-1]), vcat(true, u.bit1[1:end-1]), false)
        end
    elseif order == "succ"
        if c =='A'
            return mn_label(vcat(u.bit1[1:end-1], false), vcat(u.bit1[1:end-1], false), false)
        elseif c == 'C'
            return mn_label(vcat(u.bit1[1:end-1], false), vcat(u.bit1[1:end-1], true), false)
        elseif c == 'G' 
            return mn_label(vcat(u.bit1[1:end-1], true), vcat(u.bit1[1:end-1], false), false)
        else
            return mn_label(vcat(u.bit1[1:end-1], true), vcat(u.bit1[1:end-1], true), false)
        end
    end
end