include("kmer_counting_v2.jl")






Base.@kwdef mutable struct macro_node
    #length k-1
    #id = 0
    label = DNASeq#BitArray{1}(undef,64)
    prefixes = DefaultDict{Int64, Vector{DNASeq}}(DNASeq[])
    prefix_counts = Dict{Int64, Tuple}()
    prefixes_terminal = Vector{Bool}()
    prefix_terminal = false
    suffixes = DefaultDict{Int64, Vector{DNASeq}}(DNASeq[])
    suffix_counts = Dict{Int64,Tuple}()
    suffix_terminal = false
    suffixes_terminal = Vector{Bool}()
    wire_info = DefaultDict{Int64, Set{Tuple}}(Set())
    prefix_begin_info = DefaultDict{Int64,Tuple}((-1,-1))
end


### to use macro_nodes in dict:

function Base.isequal(a :: macro_node, b :: macro_node)
    return Base.isequal(a.label,b.label)
end
function Base.hash(a :: macro_node)
    return Base.hash(a.id)
end


## define terminal 

Terminal = DNASeq(BitArray{1}(undef,64),BitArray{1}(undef,64),0)

VTerminal = [Terminal]

### to get the maximum values in dictionary
Base.isless(p::Pair, q::Pair) =
               isless(p.second,q.second)
           

@inline function kmerge(c :: Char , kmer::DNASeq, reverse::Bool)
    k = kmer.len
    if reverse == false
        bit1 = BitArray{1}(undef, k+1)
        bit2 = BitArray{1}(undef, k+1)
        bit1[2:k+1] = kmer.bit1[65-k:64]
        bit2[2:k+1] = kmer.bit2[65-k:64]
        bit1[1], bit2[1] = char_to_int(c)
        lmer = kmer_seq(bit1,bit2,k+1)
        
    else 
        bit1 = BitArray{1}(undef, k+1)
        bit2 = BitArray{1}(undef, k+1)
        bit1[1:k] = kmer.bit1[65-k:64]
        bit2[1:k] = kmer.bit2[65-k:64]
        bit1[k+1], bit2[k+1] = char_to_int(c)
        lmer = kmer_seq(bit1,bit2,k+1)
    end

    lmer
end


## this function is used for wiring:


function comp_rev(counts)
    a(i,j) = (last(counts[i]) > last(counts[j])) || ((first(counts[i]) == first(counts[j])) && (first(counts[i])> first(counts[j])))
    return  a
end


function graph_creator(kmer_list :: DefaultDict, Alphabet :: Vector{Char}, C :: Int64)
    G = DefaultDict{Vector{DNASeq},macro_node}(0)
    vc = 0
    for x in kmer_list
        xkey , ~ = x
        x_prime_list = read_lmer_from_kmer(xkey,k-1)
        
        for x_prime in x_prime_list
            
            if !([x_prime[1]] in keys(G))
                x_prime_key,~ = x_prime
                u = macro_node()
                u.label = [x_prime_key]
                #u.id = id
                pid = 1
                sid = 1


                for c in Alphabet
                    
                    temp = kmerge(c,x_prime_key, false)
                    node = [DNASeq(vcat(zeros(Int64,64-k+1),temp.bit1[end-k+1:end-1]),vcat(zeros(Int64,64-k+1),temp.bit2[end-k+1:end-1]),k-1)]

                    if  temp in keys(kmer_list) #&& node!= [x_prime_key]
                        u.prefixes[pid] = string_to_DNASeq(string(c))
                        push!(u.prefixes_terminal, false)
                        vc = ceil(Int64,kmer_list[temp]/C)
                        u.prefix_counts[pid] = (kmer_list[temp],vc)
                        u.prefix_terminal = false
                        pid += 1
                        

                    end

                    temp = kmerge(c,x_prime_key, true)
                    node = [DNASeq(vcat(zeros(Int64,64-k+1),temp.bit1[end-k+2:end]),vcat(zeros(Int64,64-k+1),temp.bit2[end-k+2:end]),k-1)]

                    if  temp in keys(kmer_list) # && node != [x_prime_key] 
                        u.suffixes[sid] = string_to_DNASeq(string(c))
                        if length(u.suffixes[sid])==0
                            print("error\n",string_to_DNASeq(string(c))[1],"\n\n")
                        end
                        push!(u.suffixes_terminal, false)
                        vc = ceil(Int64,kmer_list[temp]/C)
                        u.suffix_counts[sid] = (kmer_list[temp],vc)
                        u.suffix_terminal = false
                        sid += 1
                    end
                end
                wiring_prep(u)
                setup_wiring(u)
                G[u.label] = u
            end
            
            
            
            u = macro_node()
        end
        
        
        

    end
#    pref = macro_node()
#    pref.prefix_terminal = true
#    pref.suffix_terminal = false
#    suff = macro_node()
#    suff.prefix_terminal = false
#    suff.suffix_terminal = true

    G

end

function wiring_prep(u :: macro_node)
    u.prefixes[length(u.prefixes)+1] = VTerminal
    push!(u.prefixes_terminal, true)
    u.prefix_counts[length(u.prefixes)] = (-1,-1)
    u.suffixes[length(u.suffixes)+1] = VTerminal
    push!(u.suffixes_terminal, true)
    u.suffix_counts[length(u.suffixes)] = (-1,-1)
end

 
function setup_wiring(u :: macro_node)
    sc, pc = 0,0
    null_sid, null_pid = -1,-1
    for (i,~) in u.suffixes
        sc += last(u.suffix_counts[i])
    end
    for (i,~) in u.prefixes
        pc += last(u.prefix_counts[i])
    end

    for (i,j) in u.suffixes
        if length(j) == 1 && j[1].len == 0
            null_sid = i
            if last(u.suffix_counts[i]) == -1
                 u.suffix_counts[i] = (1,pc-sc)
            end
        end
    end
    for (i,j) in u.prefixes
        if length(j) == 1 & j[1].len == 0
            null_pid = i
            if last(u.prefix_counts[i]) == -1
                u.prefix_counts[i] = (1, sc-pc)
            end
        end
    end

    leftover = sc + last(u.suffix_counts[null_sid])

    last_largest_pid, prefix_begin_pos = -1,-1
    wire_idx = 1
    var_p, var_s = 0,0
    top_p, top_s = 1,1
    p_size = 0
    offset_in_suffix = zeros(Int64, length(u.suffixes))

    indices_s = collect(1:length(u.suffixes))
    indices_p = collect(1:length(u.prefixes))

    indices_s = sort(indices_s, lt = comp_rev(u.suffixes))
    indices_p = sort(indices_p, lt = comp_rev(u.prefixes))
    while leftover > 0
        largest_sid = indices_s[top_s];
        largest_pid = indices_p[top_p];

        
        count = min(last(u.prefix_counts[largest_pid]) - var_p,last(u.suffix_counts[largest_sid]) - var_s)
        push!(u.wire_info[wire_idx], (largest_sid, offset_in_suffix[largest_sid],count))
        if last_largest_pid != largest_pid
            prefix_begin_pos = wire_idx
            last_largest_pid = largest_pid
        end
        p_size += 1
        wire_idx += 1
        leftover -= count
        var_p += count
        var_s += count
        offset_in_suffix[largest_sid] += count
        if var_p == last(u.prefix_counts[largest_pid])
            var_p = 0
            top_p += 1
            u.prefix_begin_info[largest_pid] = (prefix_begin_pos, p_size)
            p_size = 0
        end
  
        if var_s == last(u.suffix_counts[largest_sid])
            var_s = 0
            top_s += 1
        end
  
    end

end
### test
### G = graph_creator(kmer_list,['A','C','G','T'], 5)