include("kmer_counting_v2.jl")






Base.@kwdef mutable struct macro_node
    #length k-1
    #id = 0
    label = DNASeq#BitArray{1}(undef,64)
    prefixes = DefaultDict{Int64, Vector{DNASeq}}(DNASeq[])
    prefix_counts = Dict{Int64, Tuple}()
    prefix_terminal_id = []
    prefix_terminal = false
    suffixes = DefaultDict{Int64, Vector{DNASeq}}(DNASeq[])
    suffix_counts = Dict{Int64,Tuple}()
    suffix_terminal = false
    suffix_terminal_id = []
    wire_info = wire_node()
end


### to use macro_nodes in dict:

function Base.isequal(a :: macro_node, b :: macro_node)
    return Base.isequal(a.label,b.label)
end
function Base.hash(a :: macro_node)
    return Base.hash(a.id)
end

Base.@kwdef mutable struct wire_node
    prefix_info = DefaultDict{Int64, Set{Tuple}}(Set())
    len = 0
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
                        vc = ceil(Int64,kmer_list[temp]/C)
                        u.suffix_counts[sid] = (kmer_list[temp],vc)
                        u.suffix_terminal = false
                        sid += 1
                    end
                end
                
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



        


function setup_wiring(u::macro_node)
    sc, pc = 0,0
    
    for (i,~) in u.prefixes
        pc += last(u.prefix_counts[i])
    end
    for (i,~) in u.suffixes
        sc += last(u.suffix_counts[i])
    end
    
    if pc > sc
        u.suffixes[length(u.suffixes)+1] = VTerminal
        
        push!(u.suffix_terminal_id,length(u.suffixes))
        u.suffix_terminal = true
        u.suffix_counts[length(u.suffixes)] = (1,pc - sc)
        
    else
        
        u.prefixes[length(u.prefixes)+1] = VTerminal
        push!(u.prefix_terminal_id,length(u.prefixes))
        u.prefix_terminal = true
        u.prefix_counts[length(u.prefixes)] = (1, sc - pc)
    end

    pref_dict = Dict(keys(u.prefix_counts) .=> getfield.(values(u.prefix_counts),2))
    suff_dict = Dict(keys(u.suffix_counts) .=> getfield.(values(u.suffix_counts),2))
    offset_in_suffix = 0

    if length(u.suffixes)>length(u.prefixes)
        pid,~ = maximum(pref_dict)
        sid,~ = maximum(suff_dict)
        leftover = pref_dict[pid] - suff_dict[sid]
        offset_in_suffix = 0

        while(length(suff_dict)>length(pref_dict))

            if leftover <= 0
                push!(u.wire_info.prefix_info[pid],(sid, offset_in_suffix, suff_dict[sid]))
                #count = min(pref_dict[pid],suff_dict[sid])
                delete!(pref_dict, pid)
                delete!(suff_dict, sid)
                if length(suff_dict) != 0
                    pid,~ = maximum(pref_dict)
                    sid,~ = maximum(suff_dict)
                    leftover = pref_dict[pid] - suff_dict[sid]
                end
                offset_in_suffix = 0
            else
                push!(u.wire_info.prefix_info[pid],(sid, offset_in_suffix, suff_dict[sid]))
                #
                delete!(suff_dict,pid)
                sid,~ = maximum(suff_dict)
                leftover -= suff_dict[sid]
                #offset_in_suffix += count
            end
        end
        

        while(length(suff_dict) != 0)
            push!(u.wire_info.prefix_info[pid],(sid, offset_in_suffix, suff_dict[sid]))
            delete!(pref_dict, pid)
            delete!(suff_dict,sid)
            if length(suff_dict) != 0
                pid,~ = maximum(pref_dict)
                sid,~ = maximum(suff_dict)
                leftover = pref_dict[pid] - suff_dict[sid]
            end
            offset_in_suffix = 0
        end

    else

        pid,~ = maximum(pref_dict)
        sid,~ = maximum(suff_dict)
        leftover = pref_dict[pid] - suff_dict[sid]
        count = min(pref_dict[pid],suff_dict[sid])

        while(length(pref_dict)>length(suff_dict))
            
            if leftover >= 0
                push!(u.wire_info.prefix_info[pid],(sid, offset_in_suffix, suff_dict[sid]))
                offset_in_suffix += count
                #count = min(pref_dict[pid],suff_dict[sid])
                delete!(pref_dict, pid)
            
                if length(suff_dict) != 1
                    delete!(suff_dict, sid)
                    offset_in_suffix = 0
                end
                if length(suff_dict) != 0
                    pid,~ = maximum(pref_dict)
                    sid,~ = maximum(suff_dict)
                    leftover = pref_dict[pid] - suff_dict[sid]
                    count = min(pref_dict[pid],suff_dict[sid])

                end
                
            else
                push!(u.wire_info.prefix_info[pid],(sid, offset_in_suffix, suff_dict[sid]))
                count = min(pref_dict[pid],suff_dict[sid])
                delete!(pref_dict,pid)
                pid,~ = maximum(pref_dict)
                leftover += pref_dict[pid]
                offset_in_suffix += count
            end
        end

        while(length(suff_dict) != 0)
            push!(u.wire_info.prefix_info[pid],(sid, offset_in_suffix, suff_dict[sid]))
            delete!(pref_dict, pid)
            delete!(suff_dict,sid)
            if length(suff_dict) != 0
                pid,~ = maximum(pref_dict)
                sid,~ = maximum(suff_dict)
                leftover = pref_dict[pid] - suff_dict[sid]
            end
            offset_in_suffix = 0
        end



    end

end
### test
### G = graph_creator(kmer_list,['A','C','G','T'], 5)