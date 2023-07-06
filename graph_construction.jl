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
    prefix_info = DefaultDict{Int64, Set{Tuple}}(Set{Tuple}())
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
           
            x_prime_key,~ = x_prime
            #x_prime_key, ~ = x_prime_key
            #id += 1
            u = macro_node()
            u.label = [x_prime_key]
            #u.id = id
            pid = 1
            sid = 1
            for c in Alphabet
                
                temp = kmerge(c,x_prime_key, false)

                if temp in keys(kmer_list)
                    u.prefixes[pid] = string_to_DNASeq(string(c))
                    vc = ceil(Int64,kmer_list[temp]/C)
                    u.prefix_counts[pid] = (kmer_list[temp],vc)
                    u.prefix_terminal = false
                    pid += 1

                end

                temp = kmerge(c,x_prime_key, true)
                if temp in keys(kmer_list)
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
            G[u.label] = u
            setup_wiring(u)
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



        
function setup_wiring(u :: macro_node)
    terminal_key = 0

    if isempty(u.prefixes)
        prefix_info_2 = zeros(Int64,0)
        prefix_info_3 = zeros(Int64,0)
    else
        prefix_info_2 = zeros(Int64, maximum(keys(u.prefixes)))
        prefix_info_3 = zeros(Int64, maximum(keys(u.prefixes)))
    end
    ### initialize count for prefixes wire_info
    #count = min(maximum(last.(values(u.prefix_counts))),maximum(last.(values(u.suffix_counts))))
    #for pid in keys(u.prefixes)
    #    
    #    prefix_info_3[pid] = count
        #u.wire_info[pid][3] = count
    #end

    sc = 0
    pc = 0
    #null_sid = -1
    #null_pid = -1
    for temp in keys(u.suffixes)
        sc += u.suffix_counts[temp][2]
    end
    for temp in keys(u.prefixes)
        pc += u.prefix_counts[temp][2]
    end

    ## initialize wire_info len and offset_in_suffix
    u.wire_info.len = length(u.prefixes) + length(u.suffixes) + 1
    #for i in u.prefixes
    #    pid,~ = i
    #    prefix_info_2[pid] = length(u.suffixes)
        #u.wire_info[id][2] = length(u.suffixes)
    #end

    if pc>sc
        u.prefix_terminal = false
        u.suffix_terminal = true
        
        if isempty(u.suffixes)
            sid = 1
        else
            sid = maximum(keys(u.suffixes)) + 1
        end

        u.suffixes[sid] = VTerminal
        push!(u.suffix_terminal_id,sid)
        u.suffix_counts[sid] = (1, pc- sc)
        terminal_key = sid


    elseif sc>pc
        
        u.prefix_terminal = true
        u.suffix_terminal = false
        
        if isempty(u.prefixes)
            pid = 1
        else
            pid = maximum(keys(u.prefixes)) + 1
        end
        push!(u.prefix_terminal_id , pid)
        #print(u.prefixes[pid],"\n")
        u.prefixes[pid] = VTerminal
        u.prefix_counts[pid] = (1, sc - pc)
        push!(prefix_info_2,0)
        push!(prefix_info_3,0)

        terminal_key = pid
    end

    ## greedy matching
    
    pref_dict = Dict(keys(u.prefix_counts) .=> getfield.(values(u.prefix_counts),2))
    suff_dict = Dict(keys(u.suffix_counts) .=> getfield.(values(u.suffix_counts),2))
    #set_pid = Set(collect(1:length(u.prefixes)))
    
    
    #phase of the greedy algorithm
    ##
    fan_out = (length(u.suffixes) >= length(u.prefixes))
    #finished = isempty(pref_dict) && isempty(suff_dict)
    offset = 0
    finished = false

    if fan_out
        
        pid, pvc = maximum(pref_dict)
        sid, svc = maximum(suff_dict)
        #print(pid)
        #push!(u.wire_info.prefix_info[pid] ,(1,1,1))
        push!(u.wire_info.prefix_info[pid] ,(sid, offset, min(pref_dict[pid], suff_dict[sid])))
        finished = isempty(pref_dict) && isempty(suff_dict)
        leftover = svc - pvc
        delete!(suff_dict,sid)
        
        #pop!(set_sid,sid)

        while(!finished)

            if (length(pref_dict)-1) == length(suff_dict)
                
                delete!(pref_dict,pid)
                for i in 1:length(pref_dict)
                    pid, pvc = maximum(pref_dict)
                    sid, svc = maximum(suff_dict)
                    push!(u.wire_info.prefix_info[pid] ,(sid, offset, min(pref_dict[pid], suff_dict[sid])))
                    delete!(suff_dict,sid)
                    delete!(pref_dict,pid)
                end
                
                break

            elseif leftover >= 0
                
                delete!(pref_dict,pid)
                #pop!(set_pid,pid)
                if isempty(pref_dict) && isempty(suff_dict)
                    break 
                else
                    
                    pid, pvc = maximum(pref_dict)
                    sid, svc = maximum(suff_dict)
                    push!(u.wire_info.prefix_info[pid] ,(sid, offset, min(pref_dict[pid], suff_dict[sid])))
                    finished = isempty(pref_dict) && isempty(suff_dict)
                    leftover = svc - pvc
                    delete!(suff_dict,sid)
                    
                end
            else 
                
                
                sid, svc = maximum(suff_dict)
                push!(u.wire_info.prefix_info[pid] ,(sid, offset, min(pref_dict[pid], suff_dict[sid])))
                finished = isempty(pref_dict) && isempty(suff_dict)
                leftover += svc
                delete!(suff_dict,sid)
                
                
            end
        end
    else

        pid, pvc = maximum(pref_dict)
        sid, svc = maximum(suff_dict)
        push!(u.wire_info.prefix_info[pid] ,(sid, offset, min(pref_dict[pid], suff_dict[sid])))
        leftover = svc - pvc
        delete!(pref_dict,pid)
        offset += 1
        #finished = isempty(pref_dict) && isempty(suff_dict)
        
        #pop!(set_sid,sid)

        while(!finished)

            if length(pref_dict) == (length(suff_dict)-1)
                delete!(suff_dict,sid)
                offset = 0 
                for i in 1:length(pref_dict)
                    pid, pvc = maximum(pref_dict)
                    sid, svc = maximum(suff_dict)
                    push!(u.wire_info.prefix_info[pid] ,(sid, offset, min(pref_dict[pid], suff_dict[sid])))
                    delete!(suff_dict,sid)
                    delete!(pref_dict,pid)
                    
                end
                break

            elseif leftover <= 0
                delete!(suff_dict,sid)
                #pop!(pref_dict,pid)
                if isempty(pref_dict) && isempty(suff_dict)
                    break 
                else
                    pid, pvc = maximum(pref_dict)
                    sid, svc = maximum(suff_dict)
                    push!(u.wire_info.prefix_info[pid] ,(sid, offset, min(pref_dict[pid], suff_dict[sid])))
                    leftover = svc - pvc
                    delete!(pref_dict,pid)
                    offset += 1
                    finished = isempty(pref_dict) && isempty(suff_dict)
                end
            else 
                
                pid, pvc = maximum(pref_dict)
                push!(u.wire_info.prefix_info[pid] ,(sid, offset, min(pref_dict[pid], suff_dict[sid])))
                leftover -= pvc
                delete!(pref_dict,pid)
                offset += 1
                finished = isempty(pref_dict) && isempty(suff_dict)
            end
        end


    end



    
end


### test
### G = graph_creator(kmer_list,['A','C','G','T'], 5)
