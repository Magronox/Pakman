


"""

Deprecated



function IS(G:: DefaultDict{DNASeq,macro_node},k::Int64,l::Int64,Alphabet::Vector{Char})
    print("esesese")
    # k is for kmer in the original graph_construction
    # l is for the length of suffix/prefix at the current graph
    ## Create independent set of G
    I_set = Set()
    #idx = 1
    for (label,u) in G
        confirmed = true
        pre_neighs = []
        succ_neighs = []
        ### finding neighbors
        x = extract_pred(label,k-1)
        for c in Alphabet
            temp = kmerge(c,x,k-2,false)
            if temp in keys(G)
                push!(pre_neighs,temp)
            end
        end
        
        x = extract_succ(label,k-1)
        for c in Alphabet
            temp = kmerge(c,x,k-2,true)
            if temp in keys(G)
                push!(succ_neighs,temp)
            end
        end
        for neigh in pre_neighs
            if vcat(label.bit1,label.bit2) < vcat(neigh.bit1,neigh.bit2)
                confirmed = false
            end
        end
        for neigh in succ_neighs
            if vcat(label.bit1,label.bit2) < vcat(neigh.bit1,neigh.bit2)
                confirmed = false
            end
        end
        if confirmed == true
            push!(I_set,label)
            terminal_1_counts = (1,0)
            terminal_1_prefix = false
            terminal_2_counts = (1,0)
            terminal_key = -1
            likewise = false
            reverse = false

            for neigh in succ_neighs
                print("2\n")
                count_1 = 0
                for (k,v) in G[neigh].prefixes

                    if v == Terminal
                        terminal_1_counts = G[neigh].prefix_counts[k]
                        terminal_key = k
                        terminal_1_prefix = true
                    end
                    
                    idx = 0
                    for (k0,v0) in G[label].prefixes
                        
                        if v0 == Terminal 
                            if terminal_1_prefix
                                terminal_2_counts = G[label].prefix_counts[k0]
                                likewise = true
                                break
                            else
                                terminal_2_counts = G[label].prefix_counts[k0]
                                likewise = false
                                reverse = true
                            end

                        end
                        
                        if length(G[label].prefixes)==1
                            
                            #G[neigh].prefixes[k] = kmerge(v,v0,l)
                        else
                            if idx== 0
                                print(v,"\n")
                                print(v.bit1,"\n")
                                #G[neigh].prefixes[k] = kmerge(v,v0,l)
                                idx += 1
                                count_1 += 1
                            else
                                #G[neigh].prefixes[length(G[neigh].prefixes)+count_1] = kmerge(v,v0,l)
                                count_1 += 1
                            end
                        end
                    end
                    
                    #if length(G[label].prefixes)==1
                    #    break
                    #elseif
                    #end
                end

            end
            """
            """
            for neigh in pre_neighs
                count_1 = 0

                for (k,v) in G[neigh].suffixes
                    idx = 0
                    #if v == Terminal
                    #    terminal_1_counts = G[neigh].suffix_counts[k]
                    #    terminal_1_prefix = false #G[neigh].prefix_terminal
                    #    terminal_key = k
                    #end

                    for (k0,v0) in G[label].suffixes
                        #if v0 == Terminal
                        #    terminal_2_counts = G[neigh].suffix_counts[k0]
                        #    terminal_2_prefix = false #G[neigh].prefix_terminal
                        #    likewise = true
                        #    break
                        #elseif v == Terminal
                        #    break
                        #end

                        if length(G[label].suffixes)==1
                            G[neigh].suffixes[k] = kmerge(v,v0)
                            G[neigh].suffix_counts[]
                        else
                            if idx== 0
                                G[neigh].suffixes[k] = kmerge(v,v0)
                                idx += 1
                                count_1 += 1

                            else
                                G[neigh].suffixes[length(G[neigh].prefixes)+count_1] = kmerge(v,v0)
                                count_1 += 1
                            end
                        end
                    end
                            
                end
"""
"""
                #if reverse == true

                #if likewise
                #    final_terminal_counts = terminal_1_counts + terminal_2_counts
                #    G[neigh].prefix_terminal = false
                #    G[neigh].suffix_terminal = true
                #    G[neigh].suffixes[terminal_key] = Terminal
                #    G[neigh].suffix_counts[terminal_key] = final_terminal_counts
                #else

                

        #    end
        end
    end
    #for label in I_set
    #    delete!(G,label)
    #end
    
    #return G, length(G)-length(I_set)
end








using Distributed

n = 4
addprocs(n)

@everywhere begin
    using Random, DataStructures, DistributedArrays, SharedArrays
end

primitive type fourGbit 8 end

@everywhere function char_to_bits(char :: Char)
    if     char == 'A' return 0
    elseif char == 'C' return 1
    elseif char == 'G' return 2
    elseif char == 'T' return 3
    elseif char == '*' return 42
    else return 100
    end
end
 


@everywhere function bits_to_char(bits::Vector{Int64})
    if     bits == 0  return 'A'
    elseif bits == 1  return 'C'
    elseif bits == 2  return 'G'
    elseif bits == 3  return 'T'
    elseif bits == 42 return '*'
    else              return '*'
    end
end

   
@everywhere function string_inverter(s:: Vector)
    s_temp = ""
    for i in 1:length(s)
        s_temp = string(s_temp,s[i])
    end
    s_temp
end
@everywhere @inline function lmert_convert(s :: String)
    temp = 1
    for i in 1:length(s)
        temp += char_to_bits(s[i])*2^(2*(length(s)-i))
    end
    temp
end



slen = 1000
input = randstring("ACGT",slen)

@everywhere k = 32 # kmer length
@everywhere l = 10 # lmer length

input = only.(split(input,""))
#input = char_to_bits.(input)
#input = DArray((input).reshape())
input_dist = distribute(input)
kmap = DefaultDict{UInt64, Int64}(0)


#### because we can't have shared dictionary I use two SharedArrays
@everywhere lmer_frequency = zeros(Int64,4^l)
@distributed for i in l:length(localpart(input_dist))
    lmer_frequency[lmert_convert(string_inverter(localpart(input_dist[i-9:i])))] += 1
end
    
lmer_freq = @distributed (+) for i in 1:nworkers()
    lmer_frequency
end
    
psize = nworkers()

@everywhere kmer = []
@distributed for i in k:length(localpart(input_dist))
    for j in l:lmert_convert(string_inverter(localpart(input_dist[i-k+1:i])))
    lmer_frequency[lmert_convert(string_inverter(localpart(input_dist[i-9:i])))] += 1
    end
end
 
lockmap_v = SharedArray{UInt64}(4)
kmap_k = SharedArray{}

mnlength = 32

####test
@everywhere using DistributedArrays

a = distribute([[] for _ in procs()]) end

@sync @distributed for i = 1:10
  b = fill(i, 5)
  append!(localpart(a)[1], b) 
  print(localpart(a)," id", "\n") # I swapped push! for append!
end


struct DNASeq
    bit1::BitArray{1}
    bit2::BitArray{1}
  end

## read step

kmer = UInt64(0)
@distributed for i in 1:mnlength
    kmer = (kmer << 2) + char_to_bits(input[i])
    print(bitstring(kmer),"\n")
end

for i in mnlength+1:length(input)
    kmer = (kmer << 2) + char_to_bits(input[i])
    kmap[kmer] += 1
end

n_kmers = length(kmaps)

c_values = SharedArray{Int64}(4^l)
c_keys = SharedArray{Int64}(4^l)

@everywhere count_1 = 0
@everywhere count_2 = 0
@everywhere count_3 = 0

futures = Array{Future}( undef, (nworkers()-1)*(Int(ceil(slen/nworkers()))-k+1))

for i in 2:nworkers()
    for j in 1:Int(ceil(slen/nworkers())-k+1)
        futures[(i-2)*(Int(ceil(slen/nworkers())-k+1))+j] = @spawnat i counts[String(localpart(b)[j:j+k-1])] += 1
    end
                #r = @spawnat i (length(localpart(b)))
                #print(fetch(r),"\n")
            
end


@distributed for i in 1:slen-k+1
                if myid()==1
                    global count_1 +=1
                elseif myid()==2
                    print("a")
                    global count_2 += 1
                elseif myid()==3
                    global count_3 +=1
                else print("wow")
                end
                #print("id",myid(),a[i:i+k-1],"\n")
                    counts[String(a[i:i+k-1])] += 1

            end

@distributed for i in 1:slen-k+1
                if myid()==1
                    global count_1 +=1
                elseif myid()==2
                    print("a")
                    global count_2 += 1
                elseif myid()==3
                    global count_3 +=1
                else print("wow")
                end
                #print("id",myid(),a[i:i+k-1],"\n")
                    counts[String(a[i:i+k-1])] += 1

            end
s=[1,2,3,4,5]
distribute(s)
f = ans
@everywhere for i in s
                print(i)
            end



using BenchmarkTools
r = 1:100
futures = Array{Future}(undef, nworkers())
@btime begin
    for (i,id) in enumerate(workers())
       print("i: ",i," id: ", id,"\n")
       futures[i]
    end
end


a = zeros(10)
@distributed for i = 1:10
    print(myid(),"\n")
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
        leftover = pvc - svc
        while length(suff_dict)>length(pref_dict)
            if leftover > 0 
                push!(u.wire_info.prefix_info[pid],(sid, offset, min(pref_dict[pid], suff_dict[sid])))
                delete!(suff_dict,sid)
                sid, svc = maximum(suff_dict)
                leftover -= svc
            else
                delete!(pref_dict,pid)
                offset = 0
                pid, pvc = maximum(pref_dict)
                leftover = pvc - svc
            end
        end
        offset = 0
        while length(suff_dict)!=0
            push!(u.wire_info.prefix_info[pid],(sid, offset, min(pref_dict[pid], suff_dict[sid])))
            delete!(suff_dict,sid)
            delete!(pref_dict,pid)
            if length(suff_dict)==0
                break
            end
            sid, svc = maximum(suff_dict)
            pid, pvc = maximum(pref_dict)
        end
    else
        pid, pvc = maximum(pref_dict)
        sid, svc = maximum(suff_dict)
        leftover = pvc - svc
        while length(pref_dict)>length(suff_dict)
            if leftover < 0 
                push!(u.wire_info.prefix_info[pid],(sid, offset, min(pref_dict[pid], suff_dict[sid])))
                delete!(pref_dict,pid)
                sid, svc = maximum(pref_dict)
                offset += 1
                leftover += pvc
            else
                delete!(suff_dict,sid)
                offset = 0
                sid, svc = maximum(suff_dict)
                leftover = pvc - svc
            end
        end
        offset = 0
        while length(suff_dict)!=0
            push!(u.wire_info.prefix_info[pid],(sid, offset, min(pref_dict[pid], suff_dict[sid])))
            delete!(suff_dict,sid)
            
            delete!(pref_dict,pid)
            if length(suff_dict)==0
                break
            end
            sid, svc = maximum(suff_dict)
            pid, pvc = maximum(pref_dict)
        end
    end
    else
        pid, pvc = maximum(pref_dict)
        sid, svc = maximum(suff_dict)
        leftover = pvc - svc
        while !isempty(suff_dict)
            while leftover<0 && length(pref_dict)>1
            push!(u.wire_info.prefix_id[pid],(sid, offset, min(pref_dict[pid], suff_dict[sid])))
            offset +=1
            delete!(pref_dict,pid)
            pid, pvc = maximum(pref_dict)
            leftover += pvc
            end
            if length(pref_dict) == 1
                push!(u.wire_info.prefix_id[pid],(sid, offset, min(pref_dict[pid], suff_dict[sid])))
               
            end
            delete!(suff_dict,sid)
            sid, svc = maximum(pref_dict)
        end
    end
"""
"""
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


"""

"""

       


function setup_wiring(u::macro_node)
    sc, pc = 0,0
    
    for (i,~) in u.prefixes
        pc += last(u.prefix_counts[i])
    end
    for (i,~) in u.suffixes
        sc += last(u.suffix_counts[i])
    end
    
    if pc > sc
        @assert(u.suffixes[length(u.suffixes)] == VTerminal)
        u.suffix_counts[length(u.suffixes)] = (1,pc - sc)
        
    else
        
        @assert(u.prefixes[length(u.prefixes)] == VTerminal)
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
                push!(u.wire_info[pid],(sid, offset_in_suffix, suff_dict[sid]))
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
                push!(u.wire_info[pid],(sid, offset_in_suffix, suff_dict[sid]))
                #
                delete!(suff_dict,pid)
                sid,~ = maximum(suff_dict)
                leftover -= suff_dict[sid]
                #offset_in_suffix += count
            end
        end
        

        while(length(suff_dict) != 0)
            push!(u.wire_info[pid],(sid, offset_in_suffix, suff_dict[sid]))
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
                push!(u.wire_info[pid],(sid, offset_in_suffix, suff_dict[sid]))
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
                push!(u.wire_info[pid],(sid, offset_in_suffix, suff_dict[sid]))
                count = min(pref_dict[pid],suff_dict[sid])
                delete!(pref_dict,pid)
                pid,~ = maximum(pref_dict)
                leftover += pref_dict[pid]
                offset_in_suffix += count
            end
        end

        while(length(suff_dict) != 0)
            push!(u.wire_info[pid],(sid, offset_in_suffix, suff_dict[sid]))
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

"""



"""
function iterate_pack(G :: DefaultDict, IS_ :: Set, k:: Int64)
    transfer_nodeInfo = DefaultDict{Vector{DNASeq}, Set{Tuple}}(Set{Tuple}())
    pcontig_list = Vector{}()
    
    
    for node in IS_
        #print("IS node",DNASeq_to_string(node[1]),"\n\n")

        for (i,iset) in G[node].wire_info
            pid = i
            if  G[node].prefixes[pid][end].len > 0 && G[node].prefixes_terminal[i] == false #!(pid in G[node].prefix_terminal_id)  #&& length(G[node].prefixes)>0
                pred_node = pred_neigh(node[1], G[node].prefixes[pid] )
                if pred_node == Nothing
                    print("pred node nothing\n")
                    continue
                end
                
                #print("prefix  and pred node",DNASeq_to_string(G[node].prefixes[pid][1]),"  ",DNASeq_to_string(pred_node[1]),"\n\n")
                if find_pred_ext(G,pred_node[1], node[1]) == Nothing
                    print(pred_node[1],"\n\n",node[1],"\n\n")
                    print(G[pred_node].suffixes,"\n\n\n")
                end
                ssid,pred_ext = find_pred_ext(G,pred_node[1], node[1])
                
            end
            
            for (j, offset, visit_count) in iset
                sid = j
                
                if G[node].prefixes_terminal[pid] && G[node].suffixes_terminal[sid]#(pid in G[node].prefix_terminal_id)  && (sid in G[node].suffix_terminal_id)
                    print("\n\n\n\npleeease\n\n\n")
                    #print(G[node].prefixes[pid],G[node].suffixes[sid], kmerge(G[node].prefixes[pid],node),"\n")
                    contig = kmerge(kmerge(G[node].prefixes[pid],node),G[node].suffixes[sid])
                    push!(pcontig_list,contig)
                else
                    if !(G[node].suffixes[sid] == VTerminal )
                        succ_node = succ_neigh(node[1], G[node].suffixes[sid])

                        if succ_node == Nothing
                            print("succ node nothing error\n")
                            continue
                        end
                        ppid,succ_ext = find_succ_ext(G,node[1], succ_node[1])
                    end
                    

                    
                    if !G[node].prefixes_terminal[pid]
                         #&& pred_ext!= VTerminal ## I added the second condition but seemed necessary
                        
                        if pred_node != node
                            new_pnode_type = false
                            if G[node].suffixes_terminal[sid]
                                new_pnode_type = true
                            else
                                if succ_node == node
                                    new_pnode_type = true
                                else
                                    new_pnode_type = false
                                end
                            end

                            new_ext = kmerge(pred_ext,G[node].suffixes[sid])
                            if length(new_ext)>1 && new_ext[end-1] == Terminal
                                new_ext = new_ext[1:end-1]
                                print(pred_ext,"\n\n")
                                push!(transfer_nodeInfo[node] , (pred_node,pred_ext,new_ext,min(first(G[node].prefix_counts[pid]),first(G[nnode].suffix_counts[ssid])),G[node].wire_info1,ssid,sid,new_pnode_type))
                            elseif !G[node].prefixes_terminal[pid]
                                #print()
                                print("kafte",pred_ext,"\n\n")
                            end
                        end

                    end

                    if !G[node].suffixes_terminal[sid]
                        
                        if succ_node != node
                            new_snode_type = false
                            if G[node].prefixes_terminal[pid]
                                new_snode_type = true
                            else
                                if pred_node == node
                                    new_snode_type = true
                                else
                                    new_snode_type = false
                                end
                            end
                            

                            #print("jafar",G[node].prefixes[pid][1],"\n\n\n")
                            new_ext = kmerge(G[node].prefixes[pid],succ_ext)
                            if length(new_ext)>1 && new_ext[2] == Terminal
                                new_ext = new_ext[2:end]
                            
                            end

                            push!(transfer_nodeInfo[node] , (succ_node,succ_ext,new_ext,min(first(G[node].suffixes[sid]),first(G[nnode].prefix_counts[ppid])),0,ppid,pid,new_snode_type))
                            print("succ_node",DNASeq_to_string(new_ext[1]),"\n\n\n")
                        end
                    end

                end
            end
        end

    end
    
    return pcontig_list, transfer_nodeInfo
end
"""


"""
for (i,iset) in G[node].wire_info
    pid = i
    if !(pid in keys(G[node].prefixes))
        print("kk\n",node,"\n",G[node].wire_info,"\n",pid,"\n",G[node].prefixes)
    end
    print(pid,"\n\n",length(G[node].prefixes[pid]))
    if  G[node].prefixes[pid][end].len > 0 && G[node].prefixes_terminal[i] == false #!(pid in G[node].prefix_terminal_id)  #&& length(G[node].prefixes)>0
        pred_node = pred_neigh(node[1], G[node].prefixes[pid] )
        if pred_node == Nothing
            print("pred node nothing\n")
            continue
        end
        
        #print("prefix  and pred node",DNASeq_to_string(G[node].prefixes[pid][1]),"  ",DNASeq_to_string(pred_node[1]),"\n\n")
        if find_pred_ext(G,pred_node[1], node[1]) == Nothing
            print(pred_node[1],"\n\n",node[1],"\n\n")
            print(G[pred_node].suffixes,"\n\n\n")
        end
        ssid,pred_ext = find_pred_ext(G,pred_node[1], node[1])
        
    end
    
    for (j, offset, count) in iset
        sid = j
        
        if G[node].prefixes_terminal[pid] && G[node].suffixes_terminal[sid]#(pid in G[node].prefix_terminal_id)  && (sid in G[node].suffix_terminal_id)
            print("\n\n\n\npleeease\n\n\n")
            #print(G[node].prefixes[pid],G[node].suffixes[sid], kmerge(G[node].prefixes[pid],node),"\n")
            contig = kmerge(kmerge(G[node].prefixes[pid],node),G[node].suffixes[sid])
            push!(pcontig_list,contig)
        else
            if !(G[node].suffixes[sid] == VTerminal )
                succ_node = succ_neigh(node[1], G[node].suffixes[sid])

                if succ_node == Nothing
                    print("succ node nothing error\n")
                    continue
                end
                ppid,succ_ext = find_succ_ext(G,node[1], succ_node[1])
            end
            

            
            if !G[node].prefixes_terminal[pid]
                 #&& pred_ext!= VTerminal ## I added the second condition but seemed necessary
                
                if pred_node != node
                    new_pnode_type = false
                    if G[node].suffixes_terminal[sid]
                        new_pnode_type = true
                    else
                        if succ_node == node
                            new_pnode_type = true
                        else
                            new_pnode_type = false
                        end
                    end

                    new_ext = kmerge(pred_ext,G[node].suffixes[sid])
                    if length(new_ext)>1 && new_ext[end-1] == Terminal
                        new_ext = new_ext[1:end-1]
                        print(pred_ext,"\n\n")
                        push!(transfer_nodeInfo[node] , (pred_node,pred_ext,new_ext,min(first(G[node].prefix_counts[pid]),first(G[nnode].suffix_counts[ssid])),G[node].wire_info1,ssid,sid,new_pnode_type))
                    elseif !G[node].prefixes_terminal[pid]
                        #print()
                        print("kafte",pred_ext,"\n\n")
                    end
                end

            end

            if !G[node].suffixes_terminal[sid]
                
                if succ_node != node
                    new_snode_type = false
                    if G[node].prefixes_terminal[pid]
                        new_snode_type = true
                    else
                        if pred_node == node
                            new_snode_type = true
                        else
                            new_snode_type = false
                        end
                    end
                    

                    #print("jafar",G[node].prefixes[pid][1],"\n\n\n")
                    new_ext = kmerge(G[node].prefixes[pid],succ_ext)
                    if length(new_ext)>1 && new_ext[2] == Terminal
                        new_ext = new_ext[2:end]
                    
                    end

                    push!(transfer_nodeInfo[node] , (succ_node,succ_ext,new_ext,min(first(G[node].suffixes[sid]),first(G[nnode].prefix_counts[ppid])),0,ppid,pid,new_snode_type))
                    print("succ_node",DNASeq_to_string(new_ext[1]),"\n\n\n")
                end
            end

        end
    end
end
"""


"""
if node_type == 1 && !(new_ext in values(G[nnode].suffixes))


    if index_suf == 0 && !(new_ext in values(G[nnode].suffixes))
        if new_ext[end] == Terminal
            push!(G[nnode].suffix_terminal_id,id)
        end
    #    print("prev_suff",G[node].suffixes[id],"\n")
        G[nnode].suffixes[id] = new_ext
    #    print("new_suff",DNASeq_to_string(new_ext[1]))

        G[nnode].suffix_counts[id] = (first(G[nnode].suffix_counts[id]),visit_count)
        index_suf += 1
    elseif !(new_ext in values(G[nnode].suffixes))
    #    print("prev_suff",DNASeq_to_string(G[node].suffixes[id][1]),"\n\n")
    #    print("new_ext",DNASeq_to_string(new_ext[1]))

        G[nnode].suffixes[1 + maximum(keys(G[nnode].suffixes))] = new_ext
        #push!(G[nnode].wire_info[pid],(maximum(keys(G[nnode].suffixes)),0,visit_count))
        G[nnode].suffix_counts[maximum(keys(G[nnode].suffixes))] = (first(G[nnode].suffix_counts[id]),visit_count)
        if new_ext[end] == Terminal
            G[nnode].suffixes_terminal[maximum(G[nnode].suffixes)] = true
            #push!(G[nnode].suffix_terminal_id,maximum(keys(G[nnode].suffixes)))
        end
    end

elseif !(new_ext in values(G[node].prefixes))
    if !(new_ext in values(G[nnode].prefixes))
        print("shafte",new_ext,"\n",new_ext == VTerminal,"\n\n\n")
    end

    if index_pref == 0 && !(new_ext in values(G[nnode].prefixes))
        if new_ext[end] == Terminal
            G[nnode].prefixes_terminal[id] = true
            #push!(G[nnode].prefix_terminal_id,id)
        end
        G[nnode].prefixes[id] = new_ext
        G[nnode].prefix_counts[id] = (first(G[nnode].prefix_counts[id]), visit_count)
        index_pref += 1

    else
        G[nnode].prefixes[1+maximum(keys(G[nnode].prefixes))] = new_ext
        G[nnode].prefix_counts[maximum(keys(G[nnode].prefixes))] = (first(G[nnode].prefix_counts[id]), visit_count)
        #push!(G[nnode].wire_info[pid],(maximum(keys(G[nnode].suffixes)),1,visit_count))
        if new_ext[end] == Terminal
            G[nnode].prefixes_terminal[maximum(keys(G[nnode].suffixes))] = true
            #push!(G[nnode].prefix_terminal_id,maximum(keys(G[nnode].suffixes)))
        end
    end

end
"""



"""

function rewire(G, rewire_list)
    for node in rewire_list
        u = G[node]
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
                    push!(u.wire_info[pid],(sid, offset_in_suffix, suff_dict[sid]))
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
                    push!(u.wire_info[pid],(sid, offset_in_suffix, suff_dict[sid]))
                    #
                    delete!(suff_dict,pid)
                    sid,~ = maximum(suff_dict)
                    leftover -= suff_dict[sid]
                    #offset_in_suffix += count
                end
            end
            
    
            while(length(suff_dict) != 0)
                push!(u.wire_info[pid],(sid, offset_in_suffix, suff_dict[sid]))
                delete!(pref_dict, pid)
                delete!(suff_dict,sid)

                if length(suff_dict) != 0
                    if length(pref_dict)==0
                        print(pref_dict,"\n",u.wire_info,"\n\n")
                    end
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
                    push!(u.wire_info[pid],(sid, offset_in_suffix, suff_dict[sid]))
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
                    push!(u.wire_info[pid],(sid, offset_in_suffix, suff_dict[sid]))
                    count = min(pref_dict[pid],suff_dict[sid])
                    delete!(pref_dict,pid)
                    pid,~ = maximum(pref_dict)
                    leftover += pref_dict[pid]
                    offset_in_suffix += count
                end
            end
    
            while(length(suff_dict) != 0)
                push!(u.wire_info[pid],(sid, offset_in_suffix, suff_dict[sid]))
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
end
"""
"""
        if fan_out
            
            pid, pvc = maximum(pref_dict)
            sid, svc = maximum(suff_dict)
            leftover = pvc - svc
            if leftover > 0
            while length(suff_dict)>length(pref_dict)
                print("pref ", length(pref_dict),"\n")
                print("suff ", length(suff_dict),"\n")
                if leftover > 0 
                    push!(u.wire_info[pid],(sid, offset, min(pref_dict[pid], suff_dict[sid])))
                    delete!(suff_dict,sid)
                    sid, svc = maximum(suff_dict)
                    leftover -= svc
                else
                    delete!(pref_dict,pid)
                 
                    pid, pvc = maximum(pref_dict)
                    leftover = pvc - svc
                    offset = 0

                    if length(pref_dict)!= 0 
                        print("jafar\n")
                    elseif length(suff_dict) != 0 
                        print(length(suff_dict))
                        print(u.prefixes,"\n")
                        print("suffs \n",u.suffixes,"\n\n")
                        print(u.wire_info)
                        print("\n\n",suff_dict,"\n\n")
                        print("suff wiring error\n")
                        @assert(false)
                    elseif length(pref_dict) != 0
                        print("second pref error\n")
                    else
                        return
                    end
                end
            end
            offset = 0
            while length(suff_dict)!=0
                push!(u.wire_info[pid],(sid, offset, min(pref_dict[pid], suff_dict[sid])))
                delete!(suff_dict,sid)
                delete!(pref_dict,pid)
                if length(suff_dict)==0
                    break
                end
                sid, svc = maximum(suff_dict)
                pid, pvc = maximum(pref_dict)
            end
        else
            pid, pvc = maximum(pref_dict)
            sid, svc = maximum(suff_dict)
            leftover = pvc - svc
            while length(pref_dict)>length(suff_dict)
                print("pref ", length(pref_dict),"\n")
                print("suff ", length(suff_dict),"\n")
                if leftover < 0 
                    push!(u.wire_info[pid],(sid, offset, min(pref_dict[pid], suff_dict[sid])))
                    delete!(pref_dict,pid)
                    pid, pvc = maximum(pref_dict)
                    offset += 1
                    leftover += pvc
                else
                    delete!(suff_dict,sid)
                    offset = 0
                    if length(suff_dict)!= 0
                        sid, svc = maximum(suff_dict)
                        leftover = pvc - svc
                    elseif length(pref_dict) != 0
                        print(length(pref_dict))
                        print(u.wire_info)
                        print("\n\n",pref_dict,"\n\n")
                        print("wiring error\n")
                        print("pref\n",u.prefixes,"\n\n")
                        print(u.suffixes,"\n")
                        @assert(false)
                    else
                        return
                    end
                end
            end
            offset = 0
            while length(suff_dict)!=0
                push!(u.wire_info[pid],(sid, offset, min(pref_dict[pid], suff_dict[sid])))
                delete!(suff_dict,sid)
                
                delete!(pref_dict,pid)
                if length(suff_dict)==0
                    break
                end
                sid, svc = maximum(suff_dict)
                pid, pvc = maximum(pref_dict)
            end
        end
    end
end
## Test
## slen = 1000
## input = randstring("ACGT",slen)

## DNA_seq = string_to_DNASeq(input)
## k=3
## kmer_list = read_kmer(DNA_seq, length(input),k)

## G = graph_creator(kmer_list,['A','C','G','T'], 5)
## I = IS(G,['A','C','G','T'],k)
## G_new, ~ = compact_graph!(G,k,floor(Int64,length(G)*4/5))
"""

"""
for (sid,offset_in_suffix,sz) in mn.prefix_begin_info[pid]
    
        #if inter
        
    if offset_in_prefix < internal_off
        off_in_wire = offset_in_prefix - internal_off
    end
    next_off = offset_in_suffix + off_in_wire
    
    freq_in_wire = min(freq_rem,(sz- off_in_wire))
    contig_new = kmerge(contig, mn.suffixes[sid])
    
    if mn.suffix_terminal_id == sid
        
        push!(output, contig_new)
    else
        next_mn = G[succ_neigh(mn.label,mn.suffixes[sid])]
        next_prefix_id, ~ = find_succ_ext(G, mn.label, next_mn.label)
        walk!(G,contig_new,freq_in_wire,next_off,next_mn,next_prefix_id,output)
    end
    freq_rem -= freq_in_wire
    
    internal_off = sz
    

end


end


"""


