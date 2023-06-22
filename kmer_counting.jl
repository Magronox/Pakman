


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
"""