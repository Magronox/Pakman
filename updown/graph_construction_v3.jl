include("kmer_counting_v3.jl")


####====
K = 32
coverage = 10
####====



Base.@kwdef mutable struct macro_node
    
    #label = DNASeq
    affix = Vector{DNASeq}() #DefaultDict{Int64, Vector{DNASeq}}(DNASeq[])
    isPrefix = false
    isTerminal = false
    num_wires = 0
    wire_idx = 0
    count = (-1,-1)

end

Base.@kwdef mutable struct WireNode 
    sid = 0
    offset = 0
    count = 0

end

Base.@kwdef mutable struct ModifiedNode 
    old_affix  = Vector{DNASeq} #DefaultDict{Int64, Vector{DNASeq}}(DNASeq[])
    new_affix  = Vector{DNASeq} #DefaultDict{Int64, Vector{DNASeq}}(DNASeq[])
    isPrefix   = false
    isTerminal = false
    num_wires  = 0
    wire_index = 0
    count      = (-1, -1)
end





## define terminal 

Terminal = DNASeq(BitArray{1}(undef,64),BitArray{1}(undef,64),0)

VTerminal = [Terminal]

### to get the maximum values in dictionary
Base.isless(p::Pair, q::Pair) = isless(p.second,q.second)


### macro node equality check
function Base.isequal(mn1::Vector{T}, mn2::Vector{T}) where T <: macro_node
    return mn1.affix==mn2.affix && mn1.isPrefix == mn2.isPrefix && mn1.isTerminal == mn2.isTerminal 
end
Base.:(==)(mn1 ::macro_node, mn2 :: macro_node) =  mn1.affix==mn2.affix && mn1.isPrefix == mn2.isPrefix && mn1.isTerminal == mn2.isTerminal 




function Comp_rev(nodes :: Vector{macro_node})
    a(i,j) = (nodes[i].isPrefix!=nodes[j].isPrefix) ? nodes[i].isPrefix : (last(nodes[i].count)!=last(nodes[j].count) ? last(nodes[i].count)>last(nodes[j].count) : (first(nodes[i].count)!= first(nodes[i].count) ? first(nodes[i].count)>first(nodes[j].count) : (nodes[i].affix >nodes[j].affix)))
    return a
end
 

function construct_macronodes(kmer_list :: DefaultDict{DNASeq, Int64}, min_index :: Int64)
    @show "K is equal to $K and coverage is $coverage"
    G = DefaultDict{mn_label,Vector{macro_node}}(Vector{macro_node}())
    idx = 0
    for (kmer, count) in kmer_list
        idx += 1
        if count >= min_index 
            node_1,node_2 = read_from_kmer(kmer)
            suffix_mn = macro_node()
            prefix_mn = macro_node()
            push!(prefix_mn.affix, kmer_seq(BitArray(kmer.bit1[1]),BitArray(kmer.bit2[1]),1))
            prefix_mn.isPrefix = true  
            prefix_mn.count    = (count, ceil(Int64, count/coverage))
            push!(get!(G, node_2, macro_node[]), prefix_mn)
            #push!(G[node_2], prefix_mn)
            #end
            push!(suffix_mn.affix,kmer_seq(BitArray(kmer.bit1[end]),  BitArray(kmer.bit2[end]),  1))
            suffix_mn.isPrefix = false
            suffix_mn.count    = (count, ceil(Int64, count/coverage))
            #if !(suffix_mn in G[node_2])
            push!(get!(G, node_1, macro_node[]), suffix_mn)
            #end
            
            
            
        end
    end
    return G
end

function finish_wiremaps(G :: DefaultDict{mn_label,Vector{macro_node}})
    wire_dict = DefaultDict{mn_label,Vector{WireNode}}(Vector{WireNode}())
    idx = 0
    for (label,node) in G
        idx += 1
        if idx % 100000 == 0
            print("Wiring\n")
        end
        prefix_mn = macro_node()                      
        suffix_mn = macro_node()
        prefix_mn.isPrefix = true
        suffix_mn.isPrefix = false
        prefix_mn.isTerminal = true
        suffix_mn.isTerminal = true
        push!(get!(G, label, macro_node[]), suffix_mn)
        push!(get!(G, label, macro_node[]), prefix_mn)
        for i in 1:length(node)
            push!(get!(wire_dict, label, WireNode[]), deepcopy(WireNode()))
        end
    end
    return wire_dict
end

function save_build_info(input_file:: String,output_file_1 :: String)
    kmer_list =  DefaultDict{DNASeq, Int64}(0)
    idx = 0
    print("Started Reading\n")
    for (id,seq) in FastaReader(input_file)
        if id== "ERCC-00138"
            continue
        end
        idx += 1
        input = string_to_DNASeq(seq)
        #print(seq,"\n")
        temp_kmers = read_kmer(input,32)
        for (t,v) in temp_kmers
            kmer_list[t] += v
        end
        if idx %100000 == 0
            print("reading\n")
        end
    end
    
    print("Finished Reading kmers\n")
    @time G = construct_macronodes(kmer_list,2)
    @time wire_dict = finish_wiremaps(G)
    @time RewireMN(G,wire_dict)

    open(output_file_1, "w") do file_
        for (i,j) in G
            print(file_,"Vertex, $(seq_to_int64(kmer_seq(i.bit1,i.bit2,31))), 0\n")
            for node in j
                if length(node.affix)>0
                    print(file_,"Edge, $(seq_to_int64(kmer_seq(i.bit1,i.bit2,31))), $(node.isTerminal),$(first(node.count)), $(last(node.count)), $(node.wire_idx), $(node.num_wires),1, $(node.affix)\n")
                else
                    print(file_,"Edge, $(seq_to_int64(kmer_seq(i.bit1,i.bit2,31))), $(Int(node.isTerminal)), $(first(node.count)), $(last(node.count)), $(node.wire_idx), $(node.num_wires), 1, $(node.affix)\n")
                
                end
            end
        end
    end

end

function save_wire_info(input_file:: String,output_file_1 :: String, output_file_2 :: String)
    kmer_list =  DefaultDict{DNASeq, Int64}(0)
    idx = 0
    n_edges = 0
    print("Started Reading\n")
    for (id,seq) in FastaReader(input_file)
        if id== "ERCC-00138"
            continue
        end
        idx += 1
        input = string_to_DNASeq(seq)
        #print(seq,"\n")
        temp_kmers = read_kmer(input,32)
        for (t,v) in temp_kmers
            kmer_list[t] += v
        end
        if idx %100000 == 0
            print("reading\n")
        end
    end
    
    print("Finished Reading kmers\n")
    @time G = construct_macronodes(kmer_list,2)
    @time wire_dict = finish_wiremaps(G)
    @time RewireMN(G,wire_dict)

    open(output_file_1, "w") do file_
        for (i,j) in G
            print(file_,"Vertex, $(seq_to_int64(kmer_seq(i.bit1,i.bit2,31))), 0\n")
            
            for node in j
                n_edges += 1
                if length(node.affix)>0 
                    #print(file_,"nodeii $(node.wire_idx) \n")
                    print(file_,"Edge, $(seq_to_int64(kmer_seq(i.bit1,i.bit2,31))), $(wire_dict[i][node.wire_idx+1].sid), $(wire_dict[i][node.wire_idx+1].offset), $(wire_dict[i][node.wire_idx+1].count)\n")
                else
                    print(file_,"Edge, $(seq_to_int64(kmer_seq(i.bit1,i.bit2,31))), $(wire_dict[i][node.wire_idx+1].offset), $(wire_dict[i][node.wire_idx+1].count)\n")
                
                end
                
            end
        end
    end
    open(output_file_2, "w") do file_
        print(file_, "number of nodes $(length(G))\nnumber of edges $n_edges")
    end

end

function ModifyMN()
    for mod in ModifiedNodes
        found = false
        for mn_node in values(G)
            if mod.isPrefix != mn_node.isPrefix
                continue   
            end    
            if mod.old_affix != mn_node.affix
                continue
            end
            found = true
   
            if mn_node.isTerminal
                tmp = macro_node()
                tmp.affix      = mod.new_affix
                tmp.isPrefix   = mod.isPrefix
                tmp.isTerminal = mod.isTerminal
                tmp.num_wires  = mod.num_wires
                tmp.wire_index = mod.wire_index
                tmp.count      = mod.count
                macroNodes.push_back(tmp)      
                wire_dict[key] = WireNode()    
      
            else                                                      
                mn_node.affix      = mod.new_affix                  
                mn_node.isPrefix   = mod.isPrefix
                mn_node.isTerminal = mod.isTerminal
                mn_node.num_wires  = mod.num_wires
                mn_node.wire_index = mod.wire_index
                mn_node.count      = mod.count
            end
      
            break
        end

        if !found
            tmp = macro_node()                                            
            tmp.affix      = mod.new_affix
            tmp.isPrefix   = mod.isPrefix
            tmp.isTerminal = mod.isTerminal
            tmp.num_wires  = mod.num_wires
            tmp.wire_index = mod.wire_index
            tmp.count      = mod.count
            G[key] = tmp    
            wire_dict[key] = WireNode()                       
        end

    end
end

function RewireMN(G, wire_dict)

    for (mn,nodes) in G
        wireNodes = wire_dict[mn]
        pc = 0
        sc = 0
        #index = -1
        top_prefix = 0
        top_suffix = 0
        null_prefix_id = 0
        null_suffix_id = 0
        p_count = 0
        for index in 1:length(nodes)
            node = nodes[index]
            if node.isPrefix                            
                if length(node.affix) == 0
                    null_prefix_id = index -1
                else 
                    pc += last(node.count)
                end
                p_count += 1
            else                               
                if top_suffix == 0
                    top_suffix = index -1
                end 
                if length(node.affix) == 0
                    null_suffix_id = index -1
                else 
                    sc += last(node.count)
                end
            end
        end


        nodes[null_prefix_id+1].count = (1, max(sc - pc, 0))
        nodes[null_suffix_id+1].count = (1, max(pc - sc, 0))
        

        indices = collect(1:length(nodes)) 
        indices = sort(indices, lt = Comp_rev(nodes))
        wire_idx = 0
        num_wires = 0
        var_p = 0
        var_s = 0
        offset_in_suffix = 0
        prefix_begin = 9223372036854775807
        last_prefix_id = 9223372036854775807
        leftover = sc + last(nodes[null_suffix_id+1].count)
        top_suffix = p_count 
        while (leftover > 0) 
            
            prefix_id = indices[top_prefix+1] - 1
            suffix_id = indices[top_suffix+1] - 1
            prefix_count = last(nodes[prefix_id+1].count) - var_p
            suffix_count = last(nodes[suffix_id+1].count) - var_s
            count = min(prefix_count, suffix_count)
        
            wireNodes[wire_idx+1].sid    =  suffix_id
            wireNodes[wire_idx+1].offset = offset_in_suffix
            wireNodes[wire_idx+1].count  = count
            if last_prefix_id != prefix_id
                prefix_begin = wire_idx
                last_prefix_id = prefix_id
            end
            #@show last_prefix_id,prefix_id,wire_idx,prefix_begin
            wire_idx += 1
            num_wires += 1
            var_p    += count
            var_s    += count
            
            leftover -= count
            offset_in_suffix += count
        
            if var_p == last(nodes[prefix_id+1].count)
                #@show prefix_begin
               nodes[prefix_id+1].num_wires  = num_wires
               nodes[prefix_id+1].wire_idx = prefix_begin
               var_p      = 0
               num_wires  = 0
               top_prefix += 1
            end
        
            if var_s == last(nodes[suffix_id+1].count)
               offset_in_suffix = 0
               var_s      = 0
               top_suffix += 1
            end

        end
        #print("done\n")
    end
end
