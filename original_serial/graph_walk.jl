include("graph_compaction.jl")



function run_walk(G :: DefaultDict, pcontig_list:: Vector, analysis::Bool=false)
    output = Set()
    begin_kmer_list = copy(keys(G))
    pid = -1
    for node in begin_kmer_list
        #print(" node ", DNASeq_to_string(node[1]),"\n")
        mn = G[node]
        for i in 1:length(mn.prefixes_terminal)
            if mn.prefixes_terminal[i]
                pid = i
                #print("pid ",i,"\n")
                break
            end

        end

        if pid != -1
            freq = last(mn.prefix_counts[pid])
            contig = mn.prefixes[pid]
            contig = kmerge(contig, node)
           
            if analysis
                print("walk node ", DNASeq_to_string(node[1]),"\n")
            end
            #print("calling walk with contig ",DNASeq_to_string(contig[1]), " , freq ", freq, " for node ", DNASeq_to_string(mn.label[1]),"\n\n")
            walk!(G, contig, freq, 0, mn, pid, output,analysis)
            #print("new contig\t",DNASeq_to_string(contig[1]),"\n")

        else
            print("Error Prefixes\n")
        end
    end
    
    "for contig in pcontig_list
        walk!(G, contig, freq, 0, mn, pid, output)
    end"
    #push!(output,pcontig_list)
    
    return output
end
contig_list = []
function walk!(G:: DefaultDict, pcontig:: Vector{DNASeq}, freq, offset_in_prefix :: Int64 , mn :: macro_node, pid :: Int64, output :: Set, analysis :: Bool = false)
    #pc_size = (length(contig)-1)*64 + contig[1].len
    
    freq_rem = freq
    internal_off = 0
    off_in_wire = -1
    if length(pcontig)>(500000/64)
        print("Fails with contig size ",(length(pcontig)-1)*64+pcontig[1].len)
    end

    #if isempty(mn.wire_info[pid])
    #    print("pid ",pid,"\n",mn.wire_info,"\n", mn.prefix_begin_info,"\n\n")
    #end
    #print("started \t pid ", pid, " and offset in prefix ", offset_in_prefix," wire info",mn.wire_info,"\t")
    
    sdf = 0
    for (id,offset_in_suffix,count) in mn.wire_info[pid]
        sdf += 1
        internal_off += mn.prefix_begin_info[pid][2]-1
        sz = copy(count)
        
        if (internal_off + sz) <= offset_in_prefix || internal_off > (offset_in_prefix + freq_rem)
            continue
        end

        off_in_wire = offset_in_prefix <= internal_off ? 0 : offset_in_prefix - internal_off 
        next_off = offset_in_suffix + off_in_wire#id - 1
        freq_in_wire = min(freq_rem, sz - off_in_wire )
        #print("freq in wire ", freq_in_wire, " off in wire ",off_in_wire,"\n")
        temp = false
        
        pcontig = kmerge(pcontig, mn.suffixes[id])
       
        
        if mn.suffixes_terminal[id]
            
            print("contig found\n")
            push!(output, pcontig)
            #push!(pcontig, contig)
            continue
        else
            succ_ext = mn.suffixes[id]
            key = mn.label
            succ_node = succ_neigh(key[1],succ_ext)
            #ppid, succ_ext = find_succ_ext(G,key[1], succ_node[1])
            if !(succ_node in keys(G))
                print("Error key not found\n")
                continue
            end
            
            next_mn = G[succ_node]
            next_prefix_id, ~ = find_succ_ext(G, mn.label, next_mn.label)
            @assert(succ_node == next_mn.label)
            if freq_in_wire > 0 #&& !(pcontig in output)
                if analysis
                    print("continuing the walk with node ", DNASeq_to_string(next_mn.label[1]),"\n")
                end
                walk!(G, pcontig, freq_in_wire, next_off, next_mn, next_prefix_id, output, analysis);
            end 
               
            
        end
       
        freq_rem -= freq_in_wire

    end
end


####TEST:
## slen = 1000
## input = randstring("ACGT",slen)

## DNA_seq = string_to_DNASeq(input)
## k=3
## kmer_list = read_kmer(DNA_seq, length(input),k)

## G = graph_creator(kmer_list,['A','C','G','T'], 5)
## contig_list = []
## G_new,ls = compact_graph!(G,3)
## output = run_walk(G, ls)




"""
slen = 10
input = randstring("ACGT",slen)

DNA_seq = string_to_DNASeq(input)
k=3
kmer_list = read_kmer(DNA_seq, length(input),k)

G = graph_creator(kmer_list,['A','C','G','T'], 5)
contig_list = []
G_new,ls = compact_graph!(G,100)
output = run_walk(G, [])

for i in contig_list
    for j in i
        print(DNASeq_to_string(j))
    end
    print("\n")
end

for i in input
    print(i)
end

"""
