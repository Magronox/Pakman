include("graph_compaction.jl")



function run_walk(G :: DefaultDict)
    output = Set()
    begin_kmer_list = copy(keys(G))
    for mn in begin_kmer_list
        node = G[mn]
        pid = node.prefix_terminal_id
        if pid != -1
            freq = node.prefix_counts[pid][2]
            contig = node.prefixes[pid]
            contig = kmerge(contig, mn)
            
            walk!(G, contig, freq, 0, node, pid, output)
        end

    end
    return output
end

function walk!(G:: DefaultDict, contig:: Vector{DNASeq}, freq, offset_in_prefix :: Int64 , mn :: macro_node, pid :: Int64, output :: Set)
    #pc_size = (length(contig)-1)*64 + contig[1].len
    freq_rem = freq
    internal_off = 0
    off_in_wire = -1
    for (sid,offset_in_suffix,sz) in mn.wire_info.prefix_info[pid]
        
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



####TEST:
## slen = 1000
## input = randstring("ACGT",slen)

## DNA_seq = string_to_DNASeq(input)
## k=3
## kmer_list = read_kmer(DNA_seq, length(input),k)

## G = graph_creator(kmer_list,['A','C','G','T'], 5)
## I = IS(G,['A','C','G','T'],k)
## G_new, ~ = compact_graph(G,k,5)
## output = run_walk(G)