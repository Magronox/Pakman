include("graph_compaction.jl")



function run_walk(G :: DefaultDict, pcontig_list:: Vector)
    output = Set()
    begin_kmer_list = copy(keys(G))
    pid = -1
    for node in begin_kmer_list
        mn = G[node]
        for i in 1:length(mn.prefixes_terminal)
            if mn.prefixes_terminal[i] == true
                pid = i
            end
        end

        if pid != -1
            freq = mn.prefix_counts[pid][2]
            contig = mn.prefixes[pid]
            contig = kmerge(contig, node)
            
            walk!(G, contig, freq, 0, mn, pid, output)
        end

    end
    return output
end

function walk!(G:: DefaultDict, pcontig:: Vector{DNASeq}, freq, offset_in_prefix :: Int64 , mn :: macro_node, pid :: Int64, output :: Set)
    #pc_size = (length(contig)-1)*64 + contig[1].len
    freq_rem = freq
    internal_off = 0
    off_in_wire = -1


    for (id,offset_in_suffix,count) in mn.wire_info[pid]
        internal_off += count
        sz = copy(count)
        if internal_off + sz <= offset_in_prefix || internal_off > offset_in_prefix + freq_rem
            continue
        end
        off_in_wire = offset_in_prefix <= internal_off ? 0 : offset_in_prefix - internal_off
        next_off = offset_in_suffix + off_in_wire
        freq_in_wire = min(freq_rem, sz - off_in_wire)
        contig = kmerge(pcontig, mn.suffixes[id])
        if mn.suffixes_terminal[id]
            push!(contig_list, contig)
            push!(output, contig)
        else
            succ_ext = mn.suffixes[id]
            key = mn.label
            succ_node = succ_neigh(key[1],succ_ext)
            #ppid, succ_ext = find_succ_ext(G,key[1], succ_node[1])
            if !(succ_node in keys(G))
                print("Error key not found\n")
            end
        end
        next_mn = G[succ_neigh(mn.label[1],mn.suffixes[id])]
        next_prefix_id, ~ = find_succ_ext(G, mn.label, next_mn.label)
        @assert(succ_node == next_mn)

        walk(G, pcontig, freq_in_wire, next_off, next_prefix_id, next_mn, output);
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
## I = IS(G,['A','C','G','T'],k)
## G_new, ~ = compact_graph(G,k,5)
## output = run_walk(G)
