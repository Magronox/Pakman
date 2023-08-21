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
                break
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
contig_list = []
function walk!(G:: DefaultDict, pcontig:: Vector{DNASeq}, freq, offset_in_prefix :: Int64 , mn :: macro_node, pid :: Int64, output :: Set)
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

    for (id,offset_in_suffix,count) in mn.wire_info[pid]
       
        internal_off += mn.prefix_begin_info[pid][2]
        sz = copy(count)
        
        if internal_off + sz <= offset_in_prefix || internal_off > offset_in_prefix + freq_rem
            continue
        end
        off_in_wire = offset_in_prefix <= internal_off ? 0 : offset_in_prefix - internal_off 
        next_off = offset_in_suffix + id
        freq_in_wire = min(freq_rem, sz - off_in_wire )
        contig = kmerge(pcontig, mn.suffixes[id])
        #print(@which kmerge(pcontig, mn.suffixes[id]))
        #assert(false)
        #print(DNASeq_to_string(contig[en]),"\n\n")
        if pcontig in output
            continue
        end
        if mn.suffixes_terminal[id]
            print("contig found\n")
            push!(contig_list, pcontig)
            push!(output, pcontig)
            continue
        else
            succ_ext = mn.suffixes[id]
            key = mn.label
            succ_node = succ_neigh(key[1],succ_ext)
            #ppid, succ_ext = find_succ_ext(G,key[1], succ_node[1])
            if !(succ_node in keys(G))
            #    print("Error key not found\n")
            continue
            end
            next_mn = G[succ_neigh(mn.label[1],mn.suffixes[id])]
            next_prefix_id, ~ = find_succ_ext(G, mn.label, next_mn.label)
            @assert(succ_node == next_mn.label)
            if freq_in_wire > 0 && !(pcontig in output)
                walk!(G, contig, freq_in_wire, next_off, next_mn, next_prefix_id, output);
            end 
               
            
        end
       
        freq_rem -= freq_in_wire
 


    end
end


function identify_begin_kmers(G, BeginMN)

    for (i,j) in G
        for (ii,jj) in j.prefixes
            if j.prefix_counts[ii]>0 && j.prefixes_terminal[ii]
                push!(begin_kmer_list, kmerge())
            end
        end
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



""" contig list is defined globally, so you need to initialize it as an empty vector before running the algorithm. 
This is because of the recursive nature of the walk algorithm."""

"""
slen = 5000
input = randstring("ACGT",slen)

DNA_seq = string_to_DNASeq(input)
k=10
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