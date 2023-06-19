

function compact_graph( G :: DefaultDict{DNASeq,macro_node}, phi :: Int64)
    num_mn = length(G)
    while(num_mn > phi)
        new_G, num_mn = IS(G)

    end
end

function IS(G)
    ## Create independent set of G


end