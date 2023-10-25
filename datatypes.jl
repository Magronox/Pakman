#### Adjusted Pakman -- removing affixes, fixing macro nodes

using StaticArrays, DataStructures, Random

const k = 32 ### the code is optimizied for k = 32, otherwise it needs to be adjusted.
const alphabet = ['A','C','G','T']


Base.@kwdef mutable struct DNASeq
    bit1 = MVector{64, Bool}
    bit2 = MVector{64, Bool}
    len = Int64
    function DNASeq(bit1 ::Vector, bit2 :: Vector, len :: Int64)
        Bit1 = MVector{64, Bool}(zeros(Bool, 64))
        Bit2 = MVector{64, Bool}(zeros(Bool, 64))
        Bit1[end-len+1 : end] = bit1[end-len+1 : end]
        Bit2[end-len+1 : end] = bit2[end-len+1 : end]
        new(Bit1, Bit2, len)
    end
end

Base.@kwdef mutable struct mn_label
    bit1 = MVector{31, Bool}(zeros(Bool, 31))
    bit2 = MVector{31, Bool}(zeros(Bool, 31))
    isTerminal = false
end

Base.@kwdef mutable struct macro_node
    label = mn_label()
    prefix_edge_ids = Set{Int64}()
    suffix_edge_ids = Set{Int64}()
    macro_node(label,prefix_edge_ids,suffix_edge_ids) = new(label, prefix_edge_ids, suffix_edge_ids)
end

#Base.@kwdef mutable struct macro_node
    #length k-1
    #id = 0
#    label :: mn_label#BitArray{1}(undef,64)
    #prefixes = DefaultDict{Int64, Vector{DNASeq}}(DNASeq[])
    #prefix_counts = Dict{Int64PπP, Tuple}()
    #prefixes_terminal = DefaultDict{Int64,Bool}(false)
    #suffixes = DefaultDict{Int64, Vector{DNASeq}}(DNASeq[])
    #suffix_counts = Dict{Int64,Tuple}()
    #suffixes_terminal = DefaultDict{Int64,Bool}(false)
    #wire_info = DefaultDict{Int64, Set{Tuple}}(Set())
    #prefix_begin_info = DefaultDict{Int64,Tuple}((-1,-1))
#end

Base.@kwdef mutable struct edge
    edge_id = -1
    pred_label = mn_label()
    pred_suffix = Vector{DNASeq}()
    succ_label = mn_label()
    succ_prefix = Vector{DNASeq}()
    counts = Dict{Int64, Tuple}()
    succ_counts = Dict{Int64, Tuple}()
    isPredTerminal = falses
    isSuccTerminal = false
    wire_info = DefaultDict{Int64, Set{Tuple}}(Set())
    prefix_begin_info = DefaultDict{Int64, Tuple}((-1,-1))
end


function mn_label_(seq :: DNASeq)
    var = mn_label()
    if seq.len == 31
        print("ll\n",MVector{31}(seq.bit1[end-k+2:end]), length(seq.bit1))
        var.bit1 = MVector{31}(seq.bit1[end-k+2:end])
        var.bit2 = MVector{31}(seq.bit2[end-k+2:end])
        var.isTerminal = false
    elseif seq.len == 0
        var.isTerminal = true
    else
        @assert(false)
    end
    var
end
"""
function mn_label(bit1 :: BitArray, bit2 :: BitArray, isTerminal :: Bool )
    var = mn_label()
    if ((length(bit2)==length(bit1))&& length(bit1) == k-1) 
        var.bit1 = bit1
        var.bit2 = bit2
        var.isTerminal = isTerminal
    else
         @assert(false)
    end
    var
end"""