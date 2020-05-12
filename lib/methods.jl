function Base.convert(::Type{Tuple}, tax::GBIFTaxon)
    txt = Union{String,Missing}[]
    idx = Union{Integer,Missing}[]
    for l in [:kingdom, :phylum, :class, :order, :family, :genus, :species]
        if !ismissing(getfield(tax, l))
            push!(txt, getfield(tax, l).first)
            push!(idx, getfield(tax, l).second)
        else
            push!(txt, missing)
            push!(idx, missing)
        end
    end
    return (txt..., idx...)
end
