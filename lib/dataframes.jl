entities_scaffold() = DataFrame(
    id = String[],
    origin = String[], 
    name = String[],
    match = Union{String,Missing}[]
)

host_scaffold() = DataFrame(
    id = String[],
    kingdom = Union{String,Missing}[],
    phylum = Union{String,Missing}[],
    class = Union{String,Missing}[],
    order = Union{String,Missing}[],
    family = Union{String,Missing}[],
    genus = Union{String,Missing}[],
    species = Union{String,Missing}[],
    kingdom_id = Union{Integer,Missing}[],
    phylum_id = Union{Integer,Missing}[],
    class_id = Union{Integer,Missing}[],
    order_id = Union{Integer,Missing}[],
    family_id = Union{Integer,Missing}[],
    genus_id = Union{Integer,Missing}[],
    species_id = Union{Integer,Missing}[]
)

associations_scaffold() = DataFrame(
    id = String[],
    host = String[],
    virus = String[],
    source = String[],
    index = Any[],
    method = Union{String,Missing}[]
)