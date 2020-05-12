using DataFrames
import CSV
using GBIF

#--- We need the formatted ICTV file before we start working on this

# We are going to help the CSV parser a little bit by giving it information
# about the column types.
ictv_path = joinpath("01_scaffolded_data", "ictv.csv")
ictv =
    CSV.read(ictv_path, types = [String, String, String, Union{String,Missing}])

#--- Load some extra helper functions

# These are empty dataframes
include(joinpath(pwd(), "lib", "dataframes.jl"))
# These are GBIF helpers
include(joinpath(pwd(), "lib", "methods.jl"))

#--- Specify the paths and load the files

hp3_path = joinpath("00_raw_data", "HP3")
hp3_files = ["associations", "hosts", "viruses"]

# We are going to load everything in one go, and the only trick here is that HP3
# uses both NA and Unassigned as markers for lack of information. These probably
# have different semantics, but for our purpose these correspond to Missing.
hp3_assoc, hp3_hosts, hp3_viruses = [
    CSV.read(
        joinpath(hp3_path, "$(hp3f).csv"),
        missingstrings = ["NA", "Unassigned"],
    ) for hp3f in hp3_files
]

#--- Prepare some empty dataframes

hp3_entities_virus = entities_scaffold()
hp3_entities_host = entities_scaffold()
hp3_host_taxa = host_scaffold()
hp3_associations = associations_scaffold()

#--- Map the viruses

# This part follows the same logic as the genbank file, so it's a little more
# lightly commented. The only difference is that we pick the most resolved level
# for every row in the HP3 dataset.
for virus in eachrow(hp3_viruses)
    id = string(hash(virus.vVirusNameCorrected * "HP3"))
    origin = "HP3"
    name = virus.vVirusNameCorrected
    match = missing
    if !ismissing(virus.vGenus)
        if virus.vGenus in ictv.name
            match = ictv.id[findfirst(ictv.name .== virus.vGenus)]
            push!(hp3_entities_virus, (id, origin, name, match))
            continue
        end
    end
    if !ismissing(virus.vSubfamily)
        if virus.vSubfamily in ictv.name
            match = ictv.id[findfirst(ictv.name .== virus.vSubfamily)]
            push!(hp3_entities_virus, (id, origin, name, match))
            continue
        end
    end
    if !ismissing(virus.vFamily)
        if virus.vFamily in ictv.name
            match = ictv.id[findfirst(ictv.name .== virus.vFamily)]
            push!(hp3_entities_virus, (id, origin, name, match))
            continue
        end
    end
    if !ismissing(virus.vOrder)
        if virus.vOrder in ictv.name
            match = ictv.id[findfirst(ictv.name .== virus.vOrder)]
            push!(hp3_entities_virus, (id, origin, name, match))
            continue
        end
    end
end

#--- Populate the tables with host information

# Same as above, this is very close to the genbank script
for (idx, host_row) in enumerate(eachrow(hp3_hosts))
    host_name = "$(host_row.hGenus) $(host_row.hSpecies)"
    # Prepare the entity match
    entity_hash = string(hash(host_row.hHostNameFinal * "HP3"))
    entity_row = idx
    match_gbif = missing
    try
        match_gbif = taxon(host_name, strict = false)
    catch
        continue
    end
    # Add to the taxonomy table
    match_gbif_hash = ismissing(match_gbif) ? missing : string(hash(match_gbif))
    if !isnothing(match_gbif)
        push!(hp3_host_taxa, (match_gbif_hash, convert(Tuple, match_gbif)...))
    end
    # Return everything
    if !(host_row.hHostNameFinal in hp3_entities_host.name)
        push!(
            hp3_entities_host,
            (entity_hash, "HP3", host_row.hHostNameFinal, match_gbif_hash),
        )
    end
end

#--- Prepare an association table

# Again - this is a simple mapping - and this is relatively quick
for (i, row) in enumerate(eachrow(hp3_assoc))
    association_id = string(hash(row))
    host_id = string(hash(row.hHostNameFinal * "HP3"))
    virus_id = string(hash(row.vVirusNameCorrected * "HP3"))
    source = "HP3"
    index = i
    method = replace(replace(row.DetectionMethod, ", " => "_"), " " => "_")
    push!(
        hp3_associations,
        (association_id, host_id, virus_id, source, index, method),
    )
end

#--- Make the path and write them

data_path = joinpath("01_scaffolded_data", "HP3")
ispath(data_path) || mkdir(data_path)

CSV.write(joinpath(data_path, "entities_virus.csv"), unique(hp3_entities_virus))
CSV.write(joinpath(data_path, "entities_host.csv"), unique(hp3_entities_host))
CSV.write(joinpath(data_path, "taxonomy_host.csv"), unique(hp3_host_taxa))
CSV.write(joinpath(data_path, "associations.csv"), unique(hp3_associations))
