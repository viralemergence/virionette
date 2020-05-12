import Pkg; Pkg.activate(".")

using DataFrames
import CSV
using GBIF

## Read the raw data
anth_data_path = joinpath(pwd(), "data", "raw", "Anthony", "GB_CoV_VRL_noSeqs2.csv")
anth_raw = CSV.read(anth_data_path; missingstring="NA")

## Load the ICTV master data
ictv_path = joinpath(pwd(), "data", "scaffold", "ictv.csv")
ictv = CSV.read(ictv_path, types=[String, String, String, Union{String,Missing}])

## Load the dataframes templates and other functions
include(joinpath(pwd(), "data", "scaffold", "lib", "dataframes.jl"))
include(joinpath(pwd(), "data", "scaffold", "lib", "methods.jl"))

## Map the viruses
anthony_entities_virus = entities_scaffold()

## Virus mapping
match_virus = r"(\w+)vir(us|idae|inae|ales)"
unique_virus_names = unique(anth_raw.gbGenus)

for virus_name in unique_virus_names
    regex_match = match(match_virus, virus_name)
    
    if isnothing(regex_match)
        push!(anthony_entities_virus, (
            string(hash(virus_name*"Anthony")), "Anthony", virus_name, missing
        ))
    else
        matching_idx = findfirst(lowercase.(ictv.name) .== lowercase(regex_match.match))
        if isnothing(matching_idx)
            push!(anthony_entities_virus, (
                string(hash(virus_name*"Anthony")), "Anthony", virus_name, missing
            ))
        else
            push!(anthony_entities_virus, (
                string(hash(virus_name*"Anthony")), "Anthony", virus_name, ictv.id[matching_idx]
            ))
        end
    end
end

## Hosts
anthony_hosts = host_scaffold()
anthony_entities_host = entities_scaffold()

unique_anth_hosts = filter(!ismissing, unique(anth_raw.gbHost))


for (i, host) in enumerate(unique_anth_hosts)
    entity_hash = string(hash(host*"Anthony"))
    match_gbif = missing
    try
        match_gbif = taxon(host, strict=false)
    catch
        push!(anthony_entities_host, (entity_hash, "Anthony", host, missing))
        continue
    end
    # Add to the taxonomy table
    match_gbif_hash = ismissing(match_gbif) ? missing : string(hash(match_gbif))
    if !isnothing(match_gbif)
        if !(match_gbif_hash in anthony_hosts.id)
            push!(anthony_hosts,
                (match_gbif_hash, convert(Tuple, match_gbif)...)
            )
        end
    end
    # Add to the entity table
    push!(anthony_entities_host, (entity_hash, "Anthony", host, match_gbif_hash))
end

## Associations
anthony_associations = associations_scaffold()

for row in eachrow(anth_raw)
    association_id = string(hash(row))
    host_hash = string(hash(row.gbHost*"Anthony"))
    virus_hash = string(hash(row.gbGenus*"Anthony"))
    source = "Anthony"
    index = row.gbAccession
    push!(anthony_associations,
    (
        association_id,
        host_hash,
        virus_hash,
        source,
        index,
        "GenBank"
    ))
end

## Make the path
data_path = joinpath("data", "scaffold", "Anthony")
ispath(data_path) || mkdir(data_path)

## Write the files
CSV.write(joinpath(data_path, "entities_virus.csv"), unique(anthony_entities_virus))
CSV.write(joinpath(data_path, "entities_host.csv"), unique(anthony_entities_host))
CSV.write(joinpath(data_path, "taxonomy_host.csv"), unique(anthony_hosts))
CSV.write(joinpath(data_path, "associations.csv"), unique(anthony_associations))
