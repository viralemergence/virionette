import Pkg;
Pkg.activate(".");

using DataFrames
import CSV
using GBIF

#--- Load the ICTV master data

ictv_path = joinpath(pwd(), "data", "scaffold", "ictv.csv")
ictv =
    CSV.read(ictv_path, types = [String, String, String, Union{String,Missing}])

#--- Load the dataframes templates and other functions

include(joinpath(pwd(), "data", "scaffold", "lib", "dataframes.jl"))
include(joinpath(pwd(), "data", "scaffold", "lib", "methods.jl"))

#--- Specify the paths and load the files

hp3_path = joinpath("data", "raw", "HP3")
hp3_files = ["associations", "hosts", "viruses"]
hp3_assoc, hp3_hosts, hp3_viruses = [
    CSV.read(
        joinpath(hp3_path, "$(hp3f).csv"),
        missingstrings = ["NA", "Unassigned"],
    ) for hp3f in hp3_files
]

#--- Clean host file

clean_host_file = CSV.read(joinpath(hp3_path, "HP3-associations_all mammal_compatible.csv"))
select!(clean_host_file, Not(:Column1))
hp3_assoc = clean_host_file

#--- Map the viruses

hp3_entities_virus = entities_scaffold()

# Genus, Subfamily, Family, Order

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

hp3_host_taxa = host_scaffold()
hp3_entities_host = entities_scaffold()

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

hp3_associations = associations_scaffold()

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

#--- Make the path

data_path = joinpath("data", "scaffold", "HP3")
ispath(data_path) || mkdir(data_path)

#--- Write the files

CSV.write(joinpath(data_path, "entities_virus.csv"), unique(hp3_entities_virus))
CSV.write(joinpath(data_path, "entities_host.csv"), unique(hp3_entities_host))
CSV.write(joinpath(data_path, "taxonomy_host.csv"), unique(hp3_host_taxa))
CSV.write(joinpath(data_path, "associations.csv"), unique(hp3_associations))
