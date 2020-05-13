using CSV
using DataFrames

#--- Prepare the folders

scaffold_path = joinpath("01_scaffolded_data")
flat_path = joinpath("02_flat_files")

ispath(flat_path) || mkdir(flat_path)

#--- Copy the ICTV file

cp(joinpath(scaffold_path, "ictv.csv"), joinpath(flat_path, "virus.csv"))

#--- Type annotations for the various files

taxo_types = [String, fill(Union{String,Missing}, 7)..., fill(Union{Int64,Missing}, 7)...]
entity_types = [String, String, String, Union{String,Missing}]
associations_types = [String, String, String, String, String, String]

#--- Read and merge the host taxonomies

genbank_host_taxo = CSV.read(joinpath(scaffold_path, "Genbank", "taxonomy_host.csv"); types=taxo_types)
hp3_host_taxo = CSV.read(joinpath(scaffold_path, "HP3", "taxonomy_host.csv"); types=taxo_types)

flat_host_taxo = unique(vcat(genbank_host_taxo, hp3_host_taxo))
sort!(flat_host_taxo, :genus)

CSV.write(joinpath(flat_path, "hosts.csv"), flat_host_taxo)

## Read and merge the entitites (hosts)
genbank_host_entities = CSV.read(joinpath(scaffold_path, "Genbank", "entities_host.csv"); types=entity_types)
hp3_host_entities = CSV.read(joinpath(scaffold_path, "HP3", "entities_host.csv"); types=entity_types)

flat_host_entities = unique(vcat(genbank_host_entities, hp3_host_entities))
sort!(flat_host_entities, :origin)

CSV.write(joinpath(flat_path, "entities_hosts.csv"), flat_host_entities)

#--- Read and merge the entitites (viruses)

genbank_virus_entities = CSV.read(joinpath(scaffold_path, "Genbank", "entities_virus.csv"); types=entity_types)
hp3_virus_entities = CSV.read(joinpath(scaffold_path, "HP3", "entities_virus.csv"); types=entity_types)

flat_virus_entities = unique(vcat(genbank_virus_entities, hp3_virus_entities))
sort!(flat_virus_entities, :name)

CSV.write(joinpath(flat_path, "entities_virus.csv"), flat_virus_entities)

#--- Read and merge the associations files

genbank_associations = CSV.read(joinpath(scaffold_path, "Genbank", "associations.csv"); types=associations_types)
hp3_associations = CSV.read(joinpath(scaffold_path, "HP3", "associations.csv"); types=associations_types)

flat_associations = unique(vcat(genbank_associations, hp3_associations))
sort!(flat_associations, :virus)

CSV.write(joinpath(flat_path, "associations.csv"), flat_associations)
