using DataFrames
import CSV
using GBIF

#--- Read the raw genbank data

gb_data_path = joinpath("00_raw_data", "Genbank", "GB_CoV_VRL_noSeqs.csv")
gb_raw = CSV.read(gb_data_path; missingstring = "NA")

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

#--- Prepare the empty dataframes

gb_entities_virus = entities_scaffold()
gb_entities_host = entities_scaffold()
gb_hosts = host_scaffold()
gb_associations = associations_scaffold()

#--- GET. THEM. VIRUSES.

# This is a regex to pick the correct name of the virus from the genbank
# columns. It may not be perfect, but it matches on all rows so far.
match_virus = r"(\w+)vir(us|idae|inae|ales)"

# It doesn't really make sense to do the matching on every row, so we will focus
# only on the unique names.
unique_virus_names = unique(gb_raw.gbGenus)

# We then fill in the viral entity dataframe
for virus_name in unique_virus_names

    # First, we look for the virus name
    regex_match = match(match_virus, virus_name)

    # If it doesn't match, we add an entity without a taxonomy
    if isnothing(regex_match)
        push!(
            gb_entities_virus,
            (
                string(hash(virus_name * "Genbank")),
                "Genbank",
                virus_name,
                missing,
            ),
        )
    else
        # If it matches, we have two choices
        matching_idx =
            findfirst(lowercase.(ictv.name) .== lowercase(regex_match.match))

        # If the match is not in ICTV, we add an entity without a taxonomy
        if isnothing(matching_idx)
            push!(
                gb_entities_virus,
                (
                    string(hash(virus_name * "Genbank")),
                    "Genbank",
                    virus_name,
                    missing,
                ),
            )
        else
            # But if it does, we add the entity with a taxonomy
            push!(
                gb_entities_virus,
                (
                    string(hash(virus_name * "Genbank")),
                    "Genbank",
                    virus_name,
                    ictv.id[matching_idx],
                ),
            )
        end
    end
end

#--- Work on the hosts using GBIF

unique_gb_hosts = filter(!ismissing, unique(gb_raw.gbHost))

# This part is potentially the longest, as it will do a GBIF query on all the
# hosts in the genbank data - this is at most a few minutes, so no stress
for (i, host) in enumerate(unique_gb_hosts)
    entity_hash = string(hash(host * "Genbank"))
    match_gbif = missing

    # The first step is to try and get a GBIF match
    try
        match_gbif = taxon(host, strict = false)
    catch
        # If it does not match, we push an entity with no taxonomy
        push!(gb_entities_host, (entity_hash, "Genbank", host, missing))
        continue
    end

    # If it does match, the loop continues, and so we will add something to the
    # entity list - this part is creating a hash for the taxonomic entity.
    match_gbif_hash = ismissing(match_gbif) ? missing : string(hash(match_gbif))

    # If there is a match on GBIF, we add an entry to the host taxonomy - this
    # is not the entity file, this is the actual taxonomy with the full
    # hierarchy, based on the taxonomic backbone.
    if !isnothing(match_gbif)
        if !(match_gbif_hash in gb_hosts.id)
            push!(gb_hosts, (match_gbif_hash, convert(Tuple, match_gbif)...))
        end
    end

    # Finally we add to the host entity match file
    push!(gb_entities_host, (entity_hash, "Genbank", host, match_gbif_hash))
end

#--- Create the associations dataframe

# This loop is rather straightforward - we hash everything, and add to the
# association table
for row in eachrow(gb_raw)
    association_id = string(hash(row))
    host_hash = string(hash(row.gbHost * "Genbank"))
    virus_hash = string(hash(row.gbGenus * "Genbank"))
    source = "Genbank"
    index = row.gbAccession # This is the unique ID of the SEQUENCE!
    push!(
        gb_associations,
        (association_id, host_hash, virus_hash, source, index, "GenBank"),
    )
end

#--- Write the file

data_path = joinpath("01_scaffolded_data", "Genbank")
ispath(data_path) || mkdir(data_path)
CSV.write(joinpath(data_path, "entities_virus.csv"), unique(gb_entities_virus))
CSV.write(joinpath(data_path, "entities_host.csv"), unique(gb_entities_host))
CSV.write(joinpath(data_path, "taxonomy_host.csv"), unique(gb_hosts))
CSV.write(joinpath(data_path, "associations.csv"), unique(gb_associations))
