using CSV
using DataFrames

#--- Prepare the folders

data_path = joinpath("02_flat_files")

#--- Type annotations for the various files

taxo_types = [String, fill(Union{String,Missing}, 7)..., fill(Union{Int64,Missing}, 7)...]
entity_types = [String, String, String, Union{String,Missing}]
associations_types = [String, String, String, String, String, String]
ictv_types = [String, String, String, Union{String,Missing}]

#--- Load the files

hosts = CSV.read(joinpath(data_path, "hosts.csv"); types=taxo_types)
viruses = CSV.read(joinpath(data_path, "virus.csv"); types=ictv_types)
associations = CSV.read(joinpath(data_path, "associations.csv"); types=associations_types)
hosts_entities = CSV.read(joinpath(data_path, "entities_hosts.csv"); types=entity_types)
viruses_entities = CSV.read(joinpath(data_path, "entities_virus.csv"); types=entity_types)

#--- Merge the host entities

hosts_merged = join(hosts_entities, hosts; on=:match => :id)
rename!(hosts_merged, :name => :host_entity_name)
rename!(hosts_merged, :match => :host_match)
rename!(hosts_merged, :id => :host_entity_id)

#--- Merge the virus entities

viruses_merged = join(viruses_entities, viruses; on=:match => :id, makeunique=true)
rename!(viruses_merged, :id => :virus_entity_id)
rename!(viruses_merged, :name => :virus_entity_name)
rename!(viruses_merged, :match => :virus_match)
rename!(viruses_merged, :name_1 => :virus_name)

#--- Merge the associations

#rename!(associations, :id => :association_id)
associations_with_viruses = join(viruses_merged, associations; on=[:virus_entity_id=>:virus, :origin=>:source], makeunique=true)
associations_merged = join(associations_with_viruses, hosts_merged; on=[:origin, :host=>:host_entity_id], makeunique=true)

## Get the whole thing
compl = associations_merged[.!ismissing.(associations_merged.species),:]
select!(compl, Not(r"_id"))
select!(compl, Not(r"_entity_name"))
select!(compl, Not(r"_match"))
for c in [:method, :index, :kingdom, :phylum, :host]
    select!(compl, Not(c))
end
compl = compl[(compl.rank.=="Genus").|(compl.rank.=="Subgenus"),:]

vfam = Union{Missing,String}[]
vgen = Union{Missing,String}[]
vsubgen = Union{Missing,String}[]

# NOTE this is ugly and really fragile, and should instead use something like
# `ancestor_of`, and do some light recursion, but...
function family_of(genus, viruses)
    v = first(viruses[viruses.name .== genus, :])
    @assert v.rank == "Genus"
    a = viruses[viruses.id .== v.ancestor, :]
    if a.rank == "Family"
        return a.name
    end
    anc = first(a.ancestor)
    if !ismissing(anc)
        b = first(viruses[viruses.id .== anc, :])
        if b.rank == "Family"
            return b.name
        end
    end
    return missing
end

for row in eachrow(compl)
    if row.rank == "Genus"
        push!(vgen, row.virus_name)
        push!(vsubgen, missing)
    else
        push!(vgen, viruses.name[findfirst(viruses.id .== row.ancestor)])
        push!(vsubgen, row.virus_name)
    end
    vgen_name = row.rank == "Genus" ? row.virus_name : viruses.name[findfirst(viruses.id .== row.ancestor)]
    push!(vfam, family_of(vgen_name, viruses))

end

select!(compl, Not(:virus_name))
select!(compl, Not(:rank))
select!(compl, Not(:ancestor))
rename!(compl, :family => :host_family)
rename!(compl, :genus => :host_genus)
rename!(compl, :class => :host_class)
rename!(compl, :order => :host_order)
rename!(compl, :species => :host_species)
compl.virus_family = vfam
compl.virus_genus = vgen
compl.virus_subgenus = vsubgen

sort!(compl, :virus_family)

#--- Correct the taxonomic names

corrections = CSV.read(joinpath("00_raw_data", "CORRECTIONS", "mammal_tree_name_matching.csv"))

for namechange in eachrow(corrections)
    old = replace(namechange.old_name, "_" => " ")
    new = replace(namechange.tree_name, "_" => " ")
    @info "Replaced $(old) by $(new) $(sum(compl.host_species.==old)) times"
    compl.host_species[compl.host_species.==old] .= new
end

compl = unique(compl)
select!(compl, Not(:id))
select!(compl, Not(:virus_subgenus))

compl = unique(compl)

#--- Do the last bit of cleaning

# Only mammals
virionette = compl[compl.host_class.=="Mammalia",:]

# No human
virionette = virionette[virionette.host_species.!="Homo sapiens",:]

# No Genbank toroviruses
virionette = virionette[.!((virionette.origin.=="Genbank").&(virionette.virus_genus.=="Torovirus")),:]

#--- Write the final flat file

net_path = joinpath("03_interaction_data")
ispath(net_path) || mkdir(net_path)
sort!(virionette, [:virus_family, :virus_genus, :host_species])
CSV.write(joinpath(net_path, "virionette.csv"), virionette)
