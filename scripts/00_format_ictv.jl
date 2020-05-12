using DataFrames
import CSV

#--- Read the ICTV file with all columns

# The ICTV master file is stored as a csv, but it has a bunch of columns we
# don't really need.
ictv_master = CSV.read(joinpath(pwd(), "00_raw_data", "ictv_master.csv"))

#--- Drop the columns we don't need

# Some of these columns are internal metadata, or taxonomic levels we do not
# need, so for the sake of simplicity in the flat files, we will remove them.
# This bit is a list of columns to remove.
ictv_columns_to_drop = [
    Symbol("Sort"),
    Symbol("Species"),
    Symbol("Type Species?"),
    Symbol("Genome Composition"),
    Symbol("Last Change"),
    Symbol("MSL of Last Change"),
    Symbol("Proposal for Last Change "),
    Symbol("Taxon History URL"),
    :Column23,
]

# And then we remove all of these columns one by one
for column_to_drop in ictv_columns_to_drop
    select!(ictv_master, Not(column_to_drop))
end

#--- Create the new dataframe

# We create an empty data frame, which has the unique identifier of the virus,
# its name, its taxonomic rank, and then its ancestor. The ancestor column has
# the identifier of the previous taxonomic level - having an ancestor that is
# Missing means that this is either a top-level taxon, or a taxon that is sort
# of floating in there.
ictv = DataFrame(
    id = Symbol[],
    name = String[],
    rank = Symbol[],
    ancestor = Union{Symbol,Missing}[],
)

#--- Convert the ICTV database

# We start by listing all unique taxa, by only keeping unique rows in the ICTV
# master file. This is required because we have dropped a few columns, for
# example viral species, and all metadata, so there are some pseudo-duplicates
# now.
ictv_taxa = unique(ictv_master)

# All the column names that remain are taxonomic ranks - so we can get the list
# of ranks by taking the column names.
all_ranks = names(ictv_master)

# This next bit will go through the entire list of unique taxa, and add them to
# our own dataframe.
for taxon in eachrow(ictv_taxa)

    # This next block is subsetting the ranks that are not missing (i.e. ICTV
    # knows what they are), and gets their values.
    taxon_values = values(taxon)
    no_missing = findall(!ismissing, taxon_values)
    this_ranks = all_ranks[no_missing]
    this_names = taxon_values[no_missing]

    # Then for every rank/name combination...
    for (depth, name) in enumerate(this_names)

        # ... we get its ancestor
        ancestor = missing
        if depth >= 2
            ancestor =
                Symbol(hash(this_names[depth-1] * string(this_ranks[depth-1])))
        end

        # ... then we give it a unique identifier
        name_hash = Symbol(hash(name * string(this_ranks[depth])))

        # ... and we add it to our dataframe only if it is not already there.
        if !(name_hash in ictv.id)
            push!(ictv, (name_hash, name, this_ranks[depth], ancestor))
        end
    end
end

# Let's rank by name just because
sort!(ictv, :name)

#--- Make the path and write the file

data_path = joinpath("01_scaffolded_data")
ispath(data_path) || mkdir(data_path)
CSV.write(joinpath(data_path, "ictv.csv"), unique(ictv))
