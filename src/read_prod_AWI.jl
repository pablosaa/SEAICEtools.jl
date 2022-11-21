# File part of SEAICEtools.jl
# this file include function to read the Lead Fraction (LF) products provided
# by Luisa von Albedyl from AWI
#
# See LICENSE


# --------------------------------------------------------------------------------
# Function to load the latitude and longitude grid data for Sea Ice products by
# the University of Bremen.
module load

using NCDatasets
using Dates

"""
Function to find file that follows pattern given by PATH and year, month, day

julia> file_name = FilePattern("/data/SAR", 2020, 4, 15)

julia> file_name = FilePattern("/data/SAR", 2020, 4, 15, moved=false)

The first case returns:

the second case returns:

"""
function FilePattern(PATH::String, PROD::String, yy::Number, mm::Number, dd::Number; moved=true)

    DIRSAR = joinpath(PATH, PROD, "$(yy)")

    !isdir(DIRSAR) && @error "I could not find directory $(DIRSAR)"
    tmp = readdir(DIRSAR, join=true)
    
    today = Date(yy, mm, dd) |> x->Dates.format(x, "yyyymmddT")
    intmp = contains.(tmp, today) |> findall
    tmp = tmp[intmp]

    with_pairs = contains.(tmp, "pairs_"*today)
    intmp = if moved
        .!with_pairs |> findfirst
    else
        with_pairs |> findfirst
    end
    return !isnothing(intmp) ? tmp[intmp] : intmp
end

function LatLon_AWI_product(filen::String)

    Missing2NaN(X) = replace(X, missing=>NaN)
    
    nc = NCDataset(filen, "r")
    xx = Missing2NaN(nc["x"][:,:])
    yy = Missing2NaN(nc["y"][:,:])
    close(nc)
    
    return xx, yy
end

function Data_Divergence_LeadFraction(filen::String, idx_sector::Vector{CartesianIndex{2}}; moved=true)
    vars = if moved
       Dict(:DIV=>"filtered_div_moved",
            :LF=>"lead_fraction_moved")
    else
        Dict(:DIV=>"filtered_divergence",
             :LF=>"lead_fraction")
    end

    #vars = merge(vars, Dict(:x=>"x", :y=>"y"))
    
    Missing2NaN(X) = replace(X, missing=>NaN)
    
    sar = NCDataset(filen, "r") do nc
        Dict(k => Missing2NaN(nc[v][idx_sector]) for (k,v) in vars)
    end

    return sar
end

end # end module "load"

# end of file
