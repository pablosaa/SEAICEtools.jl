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
```julia-repl
julia> file_name = FilePattern(PATH_DATA, PROD_PATH, yy, mm, dd)
julia> file_name = FilePattern("/data/SAR", 2020, 4, 15, moved=false)
julia> file_name = FilePattern("/data/SAR", 2020, 4, 15, Format4date="ddmmyyyy")
```

WHERE:
* PATH\\_DATA::String indicating base path to the data location,
* PROD\\_PATH::String indicating the path or name to the product to locate,
* yy::Int year on data file,
* mm::Int month on data file,
* dd::Int day on data file,
* Format4date::String (optional) date format to identify, default "yyyymmddT"

RETURN:
* file\\_name::String full path of file name corresponding to the given date.

The first case returns:

the second case returns:

"""
function FilePattern(PATH::String, PROD::String, yy::Number, mm::Number, dd::Number; moved=true, Format4Date="yyyymmddT")

    DIRSAR = joinpath(PATH, PROD, "$(yy)")

    !isdir(DIRSAR) && @error "I could not find directory $(DIRSAR)"
    tmp = readdir(DIRSAR, join=true)
    
    today = Date(yy, mm, dd) |> x->Dates.format(x, Format4Date)
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


"""
Read Lead Fraction data from Divergence product by Luisa Von Aldbyll.
```julia-repl
julia> LF = Data_Divergence_LeadFraction(file_name, idx_sector)
julia> LF = Data_Divergence_LeadFraction(file_name, idx_sector, LFₚᵣₒ = :accumulated_2x)
``` 
WHERE:
* file_name::String containing the full path to the sea ice divergence netCDF file,
* idx_sector::Vector{CartesianIndex{2}} Indexes for the data to extract,
* LFₚᵣₒ::Symbol (Optional) which variable to read, the following are supported:
        ** :moved (default) LF drifted corrected to second time stamp,
        ** :none  LF with no drifted corrtected, thus first time stamp,
        ** :accumulated LF accumulated from 10 days
        ** :accumulates_Tx, where T can be 0 to 9, for a specific accumulation.

OUTPUT:
* LF::Dict() With LF data as dictionary with keys :LF, :DIV, or :ACCUMULATES_Xx


"""
function Data_Divergence_LeadFraction(filen::String, idx_sector::Vector{CartesianIndex{2}}; LFₚᵣₒ = :moved)

    # ancilliary functions:
    # 1)
    Missing2NaN(X) = replace(X, missing=>NaN)

    # 2)
    function list_variables(nc::NCDataset, var::Symbol)
        tmp = Dict(Symbol(split(k,"_from")[1])=>k for k in keys(nc) if contains(k, String(var)))
        !isempty(tmp) ? tmp : @error("No variable alike $(var) was found!")
    end
    # ..end of ancilliary functions.

    # Opening given netCDF files:
    nc = NCDataset(filen, "r");
    
    vars = if LFₚᵣₒ == :moved
       Dict(:DIV=>"filtered_div_moved",
            :LF=>"lead_fraction_moved")
    elseif LFₚᵣₒ == :none
        Dict(:DIV=>"filtered_divergence",
             :LF=>"lead_fraction")
    elseif contains(String(LFₚᵣₒ), String(:accumulated))
        list_variables(nc, LFₚᵣₒ)
    end

    #vars = merge(vars, Dict(:x=>"x", :y=>"y"))
    
    sar = Dict(k => Missing2NaN(nc[v][idx_sector]) for (k,v) in vars)

    close(nc)

    return sar
end

end # end module "load"

# end of file
