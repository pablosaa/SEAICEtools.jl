# File part of SEAICEtools.jl
# this file include function to read the SIC products provided
# by the University of Bremen
#
# See LICENSE


# --------------------------------------------------------------------------------
# Function to load the latitude and longitude grid data for Sea Ice products by
# the University of Bremen.
"""
Function read the latitude and longitude coordinates for the gridded Sea Ice data.
Supported products are AMSR2 and MODIS_AMSR2 merged product.

The return value is a Array of Points i.a. Array{Points{T}}
"""
function read_LatLon_Bremen_product(LATLON_FILE::String)
    !isfile(LATLON_FILE) && error("$(LATLON_FILE) cannot be found!")
    # checking for supported variables:
    lat_var = "";
    lon_var = "";
    sic_coor = NCDataset(LATLON_FILE, "r") do ds
        if haskey(ds, "lat") & haskey(ds, "lon")
	    lat_var, lon_var = "lat", "lon"
        elseif haskey(ds, "Latitudes") & haskey(ds, "Longitudes")
            lat_var, lon_var = "Latitudes", "Longitudes"
        else
            warm("$(LATLON_FILE) does not contain valid latitude/longitude variables!")
            return nothing
        end
        Point.(Float64.(ds[lat_var][:,:]),
	       Float64.(ds[lon_var][:,:])
	       )
    end;

    return sic_coor
end
# ----/

# --------------------------------------------------------------------------------
# Function to read netCDF files and extract SIC information for given sector
"""
Function to  read netCDF files and extract SIC information for given sector.
Data supported is from University of Bremem Sea Ice repository.
Two type of products are supported:
* AMSR2 (3.124 km resolution)
* merge MODIS_AMSR2 (1km resolution)

USAGE:
To read AMSR2 Sea Ice database use for example:
> SIC = read_SIC_product(filename::String, idx_sector::Vector{CartesianIndex})
where idx_sector is a vector in element type CartesianIndex containing the indexes of the
coordinates from which SIC data will be extracted and read.

> SIC = read_SIC_product(filename::String, idx_sector::Vector{CartesianIndex}, SICPRO="mersic")
optionally uses the variable SICPRO="mersic" to read the MODIS_AMSR2 merge product for the "mersic"
variable.
Supported variables are "ASI Ice Concentration" for AMSR2, or "mersic", "asic" and "msic" for the
MODIS_AMSR2 merge product.

The return variable is a Vector containing the Sea Ice Concentration data for the specified sector
by idx_sector. In case data is not found or else, SIC is returned as nothing.


"""
function read_SIC_Bremen_product(sic_filen::String, idx_sector::Vector{CartesianIndex{2}}; SICPROD="ASI Ice Concentration")
    !isfile(sic_filen) && error("$(sic_filen) can not be found!")
    !in(SICPROD, ("ASI Ice Concentration", "mersic", "asic", "msic")) && error("$(SICPROD) not supported!")
    ncin = NCDataset(sic_filen, "r")
    if haskey(ncin, SICPROD)
        SIC = let flag = isempty(idx_sector)
            tmp = ifelse(flag, ncin[SICPROD][:,:], ncin[SICPROD][idx_sector])
            eltype(tmp) <: AbstractFloat ? tmp : Float32.(tmp)
        end
    else
        @warn("$(SICPROD) not found in file. Nothing returned!")
	SIC = nothing	
    end
    # Checking for Non data values:
    if in(SICPROD, ("mersic", "asic", "msic"))
        
        fill_value = ncin.attrib["fill_value_int8"]
        land_value = ncin.attrib["land_value_int8"]
        SIC[SIC .≥ fill_value] .= NaN32
    end
    
    return SIC
end;


# end of script
