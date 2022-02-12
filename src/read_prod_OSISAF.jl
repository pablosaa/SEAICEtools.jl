# File part of SEAICEtools.jl
# this file include function to read the SIC products provided
# by the Ocean and Sea Ice (OSI) Satellite Aplication Facillity (SAF)
#
# See LICENSE


# --------------------------------------------------------------------------------
# Function to load the latitude and longitude grid data for OSI SAF product OSI-401-b
# provided  by the Norwegian Met Office.
"""
Function read the latitude and longitude coordinates for the gridded Sea Ice data.
Supported products are OSI-401-b product.

The return value is a Array of Points i.a. Array{Points{T}}
"""
function read_LatLon_OSISAF_product(LATLON_FILE::String)

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
Data supported is from MET Office Norway repository.
Two type of products are supported:
* OSISAF 401-b (10 km resolution)
* 

USAGE:
To read OSISAF Sea Ice database use for example:
> SIC = read_SIC_OSISAF_product(filename::String, idx_sector::Vector{CartesianIndex})
where idx_sector is a vector in element type CartesianIndex containing the indexes of the
coordinates from which SIC data will be extracted and read.

> SIC = read_SIC_OSISAF_product(filename::String, idx_sector::Vector{CartesianIndex}, SICPRO="total_uncertainty")
optionally uses the variable SICPRO="ice_conc_unfiltered" to read the unfiltered product for the "ice_conc" default variable.

The return variable is a Vector containing the Sea Ice Concentration data for the specified sector by idx_sector.
In case data is not found or else, SIC is returned as nothing.


"""
function read_SIC_OSISAF_product(sic_filen::String, idx_sector::Vector{CartesianIndex{2}}; SICPROD="ice_conc")

    !isfile(sic_filen) && error("$(sic_filen) can not be found!")
    !in(SICPROD, ("ice_conc", "ice_conc_unfiltered", "total_uncertainty")) && error("$(SICPROD) not supported!")

    ncin = NCDataset(sic_filen, "r");

    SIC = if haskey(ncin, SICPROD)
        tmp = isempty(idx_sector) ? ncin[SICPROD][:] : ncin[SICPROD][idx_sector]
        sic_out = Vector{eltype(skipmissing(tmp))}(undef, size(idx_sector));
        sic_out .= NaN
        ii_nomis = findall(.!ismissing.(tmp))
        sic_out[ii_nomis] = tmp[ii_nomis]
        sic_out
    else
        @warn("$(SICPROD) not found in file. Nothing returned!")
        nothing
    end

    return SIC
end

##using NCDatasets
##using Navigation
##
##LATLON_FILE = "/home/psgarfias/Downloads/ice_conc_nh_polstere-100_multi_201912311200.nc";
##include(joinpath(homedir(), "LIM/repos/SEAICEtools.jl/src/SEAICEtools.jl"))
##
##nsa_coor = Point(71.323e0, -156.609e0)
##R_lim = 50f0
##θ₀ = -133f0
##θ₁ = 113f0
##LonLim = (-158.8, -154.2)
##LatLim = (70.95, 72)
##
##idx_box = SEAICE.extract_LonLat_Box(LonLim, LatLim, sic_coor);
##θ_all, ρ_all = SEAICE.LonLat_To_CenteredPolar(nsa_coor, sic_coor);
##idx_sector = SEAICE.Get_Sector_Indexes(θ₀, θ₁, θ_all, ρ_all, R_lim=50f0);
##
##sic_filen=LATLON_FILE;
##SICPROD = "ice_conc";
##
##lat = map(p->p.ϕ, sic_coor[idx_sector]);
##lon = map(p->p.λ, sic_coor[idx_sector]);

# end of file
