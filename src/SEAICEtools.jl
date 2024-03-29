# Script containing functions for Sea Ice analysis
# mainly thought to be used on the Arctic


module SEAICEtools

using Navigation
using Printf
using NCDatasets
using Statistics

# *********************************************************
# Functions:
"""
Returns the box limist for a given center and distance:

```julia-repl
julia> estimate_box(Pcenter, R_lim)
julia> estimate_box(Pcenter, R_lim; δR = 1e3)
```
WHERE:
* Pcenter::Point indicating center of the box,
* R_lim::Float64 latest size of the box in m,
* δR::Float64 (optional) add to R\\_lim, default 1e3 m,

RETURN:
* Tuple of latitudes and longitudes with the limits of the geographical box.

"""
function estimate_box(Pcenter::Point, R_lim; δR = 1e3)
    N_S_lim = (180e0, 0e0) .|> x->Navigation.destination_point(Pcenter, R_lim+δR, x)
    E_W_lim = (-90e0, 90e0) .|> x->Navigation.destination_point(Pcenter, R_lim+δR, x)
    return  map(x->x.λ, E_W_lim), map(x->x.ϕ, N_S_lim)
end
# ----/

# *********************************************************
"""
Function to convert Longitude and Latitude to polar coordinates (azimuth, radius)
centered on a given fixed coordinate point.

```julia-repl
julia> θ, ρ = LonLat_To_CenteredPolar(lon_p, lat_p, lon, lat)
julia> θ, ρ = LonLat_To_CenteredPolar(coor_p, lon, lat)

```
WHERE:
* coor\\_p::Point indicating the reference coordinate,
* lon\\_p and lat\\_p::Float64 : coordinates of reference coordinate,
* lon::Array longitude of the point to compute distance and azimuth,
* lat::Array latitude of the point to compute distance and azimuth,

RETURN:
* θ::Array of azimuth angles of grid relative to (lon\\_p, lat\\_p),
* ρ::Array of distance [km] from grid to coordinate (lon\\_p, lat\\_p),

"""
function LonLat_To_CenteredPolar(Origin_coor::Point{T}, Grid_coor::Matrix{Point{T}}) where T<:Real
    ρ_all = map(x -> distance(Origin_coor, x), Grid_coor);
    ρ_all /= 1f3;   # [km]
    θ_all = map(x -> bearing(Origin_coor, x), Grid_coor);

    return θ_all, ρ_all
end
function LonLat_To_CenteredPolar(Origin_coor::Point{T}, Grid_coor::Vector{Point{T}}) where T<:Real
    return let Target_coor = Grid_coor[:,:] 
        LonLat_To_CenteredPolar(Origin_coor,Target_coor)
    end
end
function LonLat_To_CenteredPolar(lat_p::T, lon_p::T, 
                                 lat::Matrix{T}, lon::Matrix{T} ) where T<:Real

    # Finding the nearest (lat, lon) to NSA:
    #lonlat =findall(isapprox.(lat_p, lat; atol=1e-2) .& isapprox.(lon_p, lon; atol=1e-2) );
    #minlat, minlon = argmin(abs.(nsa_lat .- lat)), argmin(abs.(nsa_lon .- lon));
    # Converting all latidue and longitude grid to azimuth and distance
    # with help of Navigation package:
    P_nsa = Point(lat_p, lon_p);
    P_all = Point.(lat, lon);

    return LonLat_To_CenteredPolar(P_nsa, P_all)
    ###ρ_all = map(x -> distance(P_nsa, x), P_all);
    ###ρ_all /= 1f3;   # [km]
    ###θ_all = map(x -> bearing(P_nsa, x), P_all);

    ###return θ_all, ρ_all
end
# ----/

# ---------------------------------------------------------------------
"""
Function to extract the indexes of Longitude and Latitude within a box
with limits given:
```julia-repl
julia> idx_box = extract_LonLat_Box(lonlim, latlim, grid_coor)
```
WHERE:
* lonlim::Tuple (lon\\_min, lon\\_max) [degrees] longitude of box limits,
* latlim::Tuple (lat\\_min, lat\\_max) [degrees] latitude of box limits,
* grid\\_coor::Matrix{Point} coordinated of the grid,

RETURN:
* idx_box::Vector{CartesianIndex{2}} indexes of the grid within box.

"""
function extract_LonLat_Box(lonlim::Tuple, latlim::Tuple, sic_coor::Matrix{Point{T}}) where T<:Real
	xlon, ylat = Get_LonLat_From_Point(sic_coor)
	return findall((lonlim[1] .≤ xlon .≤ lonlim[2]) .& (latlim[1] .≤ ylat .≤ latlim[2]))
end
# ----/

## ---------------------------------------------------------------------
"""
Given an Array of azimuthal angles and radius the function returns 
the Array's indexes corresponding to a sector or semi-circle W-N-E.
The sector is limited by initial and final azimuth angles.
```julia-repl
julia> indexes = Get_Sector_Indexes(θ₀, θ₁, θ_all::Matrix)
julia> indexes = Get_Sector_Indexes(θ₀, θ₁, θ_all::Matrix; R_lim=70f0)

```

NOTE: the Arrays named θ_all and D_all must be defined in memory! 

WHERE:
* θ₀::Float32, θ₁::Float32 Initial and final azimuth, respectively
* θ\\_all::Matrix{Real} azimuth angles of all points in grid,
* R\\_lim (optional) radius in [km], default R\\_lim=50.0 km

RETURN:
* indexes::Vector{Int64} indexes of grid poing laying within the sector.

"""
function Get_Sector_Indexes(θ₀::T, θ₁::T, θ_all::Matrix, ρ_all::Matrix; R_lim=50f0) where T<:Real
    θ₀ , θ₁ = sort([θ₀ , θ₁]);
    idx = θ₀ < 0.0 ?
        ((θ₀ + 360f0) .≤ θ_all .< 360f0 ) .| (θ_all .≤ θ₁) : (θ₀ .≤ θ_all .≤ θ₁);

    idx .&= (ρ_all .≤ R_lim);

    return findall(idx)
end

# --------------------------------------------------------------------------

"""
Function to extract the indexes from the grid's azimuth angles corresponding.
The extracted indexes correspond to a given fixed azimuth angle e.g. wind direction.
```julia-repl
julia> indexes = Get_Azimuthal_Indexes(windir, θ, ρ; Δθ )
```
WHERE:
* winddir::Real any fixed azimuth angle,
* θ::Vector{Real} grid azimuth angles,
* ρ::Vector{Real} grid distances,
* Δθ::Real (optionl) conical sector to consider, e.g. winddir ± Δθ, default Δθ = 0.3f0
* R\\_lim:: Real (optional) maximum distance [km] to consider, default 50 km

"""
function Get_Azimuthal_Indexes(θ_wind::T,
                               θ_sec::Vector,
                               ρ_sec::Vector;
                               Δθ = 0.3f0, R_lim = 50f0) where T<:Real
    # Δθ  is the tolerance for azimuthal angle.
    thr_θρ = ((Δθ, 0, 15), (Δθ/2, 15, 25), (Δθ/3, 25, 35),(Δθ/6, 35, R_lim))
    U = sind.(θ_sec)
    V = cosd.(θ_sec)
    δW = @. (U - sind(θ_wind))^2 + (V - cosd(θ_wind))^2 |> sqrt

    idx_radial = []
    foreach(thr_θρ) do th
	ii = findall(th[2] .< ρ_sec .≤ th[3])
	jj=findall(≤(th[1]), δW[ii] )
	push!(idx_radial, ii[jj])
    end
    return idx_radial = vcat(idx_radial...)
end
##function Get_Azimuthal_Indexes(θ_wind::T,
##                               θ_sec::Vector,
##                               ρ_sec::Vector;
##                               Δθ = 5f0, R_lim = 50f0) where T<:Real
##    # Δθ  is the tolerance for azimuthal angle.
##    if @isdefined(θ_sec) && @isdefined(ρ_sec)
##        ii = findall(ρ_sec .≤ R_lim)
##        idx_wind = findall(x -> x<Δθ, abs.(θ_wind .- θ_sec[ii]))
##        return ii[idx_wind]
##    else
##        error("Arrays θ_sec and/or ρ_sec are not defined as global variables!")
##    end
##end

# ---------------------------------------------------------------------------

"""
Function to extract SIC from a given azimuthal sector given the indexes.
```julia-repl
julia> lon, lat, SIC = SIC_from_Azimuth(winddir, θsec, Dsec, lon, lat, SIC_grid)
julia> lon, lat, SIC = SIC_from_Azimuth(winddir, θsec, Dsec, lon, lat, SIC_grid, Δθ=3.0)
julia> lon, lat, SIC = SIC_from_Azimuth(winddir, θsec, Dsec, lon, lat, SIC_grid, R_lim=70.0)
```
WHERE:
* winddir::Real fixed azimuth direction to extract data,
* θsec::Array grid azimuth angles where to extract data from,
* Dsec::Array grid distances [km] where to extract data from,
* Δθ::Real (optional) angle range to extract, e.g. winddir±Δθ, default 3f0,
* R\\_lim::Real (optional) maximum distance to consider [km], default 50.0.

RETURN:
* lon, lat::Vector of longitude and latitude of the extracted grid points,
* SIC::Vector sea ice data extracted from the input fixed azimuth direction.

"""
function SIC_from_Azimuth(θ_wind::T,
                          θ_sec::Matrix,
                          D_sec::Matrix,
                          lon::Matrix,
                          lat::Matrix,
                          SIC::Matrix;
                          Δθ = 3f0, R_lim = 50f0) where T<:Real
    # Δθ  is the tolerance for azimuthal angle.
    idx_wind = Get_Azimuthal_Indexes(θ_wind, θ_sec, D_sec, Δθ=Δθ, R_lim=R_lim)
       
    return lon[idx_wind][:], lat[idx_wind][:], SIC[idx_wind][:]
end

# -----------------------------------------------------------------------------
"""
Function to create the coordinates for a semi-circle centered at Pₒ,
Initial and final azimuthal angles: θ₀ , θ₁
Optional parameters are: R\\_lim (default 50km) and Ncirc (default 180)
```julia-repl
julia> P_circ = Create_Semi_Circle(Pₒ, θ₀, θ₁)
julia> P_circ = Create_Semi_Circle(Pₒ, θ₀, θ₁; R_lim=70e0, Ncirc=360)
```
WHERE:
* Pₒ::Point center coordinate in grid,
* θ₀::Real initial azimuth,
* θ₁::Real final azimuth,
* R\\_lim::Float64 (optional) radius [km] of the semi-circle, default 50e3,
* Ncirc::Int64 (optional) number of point to create the semi-circle, default 180.

RETURN:
* P\\_circ::Vector{Point} coordinated of the semi-circle centered at Pₒ
"""
function Create_Semi_Circle(P_nsa::Point, θ₀::T, θ₁::T;
                            R_lim=50e3, Ncirc=180) where T<:Real
    # * semi-circle to draw borders:
    
    θ₀ = ifelse(θ₀ > θ₁, θ₀ -= 360e0, θ₀)
        
    θ = Base.range(Float64(θ₀), stop=Float64(θ₁), length=Ncirc);
    ρ = map(x-> Float64(R_lim), ones(Ncirc)); # [m]; repeat([R_lim], Ncirc)  #
    y = map(x-> P_nsa, ones(Ncirc)); #repeat([P_nsa], Ncirc) #
    P_circ = map(destination_point, y, ρ, θ);

    return P_circ
end

# --------------------------------------------------------------------------

"""
Function to extract the Lontitude and Latitude from a ::Point type variable.
The input variable can be scalar or matrix:
```julia-repl
julia> lon, lat = Get_LonLat_From_Point(P₀::Point)
julia> lon, lat = Get_LonLat_From_Point(P₀::Matrix{Point})
```
WHERE:
* P₀::Point given Point to extract latitude and longitude from,
* P₀::Matrix{Point} same for Matrix of Point,

RETURN:
* lon::Real longitude of the given point P₀,
* lat::Real latitude of the given point P₀.
"""
function Get_LonLat_From_Point(P₀::Point)
    return P₀.λ, P₀.ϕ
end
function Get_LonLat_From_Point(P₀::Vector)
    Lat = map(x-> x.ϕ, P₀);
    Lon = map(x-> x.λ, P₀);

    return Lon, Lat  #map(x-> (x.λ, x.ϕ), P₀) |> x->vcat(x...) 
end
function Get_LonLat_From_Point(P₀::Matrix)
    Lat = map(x-> x.ϕ, P₀);
    Lon = map(x-> x.λ, P₀);

    return Lon, Lat
end    

# *>*>*>*>*>*> Following function not used!!! *>*>*>*>
function Extract_SIC_From_Sector(P_circ::Point)
    ## Extracting SIC for semi-circle:
    x, y = Get_LonLat_From_Point(P_circ)
    N_ci = 20;
    f(i,j) = argmin(abs.(j .- lat) + abs.(i .- lon));
    x_ci = extrema(x) |> z->range(z[1], stop=z[2], length=N_ci) |> collect ;
    y_ci = extrema(y) |> z->range(z[1], stop=z[2], length=N_ci) |> collect ;
    idx_ci = [f(i,j) for i=x_ci, j=y_ci];
    return idx_ci
end

# +++

# --------------------------------------------------------------------------
"""
Function to convert Array of Array into 2D Matrix
     filling with NaNs otherwise.
     INPUT:
         VAR::Array{T,Array{T}} -> variable to be re-arranged
         ρ::Array{T,Array{T}} -> ranges corresponding to VAR
         ρ_in::Vector{T} -> (optional) bins into re-arrange matrix
         myfunc::Function -> (optional) reduction function to apply within bins
     OUTPUT:
         VAR2D::Array{T,2}(N, R) -> VAR into 2D matrix
         where N length of VAR's outter Array, and R length-1 of ρ_in
"""
function MatOfMat_To_Profile(VAR, ρ; ρ_bin = collect(0:5:50), myfunc::Function=mean)
    Nvar = length(VAR)
    Nrng = length(ρ_bin)
        
    VAR2d = Array{Float32,2}(undef, Nvar, Nrng-1)
    for j ∈ 1:Nvar
        VAR2d[j,:] = map(2:Nrng) do i
            tmp_in = findall(ρ_bin[i-1].< ρ[j] .≤ ρ_bin[i])
            isempty(tmp_in) ? NaN : myfunc(VAR[j][tmp_in])
        end
    end
    return ρ_bin, VAR2d
end
"""
Sort the input VAR and Range into a given binned range vector.
> ρ_bin, VAR_bin, ERR_bin = Sort_To_Profile(VAR, ρ)
> ρ_bin, VAR_bin, ERR_bin = Sort_To_Profile(VAR, ρ, ρ_bin=[0:15:50])
> ρ_bin, VAR_bin, ERR_bin = Sort_To_Profile(VAR, ρ, myfunc=median, myerr=var)

INPUTS:
* VAR::Vector containing the unsorted data
* ρ::Vector ranges for every element of VAR
* ρ_bin::Vector (optional) a vector with the binned values to sort
* myfunc::Function (optional) to reduce from all values within a bin
* myerr::Function (optional) to estimate uncertainty from applying the reduction

RETURN:
ρ_bin, VAR_bin, ERR_bin::Vector size(ρ_bin)

"""
function Sort_To_Profile(VAR, ρ; ρ_bin::Vector = Vector(0:5:50), myfunc::Function=mean, myerr::Function=std)

	applyfunc(Z) = filter(!isnan, Z) |> x-> !isempty(x) ? [myfunc(x) myerr(x)] : [NaN NaN]
    idx = [findall(low .< ρ .≤ hig) for (low, hig) in zip(ρ_bin[1:end-1], ρ_bin[2:end])]
	varout = [isempty(k) ? [NaN NaN] : applyfunc(VAR[k]) for k ∈ idx] |> x->vcat(x...)
	
	return ρ_bin[2:end], varout[:,1], varout[:,2]
end
# ----/

# **************************************************************************
# Function to arrange wind direction selections along the range into 2D Arrays
#
"""
> idx_radial, ρ2d, SIC2d, errSIC2d = wind_idx_radial(wind_dir, wind_range, θ_box, ρ_box, SIC; ρ_bin = Vector(0:5:50))

INPUT:
* wind_dir::Vector,
* wind_range::Vector,
* θ_box::Vector,
* ρ_box::Vector,
* SIC::Vector;
* ρ_bin = Vector(0:5:50) (optional)

RETURN:
* idx_radial,
* ρ2d,
* SIC2d,
* errSIC2d

"""
function wind_idx_radial(wind_dir::Vector, wind_range::Vector, θ_box::Vector, ρ_box::Vector, SIC::Vector; ρ_bin = Vector(0:5:50))
	
    # creating empty output variables:
    N_ρ = length(ρ_bin)-1
    N_t = length(wind_dir)
    idx_radial = Vector{Any}(undef, 0)
    medSIC2d = Matrix{AbstractFloat}(undef, N_ρ , N_t)
    stdSIC2d = Matrix{AbstractFloat}(undef, N_ρ , N_t)

    # starting loop over time index:
    for tidx in eachindex(wind_dir)
	idx_vec = [SEAICE.Get_Azimuthal_Indexes(wd, θ_box, ρ_box, R_lim=wr) for (wd, wr) ∈ zip(wind_dir[tidx], wind_range[tidx])]
	# collecting all indexes for all wind directions corresponding to every time step:
	# idx_vec::Vector{Vector{Float, N}} where N is the number of wind directions ∀ N -> variable number of ranges.
	push!(idx_radial, idx_vec)

	# calculating the gridded SIC for a fixed distance ρ (0:5:50) for every time step:
	_, sic2d, stdsic = let idx = vcat(idx_vec...)
	    Sort_To_Profile(SIC[idx], ρ_box[idx], ρ_bin=ρ_bin)
	end

        # feeding into 2D variables
	medSIC2d[:, tidx] = sic2d
	stdSIC2d[:, tidx] = stdsic
    end
    return idx_radial, ρ_bin, medSIC2d, stdSIC2d
end
# ----/

# Sorting the variables to store:
"""
    Function to retrieve the variables to store as MAT file
    using the wind dependent indexes:
    INPUT: ii::Array{Array{Int32,1},1} Indexes as obtained from 
    the function Get_Azimuthal_Indexes()

"""
function Get_MATVAR_From_Index(ii, lon, lat, SIC)
    # Getting the lengths of data corresponding to the indexes
    # of specified wind direction:
    nlen = length.(ii);
    Ncol = maximum(nlen);
    Nrow = length(nlen);

    var_test = true
    # do the following variables exist?
    # idx: indexes of lon, lat for the sector of interest
    var_test &= @isdefined(idx)
    # ρ_all: distance of all pixels within sector
    var_test &= @isdefined(ρ_all)
    # θ_all: bearing angles of all pixels within sector
    var_test &= @isdefined(θ_all)
    # SIC: Sea Ice Concentration of the sector (lon, lat)
    var_test &= @isdefined(SIC)
    # lon and lat: longitude and latitude of the sector
    var_test &= @isdefined(lon)
    var_test &= @isdefined(lat)

    @assert var_test "Some or all Global variables not defined!"

    # TEMPLATE:
    # tmp = map(i->[ρ_all[idx][i]' fill(NaN,1,Ncol-length(i))], idx_w50) |> x->vcat(x...);
    #
    # Getting the values of radius and azimuth for every extraction:
    # dd = map(i -> [ρ_all[idx][i]' fill(NaN, 1, Ncol-length(i))], ii)
    dd = map(i -> ρ_all[idx][i], ii)
    # Getting the indexes sorted as a function of distance:
    tmp = sortperm.(dd);

    # re-arranging the input indexes:
    
    ρ_w = map(x-> ρ_all[idx][x], ii);
    θ_w = map(x -> θ_all[idx][x], ii);
    SIC_w = map(x -> SIC[x], ii);
    lon_w = map(x -> lon[x], ii);
    lat_w = map(x -> lat[x], ii);
    range, SIC2D_w = MatOfMat_To_Profile(SIC_w, ρ_w, myfunc=(round ∘ minimum))
    
    # converting into 2D matrices:
    
    #θ_w = map(i-> θ_w[i, tmp[i]]', 1:Nrow) |> x->vcat(x...);
    return ρ_w, θ_w, SIC_w, lon_w, lat_w, SIC2D_w, range
end
# ----/

# --------------------------------------------------------------
# Function to make pair of vectors for lon, lat from a given list of
# angles and distance to a center point
"""
Function to make pair of vectors for lon, lat from a given list of
angles and distance to a center point
USAGE:

> lon_lines, lat_lines = create_pair_lines(ref_coor, R_lim, wind_dir)

WHERE:
* ref_coor::Point  is the initial reference point
* R_lim::Float64 is the distance in km
* wind_dir::Vector{Float64} is the vector with several angles in degrees
RETURN:
* lon_lines, lat_lines::Vector{Float64} are pairs of points with lon and lat from
R_lim at angle wind_dir[i] to reference point ref_coor

"""
function create_pair_lines(ref_coor::Point{T}, ρ::T, WD::Vector) where T<:Real
    
    P_lines = map(WD) do θ
        destination_point(ref_coor, ρ, θ)
    end
    Nlines = length(WD)
    lon_lines = fill(ref_coor.λ, 2Nlines)
    lat_lines = fill(ref_coor.ϕ, 2Nlines)
    lon_lines[1:2:end], lat_lines[1:2:end] = Get_LonLat_From_Point(P_lines)

    return lon_lines, lat_lines
end
function create_pair_lines(ref_coor::Point{T}, ρ::Vector, WD::Vector) where T<:Real
    Nlines = length(WD)
    length(ρ) != Nlines && error("R_lim and wind_dir must be same length")

    P_lines = [destination_point(ref_coor, r, θ) for (r, θ) ∈ zip(ρ, WD)]
    lon_lines = fill(ref_coor.λ, 2Nlines)
    lat_lines = fill(ref_coor.ϕ, 2Nlines)
    lon_lines[1:2:end], lat_lines[1:2:end] = Get_LonLat_From_Point(P_lines)

    return lon_lines, lat_lines
end
# ----/

# ****************** EXTRA FILES *************************
include("read_prod_UniBremen.jl")

include("read_prod_OSISAF.jl")

include("read_prod_AWI.jl")

include("aux_math_functions.jl")
# ----/

end # .../ end of module

# end of script.

### ----/
