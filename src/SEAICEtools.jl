# Script containing functions for Sea Ice analysis
# mainly thought to be used on the Arctic

using NCDatasets
using Navigation
using Printf

# *********************************************************
# Functions:

"""
θ, ρ = LonLat_To_CenteredPolar(lon\\_p, lat\\_p, lon, lat)

Function to convert Longitude and Latitude to polar coordinates (azimuth, radius)
centered on a given fixed coordinate point.
Input:
    (lon\\_p, lat\\_p) : centered coordinates ::Float64,
            lon, lat : Arrays with domain coordinats ::Array{Float64,2}.

Return:
    θ, ρ : as Arrays centered in (lon\\_p, lat\\_p)

"""
function LonLat_To_CenteredPolar(lon_p::Float64, lat_p::Float64,
                                 lon::Array{Float64,2}, lat::Array{Float64,2} )
    # Finding the nearest (lat, lon) to NSA:
    lonlat =findall(isapprox.(lat_p, lat; atol=1e-2) .& isapprox.(lon_p, lon; atol=1e-2) );
    #minlat, minlon = argmin(abs.(nsa_lat .- lat)), argmin(abs.(nsa_lon .- lon));
    # Converting all latidue and longitude grid to azimuth and distance
    # with help of Navigation package:
    P_nsa = Point(lat_p, lon_p);
    P_all = Point.(lat, lon);
    ρ_all = map(x -> distance(P_nsa, x), P_all);
    ρ_all /= 1f3;   # [km]
    θ_all = map(x -> bearing(P_nsa, x), P_all);

    return θ_all, ρ_all
end

## ---------------------------------------------------------------------

"""
indexes = Get\\_Sector\\_Indexes(θ₀::Float32, θ₁::Float32; R_lim=50f0)

Given an Array of azimuthal angles and radius the function returns 
the Array's indexes corresponding to a sector or semi-circle W-N-E.
The sector is limited by initial and final azimuth angles.

NOTE: the Arrays named θ_all and D_all must be defined in memory! 

Input:
    θ₀, θ₁ : Initial and final azimtuh :: Float32
    R\\_lim  : Optional radius in [km] :: Float32 (default R\\_lim=50.0 km)
Return:
    indexes: Array's indexes laying within the sector :: Int64

"""
function Get_Sector_Indexes(θ₀::Float32, θ₁::Float32; R_lim=50f0)
    θ₀ , θ₁ = sort([θ₀ , θ₁]);
    idx = θ₀ < 0.0 ?
        ((θ₀ + 360f0) .≤ θ_all .< 360f0 ) .| (θ_all .≤ θ₁) : (θ₀ .≤ θ_all .≤ θ₁);

    idx .&= (ρ_all .≤ R_lim);

    return findall(idx)
end

# --------------------------------------------------------------------------

"""
        indexes = Get\\_Azimuthal\\_Indexes(Wind_Direction, θ; Δθ )

Function to extract the indexes corresponding to a given fixed angle from an
Array of azimuthal angles from a sector of interest.

"""
function Get_Azimuthal_Indexes(θ_wind::Float32,
                               θ_sec::Array{Float64,1},
                               ρ_sec::Array{Float64,1};
                               Δθ = 5f0, R_lim = 50f0)
    # Δθ  is the tolerance for azimuthal angle.
    if @isdefined(θ_sec) && @isdefined(ρ_sec)
        idx_wind = findall(isapprox.(θ_sec, θ_wind, atol=Δθ) .& (ρ_sec .< R_lim));
        return idx_wind
    else
        error("Arrays θ_sec and/or D_sec are not defined as global variables!")
    end
end

# ---------------------------------------------------------------------------

"""
Function to extract SIC from a given azimuthal indexes:
"""
function SIC_from_Azimuth(θ_wind::Float32,
                          θ_sec::Array{Float64,2},
                          D_sec::Array{Float64,2};                          
                          Δθ = 5f0, R_lim = 50f0)
    # Δθ  is the tolerance for azimuthal angle.
    idx_wind = Get_Azimuthal_Indexes(θ_wind, θ_sec)
       
    return lon[idx_wind][:], lat[idx_wind][:], SIC[idx_wind][:]
end

# -----------------------------------------------------------------------------
"""
Function to create the coordinates for drawing a semi-circle.
semi-circle centeres at P_nsa,
Initial and final azimuthal angles: θ₀ , θ₁
Optional parameters are: R_lim (default 50km) and Ncirc (default 180)
    
"""
function Create_Semi_Circle(P_nsa::Point, θ₀::Float64, θ₁::Float64;
                            R_lim=50e0, Ncirc=180)
    # * semi-circle to draw borders:
    R_lim *= 1e3;
    
    θ₀ = ifelse(θ₀ > θ₁, θ₀ -= 360e0, θ₀)
        
    θ = Base.range(θ₀, stop=θ₁, length=Ncirc);
    ρ = map(x-> R_lim, ones(Ncirc)); # [m];
    y = map(x-> P_nsa, ones(Ncirc));
    P_circ = map(destination_point, y, ρ, θ);

    return P_circ
end

# --------------------------------------------------------------------------

function Get_LonLat_From_Point(P₀::Point)
    y = map(x-> x.ϕ, P₀);
    x = map(x-> x.λ, P₀);
    return x, y
end

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
# +++

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
# +++

## ***********************************
## **** temporal helper functions ****
# get file from pattern:
function getFilePattern(path::String, product::String, yy, mm ,dd)
    base_dir = joinpath(path, product, @sprintf("%04d", yy))
    list_file = readdir(base_dir, join=true)
    pattern = @sprintf("%04d%02d%02d", yy, mm, dd)
    ofile = filter(x->all(occursin.(pattern, x)), list_file)
    ofile = isempty(ofile) ? "$dd.$mm.$yy none" : ofile[1]  
    return ofile
end
# +++

# Get the wind vector from Model/RadioSonde
function getSondeData(sonde_file::String; vars=(:WD, :WS, :T) )
    ncin = NCDataset(sonde_file)
    T_C  = copy(ncin["temp"][:,:])
    rh = copy(ncin["rh"][:,:])
    time = ncin["time"][:]
        height= copy(ncin["height"][:]) # km
    wind = copy(ncin["wspd"][:,:])  # m/s
    miss_val = ncin["wspd"].attrib["missing_value"]
    wind[wind .≈ miss_val] .= NaN
    wdir = copy(ncin["wdir"][:,:])  # deg
    miss_val = ncin["wdir"].attrib["missing_value"]
    wdir[wdir .≈ miss_val] .= NaN

    uwind = ncin["u_wind"][:,:]; # m/s
    miss_val = ncin["u_wind"].attrib["missing_value"]
    uwind[uwind .≈ miss_val] .= NaN

    vwind = ncin["v_wind"][:,:]; # m/s
    miss_val = ncin["v_wind"].attrib["missing_value"]
    vwind[vwind .≈ miss_val] .= NaN

    ## Potential temperature
    θ = ncin["potential_temp"][:,:] # K
    miss_val = ncin["potential_temp"].attrib["missing_value"]
    θ[θ .≈ miss_val] .= NaN
    
    ## New variables for IVT
    qv = ncin["sh"][:,:];  # gr/gr
    miss_val = ncin["sh"].attrib["missing_value"]
    qv[qv .≈ miss_val] .= NaN
    
    pa = ncin["bar_pres"][:,:];  # kPa
    miss_val = ncin["bar_pres"].attrib["missing_value"]
    pa[pa .≈ miss_val] .= NaN
    
    close(ncin)
    return Dict(:time=>time, :height=>height[:,1], :T=>T_C) #, U=uwind, V=vwind, RH=rh, WS=wind, WD=wdir, QV=qv, Pa=pa, θ = θ)

end


# Calculating Integrated Water Vapour Transport
function getIVT(rs)
    ΔP = rs.Pa[2:end,:] - rs.Pa[1:end-1,:]
    ΔP *= 1e3 # [Pa]
    tmp_u = rs.QV .* rs.U
    tmp_v = rs.QV .* rs.V
    z_top = 235  # index of the top layer to consider ≈ 306 hPa.
    # calculating IVT components u and v:
    IVT_u = map(1:z_top) do z
        sum(tmp_u[z:z_top,:] .* ΔP[z:z_top,:], dims=1)
    end
    IVT_v = map(1:z_top) do z
        sum(tmp_v[z:z_top,:] .* ΔP[z:z_top,:], dims=1)
    end
    
    IVT_u = vcat(IVT_u...)
    IVT_v = vcat(IVT_v...)
    IVT_vec = cat(IVT_u, IVT_v, dims=3)
    IVT_vec /= 9.81
    # calculating total IVT:
    IVT = sqrt.(IVT_u.^2 + IVT_v.^2)
    IVT /= 9.81
    
    return IVT, IVT_vec
end

# Get Maximum IVT index and values
function getMaxIVT_θ(IVT::Array{Float64,2})
    
    ii = argmax(IVT, dims=1)
    # converting the CartesianIndexes ii from 2-D to 1-D :
    ii = vcat(ii...)

    # retrieving the max of IVT into a vector:
    IVTmax = IVT[ii][:]
    
    return ii, IVTmax
end

# Calculating Richardson Number
function Ri_N(rs)
    θv = rs.θ.*(1.0 .+ 0.6*rs.QV)  # K
    Δθv = θv[2:end,:] .- θv[1,:]'  # K
    ΔZ = 1e3*(rs.height[2:end] .- rs.height[1])  # m
    ΔU = rs.U[2:end,:] .- rs.U[1,:]'  # m/s
    N2 = (Δθv./ΔZ)./θv[1:end-1,:]
    N2 *= 9.81

    tmp = (ΔU./ΔZ).^2
    Ri = N2./tmp
    return N2, Ri
end

# end of script.
