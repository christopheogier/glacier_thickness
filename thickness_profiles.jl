include("function.jl")

data_root = joinpath(dirname(@__DIR__), "data")
aletsch_data_dir = joinpath(data_root, "Aletsch")
profiles_dir = joinpath(aletsch_data_dir, "profiles_gravi_aletsch")
gravimetry_path = joinpath(profiles_dir, "gps_points_20260513.txt")
output_dir = joinpath(@__DIR__, "output", "profiles")

dem_path = joinpath(aletsch_data_dir, "B36-26_GlacierBed_Grab2021.tif")
surface_path = joinpath(aletsch_data_dir, "alti3D_2m_Aletsch_2023.tif")
historical_surface_path = joinpath(aletsch_data_dir, "aletsch_DEM_surface_2017xxxx_dsm_lv95_ln02.tif")
uncertainty_minus_path = joinpath(aletsch_data_dir, "B36-26_IceThicknessUncertaintyMinus_Grab2021.tif")
uncertainty_plus_path = joinpath(aletsch_data_dir, "B36-26_IceThicknessUncertaintyPlus_Grab2021.tif")
profile_paths = [
    joinpath(profiles_dir, "Profile_1.kml"),
    joinpath(profiles_dir, "Profile_2.kml"),
]

output_paths = profile_dem_values(
    dem_path,
    profile_paths;
    output_dir,
    source_epsg=4326,
    surface_path,
    historical_surface_path,
    profile_extension=200.0,
    uncertainty_minus_path,
    uncertainty_plus_path,
)

for output_path in output_paths
    println("Saved profile to: ", output_path)
end

profile_plot_path = joinpath(@__DIR__, "output", "plots", "aletsch_bedrock_profiles.png")
plot_profile_elevations(
    output_paths;
    output_path=profile_plot_path,
    labels=["Profile 1", "Profile 2"],
    gravimetry_path,
    gravimetry_max_offset=15.0,
    surface_label="SwissALTI3D surface 2023",
    historical_surface_label="SwissALTI3D surface 2017",
    distance_origin_stations=[111, 116],
    assumed_snow_depths=[1.5, 0.5], # Profile 1 and Profile 2 measured snow depths (m).
)
println("Saved profile plot to: ", profile_plot_path)
