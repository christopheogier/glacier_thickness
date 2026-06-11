include("function.jl")

data_root = joinpath(dirname(@__DIR__), "data")
biafo_data_dir = joinpath(data_root, "data_biafo")
profiles_dir = joinpath(biafo_data_dir, "Hewitt_1986_ice thickness")
output_dir = joinpath(@__DIR__, "output", "profiles", "biafo")

profile_paths = [
    joinpath(profiles_dir, "1986_hewitt_A.shp"),
    joinpath(profiles_dir, "1986_hewitt_B.shp"),
    joinpath(profiles_dir, "1986_hewitt_C.shp"),
]

consensus_path = joinpath(biafo_data_dir, "RGI60-14.00005_thickness_consensus.tif")
milan_path = joinpath(biafo_data_dir, "THICKNESS_RGI-13-15.6_2022September22_on_consensus_grid.tif")
frank_path = joinpath(biafo_data_dir, "RGI60-14.00001_thickness_frank2026_on_consensus_grid.tif")

consensus_profiles = profile_dem_values(consensus_path, profile_paths; output_dir, source_epsg=4326)
milan_profiles = profile_dem_values(milan_path, profile_paths; output_dir, source_epsg=4326)
frank_profiles = profile_dem_values(frank_path, profile_paths; output_dir, source_epsg=4326)

hewitt_points = [
    ([650, 1400, 2100], [450, 450, 450]),
    ([500, 1300, 1600, 1900], [450, 1400, 1300, 1200]),
    ([50, 500, 1000, 1750, 2450, 2550], [300, 350, 700, 650, 500, 400]),
]

profile_plot_path = joinpath(@__DIR__, "output", "plots", "biafo_thickness_profiles.png")
plot_biafo_thickness_profiles(
    [consensus_profiles, milan_profiles, frank_profiles];
    output_path=profile_plot_path,
    profile_labels=["Profile A", "Profile B", "Profile C"],
    model_labels=["Consensus 2019", "Milan 2022", "Frank 2026"],
    hewitt_points,
    xerr=250,
)

println("Saved Biafo profile plot to: ", profile_plot_path)
