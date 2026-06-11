include("function.jl")

data_root = joinpath(dirname(@__DIR__), "data")
aletsch_data_dir = joinpath(data_root, "Aletsch")
biafo_data_dir = joinpath(data_root, "data_biafo")
plots_dir = joinpath(@__DIR__, "output", "plots")

Aletsch_consensus_id = "11.01450"
Aletsch_Frank_id = "11.01222"
Aletsch_sgi_id = "B36-26"

surface_dem_path = joinpath(aletsch_data_dir, "alti3D_2m_Aletsch_2023.tif")
swissalti3d_csv_path = joinpath(aletsch_data_dir, "ch.swisstopo.swissalti3d-BNvIyZl9.csv")
swissalti3d_tiles_dir = joinpath(aletsch_data_dir, "swissalti3d_tiles")

# Build the SwissALTI3D surface DEM mosaic once. Later runs reuse the saved GeoTIFF.
build_swissalti3d_mosaic(
    swissalti3d_csv_path,
    surface_dem_path;
    tiles_dir=swissalti3d_tiles_dir,
)

consensus_thickness_path = joinpath(aletsch_data_dir, "RGI60-$Aletsch_consensus_id" * "_thickness_Farinotti2019.tif")
glate_only_bedrock_path = joinpath(aletsch_data_dir, "$Aletsch_sgi_id" * "_bed_glateonly_grab2021.asc")
grab2021_bedrock_path = joinpath(aletsch_data_dir, "$Aletsch_sgi_id" * "_GlacierBed_Grab2021.tif")
frank2026_bedrock_path = joinpath(aletsch_data_dir, "RGI60-$Aletsch_Frank_id" * "_bedrock_franck2026.tif")

glate_only_thickness_path = thickness_from_surface_and_bedrock(
    surface_dem_path,
    glate_only_bedrock_path,
    joinpath(aletsch_data_dir, "$Aletsch_sgi_id" * "_thickness_surface_minus_glateonly_grab2021.tif"),
)

grab2021_thickness_path = thickness_from_surface_and_bedrock(
    surface_dem_path,
    grab2021_bedrock_path,
    joinpath(aletsch_data_dir, "$Aletsch_sgi_id" * "_thickness_surface_minus_grab2021.tif"),
)

frank2026_thickness_path = thickness_from_surface_and_bedrock(
    surface_dem_path,
    frank2026_bedrock_path,
    joinpath(aletsch_data_dir, "RGI60-" * Aletsch_Frank_id * "_thickness_surface_minus_frank2026.tif");
    mask_path=joinpath(aletsch_data_dir, "RGI60-" * Aletsch_Frank_id * "_thickness_frank2026.tif"),
    force=true,
)

frank2026_plot_path = joinpath(aletsch_data_dir, "RGI60-" * Aletsch_Frank_id * "_thickness_surface_minus_frank2026_on_grab2021_grid.tif")
_warp_raster_to_match(frank2026_thickness_path, grab2021_thickness_path, frank2026_plot_path; method="bilinear")

milan2022_aletsch_thickness_path = joinpath(aletsch_data_dir, "THICKNESS_RGI-11_2021July09.tif")
milan2022_aletsch_plot_path = joinpath(aletsch_data_dir, "THICKNESS_RGI-11_2021July09_on_grab2021_grid.tif")
_warp_raster_to_match(milan2022_aletsch_thickness_path, grab2021_thickness_path, milan2022_aletsch_plot_path; method="bilinear")

thickness_paths = [
    consensus_thickness_path,
    milan2022_aletsch_plot_path,
    grab2021_thickness_path,
    frank2026_plot_path,
]

labels = [
    "Consensus 2019 (surf. ~2000-2015)",
    "Milan 2022",
    "Grab 2021 (surf. 2023)",
    "Frank 2026 (surf. 2023)",
]

output_path = joinpath(plots_dir, "aletsch_thick_overview.png")
plot_aletsch_thickness_overview(thickness_paths, labels; output_path, outline_path=grab2021_thickness_path)
println("Saved figure to: ", output_path)

# Plot one bedrock difference map: Grab 2021 minus Frank 2026.
frank2026_bedrock_plot_path = joinpath(aletsch_data_dir, "RGI60-" * Aletsch_Frank_id * "_bedrock_franck2026_on_grab2021_grid.tif")
_warp_raster_to_match(frank2026_bedrock_path, grab2021_bedrock_path, frank2026_bedrock_plot_path; method="bilinear")

bedrock_difference_path = joinpath(aletsch_data_dir, "B36-26_bedrock_difference_grab2021_minus_frank2026.tif")
bedrock_difference_raster(grab2021_bedrock_path, frank2026_bedrock_plot_path, bedrock_difference_path; force=true)

bedrock_output_path = joinpath(plots_dir, "aletsch_bedrock_difference.png")
plot_aletsch_bedrock_difference(bedrock_difference_path; output_path=bedrock_output_path, outline_path=grab2021_thickness_path)
println("Saved figure to: ", bedrock_output_path)



# Biafo comparison: use the consensus raster as the plotting extent/grid.
Biafo_consensus_id = "14.00005" # just the Biafo
Biafo_Frank_id = "14.00001" # this is the Siachen complex, which includes Biafo

biafo_consensus_thickness_path = joinpath(biafo_data_dir, "RGI60-" * Biafo_consensus_id * "_thickness_consensus.tif")
biafo_frank_thickness_path = joinpath(biafo_data_dir, "RGI60-" * Biafo_Frank_id * "_thickness_frank2026.tif")
biafo_milan_thickness_path = joinpath(biafo_data_dir, "THICKNESS_RGI-13-15.6_2022September22.tif")
biafo_milan_on_consensus_path = joinpath(biafo_data_dir, "THICKNESS_RGI-13-15.6_2022September22_on_consensus_grid.tif")
biafo_frank_on_consensus_path = joinpath(biafo_data_dir, "RGI60-" * Biafo_Frank_id * "_thickness_frank2026_on_consensus_grid.tif")
_warp_raster_to_match(biafo_frank_thickness_path, biafo_consensus_thickness_path, biafo_frank_on_consensus_path; method="bilinear")
_warp_raster_to_match(biafo_milan_thickness_path, biafo_consensus_thickness_path, biafo_milan_on_consensus_path; method="bilinear")

biafo_output_path = joinpath(plots_dir, "biafo_thickness_comparison.png")
plot_biafo_thickness_comparison_with_milan(
    biafo_consensus_thickness_path,
    biafo_milan_on_consensus_path,
    biafo_frank_on_consensus_path;
    output_path=biafo_output_path,
)
println("Saved figure to: ", biafo_output_path)
