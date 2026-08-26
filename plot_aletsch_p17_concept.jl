include("aletsch_cross_section_plotting.jl")

const DEPTH_SCALE = 1.3

function main()
    surface = read_raster(joinpath(ALETSCH, "alti3D_2m_Aletsch_2023.tif"))
    grab = read_raster(joinpath(ALETSCH, "B36-26_GlacierBed_Grab2021.tif"))
    uncertainty_minus = read_raster(joinpath(ALETSCH, "B36-26_IceThicknessUncertaintyMinus_Grab2021.tif"))
    frank = read_raster(joinpath(ALETSCH, "RGI60-11.01222_bedrock_franck2026_on_grab2021_grid.tif"))
    uncertainty_plus = read_raster(joinpath(ALETSCH, "B36-26_IceThicknessUncertaintyPlus_Grab2021.tif"))
    sections = profiles_between_gravimetry(surface, grab)
    section = only(filter(s -> s.label == "2011 F1-P17", sections))
    comparison_section = only(filter(s -> s.label == "2011 F1-P21", sections))
    output = joinpath(OUTPUT_DIR, "aletsch_p17_concept.png")
    plot_concept_cross_section(section, sections;
                               depth_scale=DEPTH_SCALE,
                               gravity_index=1,
                               comparison_section=comparison_section,
                               frank_raster=frank,
                               uncertainty_minus=uncertainty_minus,
                               uncertainty_plus=uncertainty_plus,
                               output_path=output)
    println("Saved $output")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
