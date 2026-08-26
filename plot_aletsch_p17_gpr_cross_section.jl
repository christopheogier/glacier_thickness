include("aletsch_cross_section_plotting.jl")

const EXPONENTS = (1.0, 2.0)
const PRESCRIBED_DEPTH_SCALES = (1.0, 1.35)

function main()
    surface = read_raster(joinpath(ALETSCH, "alti3D_2m_Aletsch_2023.tif"))
    grab = read_raster(joinpath(ALETSCH, "Grab2019_10m_Aletsch_final.tif"))
    sections = profiles_between_gravimetry(surface, grab)
    section = only(filter(s -> s.label == "2011 F1-P17", sections))
    output = joinpath(OUTPUT_DIR, "aletsch_p17_gpr_cross_section.png")
    plot_prescribed_gap_comparison(section;
                                   exponents=EXPONENTS,
                                   depth_scales=PRESCRIBED_DEPTH_SCALES,
                                   output_path=output)
    println("Saved $output")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
