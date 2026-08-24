include("plot_gpr_sections_between_gravimetry.jl")

const PRESCRIBED_DEPTH_SCALE = 1.0
const HIGH_DEPTH_SCALE = 1.35

function prescribed_gap_shape(section, b; depth_scale=PRESCRIBED_DEPTH_SCALE)
    valid_mask = isfinite.(section.bed)
    run_starts = findall(i -> valid_mask[i] && (i == 1 || !valid_mask[i-1]), eachindex(valid_mask))
    run_stops = findall(i -> valid_mask[i] && (i == length(valid_mask) || !valid_mask[i+1]), eachindex(valid_mask))
    length(run_starts) >= 2 || error("P17 does not contain an internal GPR gap")
    left_edge, right_edge = run_stops[1], run_starts[2]

    gap_indices = left_edge:right_edge
    grab_thickness = section.surface[gap_indices] .- section.grab_bed[gap_indices]
    finite = findall(isfinite, grab_thickness)
    isempty(finite) && error("P17 has no finite Grab thickness inside the GPR gap")
    local_bottom = finite[argmax(grab_thickness[finite])]
    bottom_index = first(gap_indices) + local_bottom - 1
    prescribed_thickness = depth_scale * grab_thickness[local_bottom]
    x0 = section.distance[bottom_index]
    z0 = section.surface[bottom_index] - prescribed_thickness

    fitted = anchored_power_values(
        section.distance,
        section.distance[left_edge], section.bed[left_edge],
        section.distance[right_edge], section.bed[right_edge],
        x0, z0, Float64(b),
    )
    (; b=Float64(b), fitted, x0, z0, prescribed_thickness,
       left_edge, right_edge, bottom_index)
end

function plot_p17_map!(layout, section)
    map_paths = [joinpath(ALETSCH, "swiss-map-raster25_2023_1269_krel_1.25_2056.tif"),
                 joinpath(ALETSCH, "swiss-map-raster25_2023_1249_krel_1.25_2056.tif")]
    maps = read_map_rgb.(map_paths)
    ax = Axis(layout, title="P17 location (LV95)", xlabel="Easting (m)",
              ylabel="Northing (m)", aspect=DataAspect())
    for map in maps
        image!(ax, map.extent[1]..map.extent[2], map.extent[3]..map.extent[4], map.rgb)
    end
    lines!(ax, section.x, section.y; color=:deeppink2,
           linewidth=2.8, label="Grab et al. (2021) cross-section")
    lines!(ax, section.map_x, section.map_y; color=:black, linewidth=3.2,
           label="Airborne GPR bed reflections (2011)")
    valid = findall(isfinite, section.map_x)
    index = valid[cld(length(valid), 2)]
    scatter!(ax, [section.x[index]], [section.y[index]]; marker=:rect,
             markersize=Vec2f(54, 40), color=:white, strokecolor=:black,
             strokewidth=2, marker_offset=(0, 38))
    text!(ax, section.x[index], section.y[index]; text="17", color=:black,
          fontsize=25, font=:bold, align=(:center, :center), offset=(0, 38))
    axislegend(ax; position=:rt, labelsize=16, framevisible=false)
    xlims!(ax, minimum(section.x)-500, maximum(section.x)+500)
    ylims!(ax, minimum(section.y)-500, maximum(section.y)+500)
    ax
end

function main()
    surface_raster = read_raster(joinpath(ALETSCH, "alti3D_2m_Aletsch_2023.tif"))
    grab_bed_raster = read_raster(joinpath(ALETSCH, "B36-26_GlacierBed_Grab2021.tif"))
    sections = profiles_between_gravimetry(surface_raster, grab_bed_raster)
    section = only(filter(s -> s.label == "2011 F1-P17", sections))
    fit_b1 = prescribed_gap_shape(section, 1.0)
    fit_b2 = prescribed_gap_shape(section, 2.0)
    deep_b1 = prescribed_gap_shape(section, 1.0; depth_scale=HIGH_DEPTH_SCALE)
    deep_b2 = prescribed_gap_shape(section, 2.0; depth_scale=HIGH_DEPTH_SCALE)

    fig = Figure(size=(1700, 950))
    ax = Axis(fig[1, 1], title="P17: power-law fits anchored at the GPR gap boundaries",
              xlabel="Profile distance (m)", ylabel="Bed elevation (m a.s.l.)",
              aspect=DataAspect())
    lines!(ax, section.distance, section.surface; color=:dodgerblue3,
           linewidth=2.3, label="2023 surface")
    lines!(ax, section.distance, section.grab_bed; color=:deeppink2,
           linewidth=2.0, label="Grab et al. (2021) cross-section")
    lines!(ax, section.distance, section.bed; color=:black,
           linewidth=3.0, label="Measured GPR bed")
    lines!(ax, section.distance, fit_b1.fitted; color=:darkorange,
           linestyle=:dash, linewidth=2.8, label="b = 1, α = 1.0")
    lines!(ax, section.distance, fit_b2.fitted; color=:purple4,
           linestyle=:dash, linewidth=2.8, label="b = 2, α = 1.0")
    lines!(ax, section.distance, deep_b1.fitted; color=:darkorange,
           linestyle=:dot, linewidth=3.2, label="b = 1, α = 1.35")
    lines!(ax, section.distance, deep_b2.fitted; color=:purple4,
           linestyle=:dot, linewidth=3.2, label="b = 2, α = 1.35")
    text!(ax, 0.5, 0.52; space=:relative,
          text="prescribed max. thickness: " *
               "$(round(Int, fit_b1.prescribed_thickness)) m (α = 1.0), " *
               "$(round(Int, deep_b1.prescribed_thickness)) m (α = 1.35)",
          color=:black, fontsize=18, font=:bold, align=(:center, :center))

    scatter!(ax, [Point2f(0.935, 0.09)]; space=:relative, marker=:rect,
             markersize=Vec2f(54, 40), color=:white, strokecolor=:black,
             strokewidth=2)
    text!(ax, 0.935, 0.09; text="17", space=:relative, color=:black,
          fontsize=26, font=:bold, align=(:center, :center))
    axislegend(ax; position=:lb, labelsize=16, framevisible=false)


    output = joinpath(OUTPUT_DIR, "aletsch_p17_gpr_cross_section.png")
    save(output, fig; px_per_unit=1.5)
    println("Saved $output")
end

main()
