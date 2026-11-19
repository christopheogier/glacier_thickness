include("plot_gpr_sections_between_gravimetry.jl")

function internal_gap_edges(section)
    valid = isfinite.(section.bed)
    starts = findall(i -> valid[i] && (i == 1 || !valid[i-1]), eachindex(valid))
    stops = findall(i -> valid[i] && (i == length(valid) || !valid[i+1]), eachindex(valid))
    length(starts) >= 2 || error("$(section.label) has no internal GPR gap")
    gaps = [(stops[i], starts[i+1]) for i in 1:length(starts)-1]
    gaps[argmax([section.distance[right] - section.distance[left] for (left, right) in gaps])]
end

function fit_grab_gap(section; b=2.0)
    left, right = internal_gap_edges(section)
    target = fill(NaN, length(section.bed))
    target[left:right] .= section.grab_bed[left:right]
    target[left], target[right] = section.bed[left], section.bed[right]
    anchors = fill(NaN, length(section.bed))
    anchors[left], anchors[right] = section.bed[left], section.bed[right]
    levels = copy(section.surface)
    levels[left], levels[right] = section.bed[left], section.bed[right]
    anchored_fit((; section.distance, bed=target, surface=levels, anchors); fixed_b=b)
end

function scale_fit_depth(section, fit, scale)
    bottom_index = argmin(abs.(section.distance .- fit.x0))
    surface_at_bottom = section.surface[bottom_index]
    base_thickness = surface_at_bottom - fit.z0
    deep_z0 = surface_at_bottom - scale * base_thickness
    fitted = anchored_power_values(section.distance, fit.xleft, fit.zleft,
                                    fit.xright, fit.zright, fit.x0, deep_z0, fit.b)
    (; fitted, max_thickness=scale * base_thickness)
end

function prescribed_gap_shape(section, b; depth_scale=1.0)
    left, right = internal_gap_edges(section)
    gap = left:right
    thickness = section.surface[gap] .- section.grab_bed[gap]
    finite = findall(isfinite, thickness)
    isempty(finite) && error("$(section.label) has no finite Grab thickness in its GPR gap")
    local_bottom = finite[argmax(thickness[finite])]
    bottom_index = first(gap) + local_bottom - 1
    prescribed_thickness = depth_scale * thickness[local_bottom]
    x0 = section.distance[bottom_index]
    z0 = section.surface[bottom_index] - prescribed_thickness
    fitted = anchored_power_values(section.distance,
        section.distance[left], section.bed[left],
        section.distance[right], section.bed[right], x0, z0, Float64(b))
    (; b=Float64(b), fitted, x0, z0, prescribed_thickness,
       left_edge=left, right_edge=right, bottom_index)
end



function profile_short_label(section)
    result = match(r"P(\d+)$", section.label)
    result === nothing ? section.label : "P" * result.captures[1]
end

function plot_concept_location_map!(layout, section, all_sections; gravity_index=1, additional_highlight=nothing)
    map_paths = [joinpath(ALETSCH, "swiss-map-raster25_2023_1269_krel_1.25_2056.tif"),
                 joinpath(ALETSCH, "swiss-map-raster25_2023_1249_krel_1.25_2056.tif")]
    maps = read_map_rgb.(map_paths)
    gravimetry = read_csv_xy.(GRAVIMETRY_PROFILE_PATHS)
    gravity = gravimetry[gravity_index]
    ax = Axis(layout, title="Aletschgletscher tongue", xlabel="Easting (m)", ylabel="Northing (m)",
              aspect=DataAspect(), titlesize=56, xlabelsize=88, ylabelsize=88,
              xticklabelsize=30, yticklabelsize=30)
    for map in maps
        image!(ax, map.extent[1]..map.extent[2], map.extent[3]..map.extent[4], map.rgb)
    end
    for candidate in all_sections
        lines!(ax, candidate.map_gap_x, candidate.map_gap_y; color=:deeppink2, linewidth=3)
        lines!(ax, candidate.map_x, candidate.map_y; color=:black, linewidth=4)
        valid = findall(isfinite, candidate.map_x)
        isempty(valid) && continue
        if candidate.label == section.label || (additional_highlight !== nothing && candidate.label == additional_highlight.label)
            label_index = valid[cld(length(valid), 2)]
            profile_number = profile_short_label(candidate)
            text!(ax, candidate.x[label_index], candidate.y[label_index];
                  text=profile_number, color=:black, fontsize=50, font=:bold,
                  align=(:center, :center), offset=(0, 38))
        end
    end
    stations = read_gravimetry_points()
    scatter!(ax, getproperty.(stations, :x), getproperty.(stations, :y);
             color=:purple4, markersize=30)

    allx = vcat([s.x for s in all_sections]..., [g.x for g in gravimetry]...)
    ally = vcat([s.y for s in all_sections]..., [g.y for g in gravimetry]...)
    xlims!(ax, extrema(allx) .+ (-500, 500))
    ylims!(ax, extrema(ally) .+ (-500, 500))
    ax
end

function plot_concept_cross_section(section, all_sections; depth_scale=1.3,
                                    gravity_index=1, comparison_section=nothing, frank_raster=nothing, uncertainty_minus=nothing,
                                    uncertainty_plus=nothing, map_first=false, output_path)
    base_fit = fit_grab_gap(section; b=2.0)
    deep_fit = scale_fit_depth(section, base_fit, depth_scale)
    fig = Figure(size=(4800, 2600), figure_padding=100)
    main_column = map_first && comparison_section !== nothing ? 3 : 1
    comparison_column = 2
    map_column = map_first ? 1 : (comparison_section === nothing ? 2 : 3)
    ax = Axis(fig[1, main_column], title="$(profile_short_label(section)) cross-section",
              xlabel="Distance (m)", ylabel="Elevation (m)",
              titlesize=56, xlabelsize=88, ylabelsize=88,
              xticklabelsize=88, yticklabelsize=88, yticks=1600:200:2800)
    lines!(ax, section.distance, section.surface; color=:dodgerblue3, linewidth=4)
    if uncertainty_minus !== nothing && uncertainty_plus !== nothing
        minus_values = [sample_bilinear(uncertainty_minus, x, y) for (x, y) in zip(section.x, section.y)]
        plus_values = [sample_bilinear(uncertainty_plus, x, y) for (x, y) in zip(section.x, section.y)]
        lower = section.grab_bed .- plus_values
        upper = section.grab_bed .- minus_values
        band!(ax, section.distance, lower, upper; color=(:grey55, 0.38))
    end
    lines!(ax, section.distance, section.grab_bed; color=:deeppink2, linewidth=3.5)
    lines!(ax, section.distance, section.bed; color=:black, linewidth=5)
    lines!(ax, section.distance, deep_fit.fitted; color=:red3,
           linestyle=:dash, linewidth=4)
    if frank_raster !== nothing
        frank_bed = [sample_bilinear(frank_raster, x, y) for (x, y) in zip(section.x, section.y)]
        lines!(ax, section.distance, frank_bed; color=:blue3, linewidth=4)
    end
    xlims!(ax, 0, 1500)
    surface_maximum = maximum(filter(isfinite, section.surface))
    ylims!(ax, 1800, surface_maximum + 40)

    if comparison_section !== nothing
        comparison_ax = Axis(fig[1, comparison_column],
            title=profile_short_label(comparison_section) * " cross-section",
            xlabel="Distance (m)", ylabel="Elevation (m)",
            titlesize=56, xlabelsize=88, ylabelsize=88,
            xticklabelsize=88, yticklabelsize=88, yticks=1600:200:2800)
        lines!(comparison_ax, comparison_section.distance, comparison_section.surface;
               color=:dodgerblue3, linewidth=4)
        if uncertainty_minus !== nothing && uncertainty_plus !== nothing
            comparison_minus = [sample_bilinear(uncertainty_minus, x, y)
                                for (x, y) in zip(comparison_section.x, comparison_section.y)]
            comparison_plus = [sample_bilinear(uncertainty_plus, x, y)
                               for (x, y) in zip(comparison_section.x, comparison_section.y)]
            comparison_lower = comparison_section.grab_bed .- comparison_plus
            comparison_upper = comparison_section.grab_bed .- comparison_minus
            band!(comparison_ax, comparison_section.distance, comparison_lower, comparison_upper;
                  color=(:grey55, 0.38))
        end
        lines!(comparison_ax, comparison_section.distance, comparison_section.grab_bed;
               color=:deeppink2, linewidth=3.5)
        lines!(comparison_ax, comparison_section.distance, comparison_section.bed;
               color=:black, linewidth=5)
        comparison_minimum = minimum(vcat(filter(isfinite, comparison_section.bed),
                                          filter(isfinite, comparison_section.grab_bed)))
        if frank_raster !== nothing
            comparison_frank_bed = [sample_bilinear(frank_raster, x, y)
                                    for (x, y) in zip(comparison_section.x, comparison_section.y)]
            lines!(comparison_ax, comparison_section.distance, comparison_frank_bed;
                   color=:blue3, linewidth=4)
            comparison_minimum = min(comparison_minimum,
                                     minimum(filter(isfinite, comparison_frank_bed)))
        end
        xlims!(comparison_ax, 0, maximum(comparison_section.distance))
        comparison_surface_maximum = maximum(filter(isfinite, comparison_section.surface))
        ylims!(comparison_ax, comparison_minimum - 40, comparison_surface_maximum + 40)
    end

    plot_concept_location_map!(fig[1, map_column], section, all_sections;
                               gravity_index=gravity_index,
                               additional_highlight=comparison_section)
    colgap!(fig.layout, 1, 90)
    colgap!(fig.layout, 2, 90)
    colsize!(fig.layout, 1, Relative(0.333))
    colsize!(fig.layout, 2, Relative(0.333))
    comparison_section !== nothing && colsize!(fig.layout, 3, Relative(0.334))
    gravity_element = [MarkerElement(color=:purple4, marker=:circle, markersize=51),
                       LineElement(color=:black, linewidth=7)]
    gravity_label = ["Gravimetry measurements\n(2026)", "Airborne GPR (2011)"]
    Legend(fig[2, map_column], gravity_element, gravity_label;
           orientation=:vertical, labelsize=72, patchsize=(108, 69),
           framevisible=false, tellwidth=false)
    rowsize!(fig.layout, 2, Fixed(300))
    rowgap!(fig.layout, 1, 25)

    save(output_path, fig; px_per_unit=1.5)
    fig
end

function plot_prescribed_gap_comparison(section; exponents=(1.0, 2.0),
                                        depth_scales=(1.0, 1.35), output_path)
    colors = Dict(1.0 => :darkorange, 2.0 => :purple4)
    styles = Dict(depth_scales[1] => :dash, depth_scales[2] => :dot)
    fits = Dict((b, scale) => prescribed_gap_shape(section, b; depth_scale=scale)
                for b in exponents for scale in depth_scales)
    fig = Figure(size=(1700, 950))
    ax = Axis(fig[1, 1], title="$(profile_short_label(section)): prescribed depth",
              xlabel="Profile distance (m)", ylabel="Bed elevation (m a.s.l.)",
              aspect=DataAspect())
    lines!(ax, section.distance, section.surface; color=:dodgerblue3,
           linewidth=2.3, label="2023 surface")
    lines!(ax, section.distance, section.grab_bed; color=:deeppink2,
           linewidth=2.0, label="Grab 2021")
    lines!(ax, section.distance, section.bed; color=:black,
           linewidth=3.0, label="GPR")
    for b in exponents, scale in depth_scales
        fit = fits[(b, scale)]
        lines!(ax, section.distance, fit.fitted; color=colors[b],
               linestyle=styles[scale], linewidth=scale == depth_scales[1] ? 2.8 : 3.2,
               label="b = $(Int(b)), α = $scale")
    end
    base = fits[(first(exponents), first(depth_scales))]
    deep = fits[(first(exponents), last(depth_scales))]
    text!(ax, 0.5, 0.52; space=:relative,
          text="prescribed max. thickness: $(round(Int, base.prescribed_thickness)) m " *
               "(α = $(first(depth_scales))), $(round(Int, deep.prescribed_thickness)) m " *
               "(α = $(last(depth_scales)))",
          color=:black, fontsize=18, font=:bold, align=(:center, :center))
    text!(ax, 0.935, 0.09; text=replace(profile_short_label(section), "P" => ""),
          space=:relative, color=:black, fontsize=26, font=:bold,
          align=(:center, :center))
    axislegend(ax; position=:lb, labelsize=16, framevisible=false)
    save(output_path, fig; px_per_unit=1.5)
    fig
end
