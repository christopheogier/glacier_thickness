include("plot_aletsch_gpr_profiles.jl")

const GRAVIMETRY_PROFILE_PATHS = [
    joinpath(@__DIR__, "output", "profiles", "Profile_1_B36-26_GlacierBed_Grab2021.csv"),
    joinpath(@__DIR__, "output", "profiles", "Profile_2_B36-26_GlacierBed_Grab2021.csv"),
]
const FULL_FIT_CSV = joinpath(OUTPUT_DIR, "aletsch_flow_cross_section_powerlaw_fits.csv")

function full_cross_section_mean_b()
    rows = readlines(FULL_FIT_CSV)[2:end]
    values = [parse(Float64, split(row, ',')[4]) for row in rows if !isempty(strip(row))]
    isempty(values) && error("No full cross-section fits found in $FULL_FIT_CSV")
    sum(values) / length(values)
end

function full_cross_section_b(profile_label)
    for row in readlines(FULL_FIT_CSV)[2:end]
        fields = split(row, ',')
        length(fields) >= 4 && fields[2] == profile_label && return parse(Float64, fields[4])
    end
    error("No full cross-section fit found for $profile_label")
end

function read_csv_xy(path)
    x, y = Float64[], Float64[]
    for row in readlines(path)[2:end]
        fields = split(row, ',')
        length(fields) >= 2 || continue
        px, py = tryparse(Float64, fields[1]), tryparse(Float64, fields[2])
        (px === nothing || py === nothing) && continue
        push!(x, px); push!(y, py)
    end
    (; x, y)
end

function read_gravimetry_points()
    path = joinpath(ALETSCH, "profiles_gravi_aletsch", "gps_points_20260513.txt")
    rows = readlines(path)
    headers = strip.(split(rows[1], ","))
    station_col = findfirst(==("Gravi_n"), headers)
    x_col = findfirst(==("Rechtswert"), headers)
    y_col = findfirst(==("Hochwert"), headers)
    elevation_col = findfirst(==("Orthom. Höhe"), headers)
    [(station=parse(Int, strip(fields[station_col])),
      x=parse(Float64, strip(fields[x_col])), y=parse(Float64, strip(fields[y_col])),
      elevation=parse(Float64, strip(fields[elevation_col])))
     for fields in (split(row, ",") for row in rows[2:end])]
end

function project_station(px, py, x, y, distance)
    best_offset, best_distance = Inf, NaN
    for i in 1:length(x)-1
        dx, dy = x[i+1]-x[i], y[i+1]-y[i]
        length_squared = dx^2 + dy^2
        t = length_squared == 0 ? 0.0 : clamp(((px-x[i])*dx + (py-y[i])*dy)/length_squared, 0.0, 1.0)
        projected_x, projected_y = x[i] + t*dx, y[i] + t*dy
        offset = hypot(px-projected_x, py-projected_y)
        if offset < best_offset
            best_offset = offset
            best_distance = distance[i] + t*sqrt(length_squared)
        end
    end
    (; distance=best_distance, offset=best_offset)
end

function read_gravimetry_section(path, index)
    rows = readlines(path)[2:end]
    x, y, bed, surface = Float64[], Float64[], Float64[], Float64[]
    for row in rows
        fields = split(row, ',')
        length(fields) >= 3 || continue
        px, py = tryparse(Float64, fields[1]), tryparse(Float64, fields[2])
        pz = tryparse(Float64, fields[3])
        ps = length(fields) >= 6 ? tryparse(Float64, fields[6]) : nothing
        (px === nothing || py === nothing) && continue
        push!(x, px); push!(y, py); push!(bed, pz === nothing ? NaN : pz)
        push!(surface, ps === nothing ? NaN : ps)
    end
    distance = zeros(length(x))
    for i in 2:length(x)
        distance[i] = distance[i-1] + hypot(x[i]-x[i-1], y[i]-y[i-1])
    end
    station_distance, station_elevation = Float64[], Float64[]
    for station in read_gravimetry_points()
        match = project_station(station.x, station.y, x, y, distance)
        if match.offset <= 15.0
            push!(station_distance, match.distance)
            push!(station_elevation, station.elevation)
        end
    end
    (; label="Gravimetry Profile $index", x, y, distance, bed, grab_bed=bed, surface,
       bed_color=:deeppink2, dot_color=:purple4, station_distance, station_elevation, northing=sum(y)/length(y))
end

function all_gpr_profile_paths()
    paths = String[]
    for (root, _, files) in walkdir(GPR_DIR), file in files
        endswith(file, "_Bedrock.txt") && push!(paths, joinpath(root, file))
    end
    sort(paths)
end

function profile_name(path)
    match_result = match(r"(\d{4})\d{4}_Aletsch_BGR_flight(\d+)_profil-(\d+)_Bedrock", basename(path))
    match_result === nothing && return splitext(basename(path))[1]
    year, flight, number = match_result.captures
    "$year F$(parse(Int, flight))-P$(parse(Int, number))"
end

function gravimetry_corridor(gravimetry)
    endpoints = map(gravimetry) do profile
        west, east = argmin(profile.x), argmax(profile.x)
        ((profile.x[west], profile.y[west]), (profile.x[east], profile.y[east]))
    end
    [endpoints[1][1], endpoints[1][2], endpoints[2][2], endpoints[2][1]]
end

function point_in_polygon(x, y, polygon)
    inside = false
    j = length(polygon)
    for i in eachindex(polygon)
        xi, yi = polygon[i]; xj, yj = polygon[j]
        if (yi > y) != (yj > y) && x < (xj-xi) * (y-yi) / (yj-yi) + xi
            inside = !inside
        end
        j = i
    end
    inside
end

function profiles_between_gravimetry(surface_raster, grab_bed_raster)
    gravimetry = read_csv_xy.(GRAVIMETRY_PROFILE_PATHS)
    corridor = gravimetry_corridor(gravimetry)
    selected = []
    for path in all_gpr_profile_paths()
        startswith(basename(path), "2011") || continue
        occursin(r"profil-(012|014|016|018)_Bedrock", basename(path)) && continue
        profile = read_profile(path)
        valid = isfinite.(profile.bed) .& (profile.bed .> 0)
        count(valid) >= 5 || continue
        mean_x = sum(profile.x[valid]) / count(valid)
        mean_y = sum(profile.y[valid]) / count(valid)
        forced_boundary_profile = profile_name(path) in ("2011 F1-P11", "2011 F1-P21")
        (point_in_polygon(mean_x, mean_y, corridor) || forced_boundary_profile) || continue
        px, py, distance = copy(profile.x), copy(profile.y), copy(profile.distance)
        measured = [valid[i] ? profile.bed[i] : NaN for i in eachindex(valid)]
        surface = [sample_bilinear(surface_raster, x, y) for (x, y) in zip(px, py)]
        grab_bed = [sample_bilinear(grab_bed_raster, x, y) for (x, y) in zip(px, py)]
        gap_bed = [!valid[i] && isfinite(grab_bed[i]) ? grab_bed[i] : NaN for i in eachindex(valid)]
        surface_extension = fill(NaN, length(measured))
        anchors = copy(measured)
        if profile_name(path) == "2011 F1-P15"
            valid_indices = findall(valid)
            first_pick = first(valid_indices)
            slope_indices = valid_indices[1:min(60, length(valid_indices))]
            fit_d, fit_z = distance[slope_indices], measured[slope_indices]
            mean_d, mean_z = sum(fit_d)/length(fit_d), sum(fit_z)/length(fit_z)
            slope = sum((fit_d .- mean_d) .* (fit_z .- mean_z)) / sum(abs2, fit_d .- mean_d)
            segment_length = hypot(px[2]-px[1], py[2]-py[1])
            ux, uy = (px[2]-px[1])/segment_length, (py[2]-py[1])/segment_length
            anchor_distance = -300.0
            breakpoint_distance = -60.0
            extension_d = collect(anchor_distance:1.0:-1.0)
            extension_x = px[1] .+ extension_d .* ux
            extension_y = py[1] .+ extension_d .* uy
            extension_surface = [sample_bilinear(surface_raster, x, y) for (x,y) in zip(extension_x,extension_y)]
            n_extension = length(extension_d)
            px = vcat(extension_x, px); py = vcat(extension_y, py)
            distance = vcat(extension_d, distance)
            surface = vcat(extension_surface, surface)
            grab_bed = vcat(fill(NaN,n_extension), grab_bed)
            gap_bed = vcat(fill(NaN,n_extension), gap_bed)
            measured = vcat(fill(NaN,n_extension), measured)
            surface_extension = fill(NaN, length(measured))
            anchors = copy(measured)
            new_first_pick = n_extension + first_pick
            anchor_elevation = surface[1]
            first_pick_distance = distance[new_first_pick]
            first_pick_elevation = measured[new_first_pick]
            breakpoint_elevation = first_pick_elevation + slope * (breakpoint_distance-first_pick_distance)
            for i in 1:new_first_pick
                d = distance[i]
                surface_extension[i] = d <= breakpoint_distance ?
                    anchor_elevation + (breakpoint_elevation-anchor_elevation) *
                    (d-anchor_distance)/(breakpoint_distance-anchor_distance) :
                    first_pick_elevation + slope * (d-first_pick_distance)
            end
            anchors[1:new_first_pick] .= surface_extension[1:new_first_pick]
        end
        valid_combined = isfinite.(measured)
        map_x = [valid_combined[i] ? px[i] : NaN for i in eachindex(valid_combined)]
        map_y = [valid_combined[i] ? py[i] : NaN for i in eachindex(valid_combined)]
        map_gap_x = [isfinite(gap_bed[i]) ? px[i] : NaN for i in eachindex(gap_bed)]
        map_gap_y = [isfinite(gap_bed[i]) ? py[i] : NaN for i in eachindex(gap_bed)]
        push!(selected, (; path, label=profile_name(path), x=px, y=py,
                         map_x, map_y, map_gap_x, map_gap_y, distance, bed=measured,
                         grab_bed, gap_bed, surface_extension, anchors, surface,
                         bed_color=:black, northing=mean_y))
    end
    sort!(selected; by=profile -> -profile.northing)
end

function nelder_mead(objective, start, step; iterations=1200)
    n = length(start); simplex = [copy(start) for _ in 1:n+1]
    for i in 1:n; simplex[i+1][i] += step[i]; end
    values = objective.(simplex)
    for _ in 1:iterations
        order = sortperm(values); simplex, values = simplex[order], values[order]
        centroid = sum(simplex[1:n]) / n
        reflected = centroid + (centroid-simplex[end]); fr = objective(reflected)
        if fr < values[1]
            expanded = centroid + 2(reflected-centroid); fe = objective(expanded)
            simplex[end], values[end] = fe < fr ? (expanded,fe) : (reflected,fr)
        elseif fr < values[n]
            simplex[end], values[end] = reflected, fr
        else
            contracted = centroid + 0.5(simplex[end]-centroid); fc = objective(contracted)
            if fc < values[end]
                simplex[end], values[end] = contracted, fc
            else
                for i in 2:n+1
                    simplex[i] = simplex[1] + 0.5(simplex[i]-simplex[1])
                    values[i] = objective(simplex[i])
                end
            end
        end
    end
    simplex[argmin(values)]
end

function anchored_power_values(distance, xleft, zleft, xright, zright, x0, z0, b)
    [x < xleft || x > xright ? NaN :
     x <= x0 ? z0 + (zleft-z0) * ((x0-x)/(x0-xleft))^b :
               z0 + (zright-z0) * ((x-x0)/(xright-x0))^b for x in distance]
end

function anchored_fit(section; fixed_b=nothing)
    anchor_bed = hasproperty(section, :anchors) ? section.anchors : section.bed
    anchors = findall(isfinite, anchor_bed)
    good = findall(isfinite, section.bed)
    length(anchors) >= 2 && length(good) >= 3 || error("Insufficient bedrock constraints")
    left, right = first(anchors), last(anchors)
    xleft, xright = section.distance[left], section.distance[right]
    zleft, zright = section.surface[left], section.surface[right]
    isfinite(zleft) && isfinite(zright) || error("Surface anchors must be finite")
    inside = [i for i in good if xleft <= section.distance[i] <= xright]
    x, z = section.distance[inside], section.bed[inside]
    width = xright-xleft
    imin = argmin(z)
    start_x0 = clamp(x[imin], xleft+0.05width, xright-0.05width)
    start_z0 = minimum(z)
    function loss(x0, z0, b)
        (xleft+0.01width < x0 < xright-0.01width && 0.35 <= b <= 6 && z0 < min(zleft,zright)) || return 1e30
        prediction = anchored_power_values(x, xleft, zleft, xright, zright, x0, z0, b)
        sum(abs2, z .- prediction) / length(z)
    end
    if fixed_b === nothing
        objective_free(p) = loss(p[1], p[2], exp(p[3]))
        p = nelder_mead(objective_free, [start_x0,start_z0,log(1.5)], [0.05width,20.0,0.15])
        x0,z0,logb = p; b=exp(logb)
        rmse = sqrt(objective_free(p))
    else
        b = Float64(fixed_b)
        objective_fixed(p) = loss(p[1], p[2], b)
        p = nelder_mead(objective_fixed, [start_x0,start_z0], [0.05width,20.0])
        x0,z0 = p
        rmse = sqrt(objective_fixed(p))
    end
    fitted = anchored_power_values(section.distance, xleft, zleft, xright, zright, x0, z0, b)
    (; a=NaN,b,x0,z0,xleft,zleft,xright,zright,rmse,fitted)
end

fit_power_shape(section, b) = anchored_fit(section; fixed_b=b)
fit_free_power_shape(section) = anchored_fit(section)

function plot_selection_map!(layout, sections)
    map_paths = [joinpath(ALETSCH, "swiss-map-raster25_2023_1269_krel_1.25_2056.tif"),
                 joinpath(ALETSCH, "swiss-map-raster25_2023_1249_krel_1.25_2056.tif")]
    maps = read_map_rgb.(map_paths)
    gravimetry = read_csv_xy.(GRAVIMETRY_PROFILE_PATHS)
    ax = Axis(layout, title="Selected profiles between gravimetry transects (LV95)",
              xlabel="Easting (m)", ylabel="Northing (m)", aspect=DataAspect())
    for map in maps
        image!(ax, map.extent[1]..map.extent[2], map.extent[3]..map.extent[4], map.rgb)
    end
    for (i, section) in enumerate(sections)
        lines!(ax, section.map_gap_x, section.map_gap_y; color=:deeppink2, linewidth=2.8,
               label=i == 1 ? "Interpolated bedrock (Grab et al., 2021)" : nothing)
        lines!(ax, section.map_x, section.map_y; color=:black, linewidth=3.2,
               label=i == 1 ? "Airborne GPR bed reflections (2011)" : nothing)
        valid = findall(isfinite, section.map_x)
        index = valid[cld(length(valid), 2)]
        text!(ax, section.x[index], section.y[index]; text=section.label,
              color=:black, fontsize=9, font=:bold, align=(:center,:bottom), offset=(0,5))
        profile_match = match(r"P(\d+)$", section.label)
        if profile_match !== nothing
            profile_number = profile_match.captures[1]
            scatter!(ax, [section.x[index]], [section.y[index]]; marker=:rect,
                     markersize=Vec2f(48, 34), color=:white, strokecolor=:black,
                     strokewidth=2, marker_offset=(0, 40))
            text!(ax, section.x[index], section.y[index]; text=profile_number,
                  color=:black, fontsize=22, font=:bold, align=(:center,:center), offset=(0,40))
        end
    end
    for (i, profile) in enumerate(gravimetry)
        lines!(ax, profile.x, profile.y; color=:purple4, linewidth=4,
               label=i == 1 ? "Gravimetry profiles" : nothing)
        index = cld(length(profile.x), 2)
        scatter!(ax, [profile.x[index]], [profile.y[index]]; marker=:rect,
                 markersize=Vec2f(48, 34), color=:white, strokecolor=:purple4,
                 strokewidth=2, marker_offset=(0, 28))
        text!(ax, profile.x[index], profile.y[index]; text="G$i",
              color=:purple4, fontsize=22, font=:bold, align=(:center,:center), offset=(0,28))
    end
    stations = read_gravimetry_points()
    scatter!(ax, getproperty.(stations, :x), getproperty.(stations, :y);
             color=:purple4, markersize= 30, strokecolor=:white, strokewidth=0.8)
    axislegend(ax; position=:rt, labelsize=16)
    allx = vcat([section.x for section in sections]..., [profile.x for profile in gravimetry]...)
    ally = vcat([section.y for section in sections]..., [profile.y for profile in gravimetry]...)
    xlims!(ax, extrema(allx) .+ (-500, 500))
    ylims!(ax, extrema(ally) .+ (-500, 500))
    ax
end

function main()
    surface_raster = read_raster(joinpath(ALETSCH, "alti3D_2m_Aletsch_2023.tif"))
    grab_bed_raster = read_raster(joinpath(ALETSCH, "B36-26_GlacierBed_Grab2021.tif"))
    gpr_sections = profiles_between_gravimetry(surface_raster, grab_bed_raster)
    isempty(gpr_sections) && error("No GPR bed profiles found between the gravimetry profiles")
    gravimetry_sections = [read_gravimetry_section(path, i) for (i, path) in enumerate(GRAVIMETRY_PROFILE_PATHS)]
    sections = vcat(gpr_sections, gravimetry_sections)
    sort!(sections; by=section -> -section.northing)
    v_fits = [fit_power_shape(section, 1.0) for section in sections]
    grab_fits = [fit_free_power_shape((; section.distance, bed=section.grab_bed, surface=section.surface, anchors=hasproperty(section, :anchors) ? section.anchors : section.bed)) for section in sections]
    ncols = 2; nrows = cld(length(sections), ncols)
    fig = Figure(size=(2300, 330*nrows))
    for (i, (section, v_fit, grab_fit)) in enumerate(zip(sections, v_fits, grab_fits))
        row, col = cld(i, ncols), mod1(i, ncols)
        ax = Axis(fig[row, col], title=section.label,
                  xlabel="Profile distance (m)", ylabel="Bed elevation (m a.s.l.)")
        profile_match = match(r"P(\d+)$", section.label)
        if profile_match !== nothing
            scatter!(ax, [Point2f(0.065, 0.09)]; space=:relative, marker=:rect,
                     markersize=Vec2f(54, 40), color=:white, strokecolor=:black, strokewidth=2)
            text!(ax, 0.065, 0.09; text=profile_match.captures[1], space=:relative,
                  color=:black, fontsize=26, font=:bold, align=(:center,:center))
        else
            gravity_match = match(r"Gravimetry Profile (\d+)$", section.label)
            if gravity_match !== nothing
                scatter!(ax, [Point2f(0.065, 0.09)]; space=:relative, marker=:rect,
                         markersize=Vec2f(54, 40), color=:white, strokecolor=:purple4, strokewidth=2)
                text!(ax, 0.065, 0.09; text="G" * gravity_match.captures[1], space=:relative,
                      color=:purple4, fontsize=26, font=:bold, align=(:center,:center))
            end
        end
        lines!(ax, section.distance, section.surface; color=:dodgerblue3, linewidth=2,
               label="2023 surface")
        if hasproperty(section, :surface_extension) && any(isfinite, section.surface_extension)
            lines!(ax, section.distance, section.surface_extension; color=:red3, linewidth=2.5,
                   label="Interpolated to surface")
        end
        if hasproperty(section, :gap_bed)
            lines!(ax, section.distance, section.gap_bed; color=:deeppink2, linewidth=1.2,
                   label="Grab et al. (2021) gap fill")
        end
        bed_label = section.bed_color == :black ? "Measured GPR bed" :
                    "Interpolated profile (Grab et al. 2021) at gravimetry location"
        lines!(ax, section.distance, section.bed; color=section.bed_color,
               linewidth=hasproperty(section, :dot_color) ? 1.2 : 2.5,
               label=bed_label)
        if hasproperty(section, :station_distance)
            scatter!(ax, section.station_distance, section.station_elevation;
                     color=:purple4, markersize=11, strokecolor=:white, strokewidth=0.8)
        end
        lines!(ax, section.distance, v_fit.fitted; color=:darkorange, linewidth=2.5,
               linestyle=:dash, label="V-shaped fit (b = 1)")
        lines!(ax, section.distance, grab_fit.fitted; color=:deeppink2, linewidth=2.5,
               linestyle=:dash, label="Grab et al. (2021) bed fit")
        text!(ax, 0.5, 0.56; space=:relative,
              text="b = 1 (RMSE = $(round(v_fit.rmse,digits=2)) m)",
              color=:darkorange, fontsize=15, font=:bold, align=(:center,:bottom))
        text!(ax, 0.5, 0.49; space=:relative,
              text="Grab fit: b = $(round(grab_fit.b,digits=2)) (RMSE = $(round(grab_fit.rmse,digits=2)) m)",
              color=:deeppink2, fontsize=15, font=:bold, align=(:center,:bottom))
    end
    plot_selection_map!(fig[1:nrows, ncols+1], gpr_sections)
    colsize!(fig.layout, ncols+1, Fixed(850))
    Label(fig[0, 1:ncols+1],
          "Aletsch 2011 airborne GPR partial cross-sections between the two gravimetry profiles",
          fontsize=24)
    output = joinpath(OUTPUT_DIR, "aletsch_gpr_partial_cross-sections.png")
    save(output, fig; px_per_unit=1.5)
    println("Saved $output ($(length(sections)) profiles; GPR picks unaltered, Grab gap fill shown)")
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
