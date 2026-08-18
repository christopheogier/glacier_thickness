using CairoMakie
import ArchGDAL

const ROOT = dirname(@__DIR__)
const ALETSCH = joinpath(ROOT, "data", "Aletsch")
const GPR_DIR = joinpath(ALETSCH, "GPR2009_2011")
const OUTPUT_DIR = joinpath(@__DIR__, "output", "plots")

# Profiles listed in GPR_cross_sections_2009_2011.txt.
const PROFILE_SPECS = [
    (year=2009, date="20090331", flight=3, number=42),
    (year=2009, date="20090331", flight=3, number=30),
    (year=2009, date="20090331", flight=3, number=28),
    (year=2009, date="20090331", flight=3, number=24),
    (year=2011, date="20110322", flight=1, number=21),
    (year=2011, date="20110322", flight=1, number=24),
    (year=2011, date="20110322", flight=1, number=25),
    (year=2011, date="20110322", flight=1, number=26),
    (year=2011, date="20110322", flight=1, number=27),
    (year=2011, date="20110322", flight=1, number=28),
    (year=2011, date="20110322", flight=1, number=29),
    (year=2011, date="20110322", flight=1, number=30),
    (year=2011, date="20110322", flight=1, number=31),
]

profile_label(s) = "$(s.year) F$(s.flight)-P$(s.number)"

function profile_path(s)
    filename = "$(s.date)_Aletsch_BGR_flight$(s.flight)_profil-$(lpad(s.number, 3, '0'))_Bedrock.txt"
    return joinpath(GPR_DIR, "L2_$(s.year)", filename)
end

"""Read a BGR bedrock text file and convert its LV03 coordinates to LV95."""
function read_profile(path)
    raw = readlines(path)
    rows = [split(strip(line)) for line in raw if !isempty(strip(line))]
    all(length(row) >= 4 for row in rows) || error("Malformed row in $path")
    x03 = parse.(Float64, getindex.(rows, 2))
    y03 = parse.(Float64, getindex.(rows, 3))
    bed = parse.(Float64, getindex.(rows, 4))

    x95 = similar(x03)
    y95 = similar(y03)
    ArchGDAL.importEPSG(21781; order=:trad) do source
        ArchGDAL.importEPSG(2056; order=:trad) do target
            ArchGDAL.createcoordtrans(source, target) do transform
                for i in eachindex(x03)
                    ArchGDAL.createpoint(x03[i], y03[i]) do point
                        ArchGDAL.transform!(point, transform)
                        x95[i], y95[i], _ = ArchGDAL.getpoint(point, 0)
                    end
                end
            end
        end
    end

    distance = zeros(length(x95))
    for i in 2:length(distance)
        distance[i] = distance[i - 1] + hypot(x95[i] - x95[i - 1], y95[i] - y95[i - 1])
    end
    return (; x=x95, y=y95, bed, distance)
end

function read_raster(path)
    ArchGDAL.read(path) do dataset
        band = ArchGDAL.getband(dataset, 1)
        return (array=Float64.(ArchGDAL.read(band)),
                nodata=ArchGDAL.getnodatavalue(band),
                gt=ArchGDAL.getgeotransform(dataset))
    end
end

function sample_bilinear(raster, x, y)
    gt = raster.gt
    det = gt[2] * gt[6] - gt[3] * gt[5]
    px = (gt[6] * (x - gt[1]) - gt[3] * (y - gt[4])) / det + 0.5
    py = (-gt[5] * (x - gt[1]) + gt[2] * (y - gt[4])) / det + 0.5
    nx, ny = size(raster.array)
    (1 <= px <= nx && 1 <= py <= ny) || return NaN
    x0, y0 = clamp(floor(Int, px), 1, nx), clamp(floor(Int, py), 1, ny)
    x1, y1 = min(x0 + 1, nx), min(y0 + 1, ny)
    values = (raster.array[x0, y0], raster.array[x1, y0],
              raster.array[x0, y1], raster.array[x1, y1])
    invalid(v) = !isfinite(v) || (raster.nodata !== nothing && isapprox(v, raster.nodata))
    any(invalid, values) && return NaN
    tx, ty = px - x0, py - y0
    return (1 - ty) * ((1 - tx) * values[1] + tx * values[2]) +
           ty * ((1 - tx) * values[3] + tx * values[4])
end

"""Return bed elevations split into measured picks and internal interpolated gaps."""
function split_bedrock(profile)
    valid = isfinite.(profile.bed) .& (profile.bed .> 0)
    measured = fill(NaN, length(valid))
    interpolated = fill(NaN, length(valid))
    measured[valid] .= profile.bed[valid]
    indices = findall(valid)
    isempty(indices) && return measured, interpolated
    for (left, right) in zip(indices[1:end-1], indices[2:end])
        if right > left + 1
            interpolated[left:right] .= range(profile.bed[left], profile.bed[right], length=right-left+1)
        end
    end
    return measured, interpolated
end

# Project every track in a panel onto one best-fit axis.  The axis is forced
# west-to-east and all profiles share the westernmost projected origin.
function common_profile_distances(group_profiles)
    x = vcat([profile.x for profile in group_profiles]...)
    y = vcat([profile.y for profile in group_profiles]...)
    xc, yc = x .- sum(x) / length(x), y .- sum(y) / length(y)
    cxx, cyy, cxy = sum(abs2, xc), sum(abs2, yc), sum(xc .* yc)
    angle = 0.5 * atan(2cxy, cxx - cyy)
    vx, vy = cos(angle), sin(angle)
    vx < 0 && ((vx, vy) = (-vx, -vy))
    projected = [[vx * px + vy * py for (px, py) in zip(profile.x, profile.y)]
                 for profile in group_profiles]
    western_origin = minimum(vcat(projected...))
    return [distance .- western_origin for distance in projected]
end

function plot_cross_sections(specs, profiles, surface, gravity)
    # Four near-coincident 2009/2011 survey pairs share a subplot.
    groups = [[4, 5], [6], [7], [8], [3, 9], [10], [2, 11], [12], [1, 13]]
    ncols, nrows = 2, cld(length(groups), 2)
    fig = Figure(size=(2150, 330 * nrows))
    for (i, group) in enumerate(groups)
        title = join(profile_label.(specs[group]), " / ")
        ax = Axis(fig[cld(i, ncols), mod1(i, ncols)], title=title,
                  xlabel="Common west-to-east distance (m)", ylabel="Elevation (m a.s.l.)")
        group_distances = common_profile_distances(profiles[group])
        for (j, profile_index) in enumerate(group)
            spec, profile = specs[profile_index], profiles[profile_index]
            measured, interpolated = split_bedrock(profile)
            surface_z = [sample_bilinear(surface, x, y) for (x, y) in zip(profile.x, profile.y)]
            plot_distance = group_distances[j]
            order = sortperm(plot_distance)
            plot_distance = plot_distance[order]
            measured, interpolated, surface_z = measured[order], interpolated[order], surface_z[order]
            color = spec.year == 2009 ? :purple4 : :black
            gap_color = spec.year == 2009 ? :darkorange2 : :red3
            j == 1 && lines!(ax, plot_distance, surface_z; color=:dodgerblue3, linewidth=2, label="2023 surface")
            lines!(ax, plot_distance, measured; color, linewidth=2, label=profile_label(spec))
            if any(isfinite, interpolated)
                lines!(ax, plot_distance, interpolated; color=gap_color, linewidth=2.5, label=profile_label(spec) * " interpolated")
            end
        end
        ax.xgridcolor = (:grey, 0.22)
        ax.ygridcolor = (:grey, 0.22)
        axislegend(ax; position=:rb, framevisible=true, labelsize=10)
    end
    plot_location_map!(fig[1:nrows, 3], specs, profiles, gravity)
    Label(fig[0, 1:3], "Aletsch GPR bedrock cross-sections, 2023 surface, and survey locations", fontsize=22)
    save(joinpath(OUTPUT_DIR, "aletsch_gpr_bedrock_profiles.png"), fig; px_per_unit=1.5)
end

function read_gravimetry(path)
    rows = readlines(path)[2:end]
    x, y = Float64[], Float64[]
    for line in rows
        isempty(strip(line)) && continue
        fields = split(line, ',')
        length(fields) >= 5 || continue
        push!(x, parse(Float64, strip(fields[4])))
        push!(y, parse(Float64, strip(fields[5])))
    end
    return (; x, y)
end

"""Read a downsampled RGB map window; GDAL resamples directly into the buffer."""
function read_map_rgb(path; max_dimension=1800)
    ArchGDAL.read(path) do dataset
        nx, ny = ArchGDAL.width(dataset), ArchGDAL.height(dataset)
        scale = min(1.0, max_dimension / max(nx, ny))
        bx, by = max(1, round(Int, nx * scale)), max(1, round(Int, ny * scale))
        channels = Array{UInt8}(undef, bx, by, 3)
        for band in 1:3
            ArchGDAL.read!(dataset, view(channels, :, :, band), band, 0, 0, nx, ny)
        end
        gt = ArchGDAL.getgeotransform(dataset)
        rgb = [RGBf(channels[i, j, 1] / 255, channels[i, j, 2] / 255,
                    channels[i, j, 3] / 255) for i in 1:bx, j in by:-1:1]
        extent = (gt[1], gt[1] + nx * gt[2], gt[4] + ny * gt[6], gt[4])
        return (; rgb, extent)
    end
end

function additional_bedrock_profiles(excluded_specs)
    excluded = Set(normpath(profile_path(spec)) for spec in excluded_specs)
    paths = String[]
    for (root, _, files) in walkdir(GPR_DIR), file in files
        endswith(file, "_Bedrock.txt") || continue
        path = normpath(joinpath(root, file))
        path in excluded || push!(paths, path)
    end
    profiles = []
    for path in paths
        profile = read_profile(path)
        valid = isfinite.(profile.bed) .& (profile.bed .> 0)
        any(valid) || continue
        x = [valid[i] ? profile.x[i] : NaN for i in eachindex(valid)]
        y = [valid[i] ? profile.y[i] : NaN for i in eachindex(valid)]
        push!(profiles, (; x, y))
    end
    profiles
end

function plot_location_map!(layout, specs, profiles, gravity; profile_numbers=nothing)
    map_paths = [joinpath(ALETSCH, "swiss-map-raster25_2023_1269_krel_1.25_2056.tif"),
                 joinpath(ALETSCH, "swiss-map-raster25_2023_1249_krel_1.25_2056.tif")]
    maps = read_map_rgb.(map_paths)
    ax = Axis(layout, title="Profile and gravimetry locations (LV95)",
              xlabel="Easting (m)", ylabel="Northing (m)", aspect=DataAspect())
    for map in maps
        image!(ax, map.extent[1]..map.extent[2], map.extent[3]..map.extent[4], map.rgb)
    end
    allx = vcat([p.x for p in profiles]..., gravity.x)
    ally = vcat([p.y for p in profiles]..., gravity.y)
    xlimits = extrema(allx) .+ (-700, 700)
    ylimits = (minimum(ally) - 700, 1_150_000)
    first_additional = true
    for profile in additional_bedrock_profiles(specs)
        in_view = any(isfinite(x) && isfinite(y) &&
                      xlimits[1] <= x <= xlimits[2] && ylimits[1] <= y <= ylimits[2]
                      for (x, y) in zip(profile.x, profile.y))
        in_view || continue
        lines!(ax, profile.x, profile.y; color=:black, linewidth=1.2,
               label=first_additional ? "Bedrock reflection from airborne-GPR (2009–2011)" : nothing)
        first_additional = false
    end
    colors = Makie.wong_colors()
    for (i, (spec, profile)) in enumerate(zip(specs, profiles))
        color = colors[mod1(i, length(colors))]
        lines!(ax, profile.x, profile.y; color, linewidth=3)
        label_fractions = (0.32, 0.50, 0.68)
        label_offsets = ((0, 14), (0, -18), (18, 8), (-18, -12))
        label_index = clamp(round(Int, label_fractions[mod1(i, 3)] * length(profile.x)), 1, length(profile.x))
        label_offset = label_offsets[mod1(i, 4)]
        text!(ax, profile.x[label_index], profile.y[label_index]; text=profile_label(spec), color=:black,
              fontsize=13, font=:bold, align=(:center, :center), offset=label_offset)
        if profile_numbers !== nothing
            number_offset = (label_offset[1], label_offset[2] + (label_offset[2] >= 0 ? 27 : -27))
            scatter!(ax, [profile.x[label_index]], [profile.y[label_index]];
                     marker=:rect, markersize=Vec2f(48, 36), color=:white,
                     strokecolor=color, strokewidth=2, marker_offset=number_offset)
            text!(ax, profile.x[label_index], profile.y[label_index];
                  text=string(profile_numbers[i]), color, fontsize=25, font=:bold,
                  align=(:center, :center), offset=number_offset)
        end
    end
    scatter!(ax, gravity.x, gravity.y; marker=:diamond, color=:magenta,
             strokecolor=:white, strokewidth=1, markersize=11, label="Gravimetry (2026)")
    axislegend(ax; position=:rt, labelsize=17)
    xlims!(ax, xlimits)
    ylims!(ax, ylimits)
    return ax
end

function main()
    mkpath(OUTPUT_DIR)
    profiles = [read_profile(profile_path(spec)) for spec in PROFILE_SPECS]
    surface = read_raster(joinpath(ALETSCH, "alti3D_2m_Aletsch_2023.tif"))
    gravity = read_gravimetry(joinpath(ALETSCH, "profiles_gravi_aletsch", "gps_points_20260513.txt"))
    plot_cross_sections(PROFILE_SPECS, profiles, surface, gravity)
    println("Saved plots in ", OUTPUT_DIR)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
