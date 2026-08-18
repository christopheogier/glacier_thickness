include("plot_gpr_sections_between_gravimetry.jl")

const COMPLETE_PROFILE_INDICES = [5, 7, 9, 11, 13]
const SURFACE_PATH = joinpath(ALETSCH, "alti3D_2m_Aletsch_2023.tif")
const GRAB_PATH = joinpath(ALETSCH, "B36-26_GlacierBed_Grab2021.tif")
const MILLAN_PATH = joinpath(ALETSCH, "THICKNESS_RGI-11_2021July09_on_grab2021_grid.tif")
const FRANK_PATH = joinpath(ALETSCH, "RGI60-11.01222_bedrock_franck2026_on_grab2021_grid.tif")
const FARINOTTI_PATH = joinpath(ALETSCH, "RGI60-11.01450_thickness_Farinotti2019.tif")

function comparison_sections(surface_raster, grab_raster)
    sections = Dict{String,Any}(section.label => section for section in
                                profiles_between_gravimetry(surface_raster, grab_raster))
    for index in COMPLETE_PROFILE_INDICES
        spec = PROFILE_SPECS[index]
        profile = read_profile(profile_path(spec))
        valid = isfinite.(profile.bed) .& (profile.bed .> 0)
        measured = [valid[i] ? profile.bed[i] : NaN for i in eachindex(valid)]
        surface = [sample_bilinear(surface_raster, x, y)
                   for (x, y) in zip(profile.x, profile.y)]
        northing = sum(profile.y[valid]) / count(valid)
        sections[profile_label(spec)] = (;
            label=profile_label(spec), x=profile.x, y=profile.y,
            distance=profile.distance, bed=measured, surface, northing)
    end
    sort!(collect(values(sections)); by=section -> -section.northing)
end

sample_section(raster, section) =
    [sample_bilinear(raster, x, y) for (x, y) in zip(section.x, section.y)]

function sample_section_utm32(raster, section)
    values = Float64[]
    ArchGDAL.importEPSG(2056; order=:trad) do source
        ArchGDAL.importEPSG(32632; order=:trad) do target
            ArchGDAL.createcoordtrans(source, target) do transform
                for (x, y) in zip(section.x, section.y)
                    ArchGDAL.createpoint(x, y) do point
                        ArchGDAL.transform!(point, transform)
                        utm_x, utm_y, _ = ArchGDAL.getpoint(point, 0)
                        push!(values, sample_bilinear(raster, utm_x, utm_y))
                    end
                end
            end
        end
    end
    values
end

function plot_comparison_map!(layout, sections)
    map_paths = [joinpath(ALETSCH, "swiss-map-raster25_2023_1269_krel_1.25_2056.tif"),
                 joinpath(ALETSCH, "swiss-map-raster25_2023_1249_krel_1.25_2056.tif")]
    maps = read_map_rgb.(map_paths)
    ax = Axis(layout, title="GPR transects", xlabel="Easting (m)",
              ylabel="Northing (m)", aspect=DataAspect())
    for map in maps
        image!(ax, map.extent[1]..map.extent[2], map.extent[3]..map.extent[4], map.rgb)
    end
    for section in sections
        valid = isfinite.(section.bed)
        map_x = [valid[i] ? section.x[i] : NaN for i in eachindex(valid)]
        map_y = [valid[i] ? section.y[i] : NaN for i in eachindex(valid)]
        lines!(ax, map_x, map_y; color=:black, linewidth=2.5)
        first_valid = findfirst(valid)
        first_valid === nothing && continue
        text!(ax, section.x[first_valid], section.y[first_valid]; text=section.label,
              color=:black, fontsize=13, font=:bold, align=(:left, :bottom))
    end
    all_x = vcat([section.x for section in sections]...)
    all_y = vcat([section.y for section in sections]...)
    xlims!(ax, minimum(all_x) - 500, maximum(all_x) + 500)
    ylims!(ax, minimum(all_y) - 500, maximum(all_y) + 500)
    ax
end

function main()
    mkpath(OUTPUT_DIR)
    surface_raster = read_raster(SURFACE_PATH)
    grab_raster = read_raster(GRAB_PATH)
    millan_raster = read_raster(MILLAN_PATH)
    frank_raster = read_raster(FRANK_PATH)
    farinotti_raster = read_raster(FARINOTTI_PATH)
    sections = comparison_sections(surface_raster, grab_raster)

    ncols = 2
    nrows = cld(length(sections), ncols)
    fig = Figure(size=(2350, 335 * nrows))
    axes = Axis[]
    for (i, section) in enumerate(sections)
        row, col = fldmod1(i, ncols)
        ax = Axis(fig[row, col], title=section.label,
                  xlabel=row == nrows ? "Distance along transect (m)" : "",
                  ylabel=col == 1 ? "Elevation (m a.s.l.)" : "")
        push!(axes, ax)
        distance = section.distance
        grab = sample_section(grab_raster, section)
        millan_thickness = sample_section(millan_raster, section)
        frank = sample_section(frank_raster, section)
        farinotti_thickness = sample_section_utm32(farinotti_raster, section)
        millan = section.surface .- millan_thickness
        farinotti = section.surface .- farinotti_thickness

        lines!(ax, distance, section.surface; color=:dodgerblue3, linewidth=2.3)
        lines!(ax, distance, section.bed; color=:black, linewidth=3)
        lines!(ax, distance, grab; color=:deeppink2, linewidth=2.2)
        lines!(ax, distance, millan; color=:forestgreen, linewidth=2.2)
        lines!(ax, distance, frank; color=:darkorange, linewidth=2.2)
        lines!(ax, distance, farinotti; color=:navy, linewidth=2.2)
    end
    linkyaxes!(axes...)

    plot_comparison_map!(fig[1:nrows, ncols + 1], sections)
    colsize!(fig.layout, ncols + 1, Fixed(850))
    Label(fig[0, 1:ncols],
          "Aletsch GPR cross-sections and published bed-elevation products",
          fontsize=24)
    elements = [LineElement(color=:dodgerblue3, linewidth=3),
                LineElement(color=:black, linewidth=3),
                LineElement(color=:deeppink2, linewidth=3),
                LineElement(color=:forestgreen, linewidth=3),
                LineElement(color=:darkorange, linewidth=3),
                LineElement(color=:navy, linewidth=3)]
    labels = ["2023 surface", "GPR transect", "Grab et al. (2021)",
              "Millan et al. (2022)", "Frank et al. (2026)",
              "Farinotti et al. (2019) consensus"]
    Legend(fig[1:nrows, ncols + 1], elements, labels; orientation=:vertical,
           labelsize=16, framevisible=true, tellheight=false, tellwidth=false,
           halign=:left, valign=:bottom, margin=(18, 18, 18, 18))

    output = joinpath(OUTPUT_DIR, "aletsch_all_models_cross_sections.png")
    save(output, fig; px_per_unit=1.5)
    println("Saved $output ($(length(sections)) unique GPR transects)")
end

main()
