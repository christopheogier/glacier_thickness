using CairoMakie
using Rasters
import ArchGDAL
using Downloads
using Printf

"""
    plot_tif_files(tif_paths; output_path="output/plots/tif_overview.png")

Plot one or more GeoTIFF files in a shared CairoMakie figure and save it.
"""
function plot_tif_files(tif_paths; output_path="output/plots/tif_overview.png", colormap=Reverse(:ice))
    isempty(tif_paths) && error("No TIFF files were provided.")

    rasters = Raster.(tif_paths)
    nplots = length(rasters)
    ncols = min(2, nplots)
    nrows = cld(nplots, ncols)

    fig = Figure(size=(650 * ncols, 520 * nrows))

    for (i, (path, raster)) in enumerate(zip(tif_paths, rasters))
        row = cld(i, ncols)
        col = mod1(i, ncols)
        title = splitext(basename(path))[1]

        ax = Axis(fig[row, col], title=title, aspect=DataAspect())
        hm = heatmap!(ax, raster; colormap)
        Colorbar(fig[row, col + ncols], hm, label="value")
        tightlimits!(ax)
    end

    mkpath(dirname(output_path))
    save(output_path, fig)
    return fig
end

"""
    download_swissalti3d_tiles(csv_path; tiles_dir, force=false)

Download all GeoTIFF URLs listed in `csv_path` into `tiles_dir` and return the
local tile paths. Existing tiles are reused unless `force=true`.
"""
function download_swissalti3d_tiles(csv_path; tiles_dir, force=false)
    urls = filter(!isempty, strip.(readlines(csv_path)))
    isempty(urls) && error("No URLs found in $csv_path.")

    mkpath(tiles_dir)
    tile_paths = String[]

    for url in urls
        filename = basename(first(split(url, '?')))
        tile_path = joinpath(tiles_dir, filename)
        if force || !isfile(tile_path)
            @info "Downloading SwissALTI3D tile" url tile_path
            Downloads.download(String(url), tile_path)
        end
        push!(tile_paths, tile_path)
    end

    return tile_paths
end

"""
    build_swissalti3d_mosaic(csv_path, output_path; tiles_dir, force=false)

Download SwissALTI3D GeoTIFF tiles from `csv_path`, mosaic them, and save the
surface DEM to `output_path`. Existing outputs are reused unless `force=true`.
"""
function build_swissalti3d_mosaic(csv_path, output_path; tiles_dir=joinpath(dirname(output_path), "swissalti3d_tiles"), force=false)
    if isfile(output_path) && !force
        return output_path
    end

    tile_paths = download_swissalti3d_tiles(csv_path; tiles_dir, force=false)
    @info "Building SwissALTI3D mosaic" ntile=length(tile_paths) output_path

    rasters = Raster.(tile_paths)
    dem = mosaic(first, rasters...)

    mkpath(dirname(output_path))
    write(output_path, Rasters.GDALsource(), dem; force=true)
    return output_path
end

"""
    thickness_from_surface_and_bedrock(surface_path, bedrock_path, output_path)

Warp `surface_path` onto the bedrock raster grid, compute `surface - bedrock`,
save the thickness GeoTIFF, and return `output_path`.
"""
function thickness_from_surface_and_bedrock(surface_path, bedrock_path, output_path; method="bilinear", force=false, mask_path=nothing)
    if isfile(output_path) && !force
        return output_path
    end

    mktempdir() do tmpdir
        warped_surface_path = joinpath(tmpdir, "surface_on_bedrock_grid.tif")
        _warp_raster_to_match(surface_path, bedrock_path, warped_surface_path; method)

        surface_on_bedrock_grid = Raster(warped_surface_path)
        bedrock = Raster(bedrock_path)
        thickness = surface_on_bedrock_grid .- bedrock

        if mask_path !== nothing
            warped_mask_path = joinpath(tmpdir, "mask_on_bedrock_grid.tif")
            _warp_raster_to_match(mask_path, bedrock_path, warped_mask_path; method="near")
            mask = Raster(warped_mask_path)
            thickness = map((h, m) -> (!ismissing(h) && !ismissing(m) && isfinite(Float64(h)) && isfinite(Float64(m)) && Float64(m) > 0) ? h : missing, thickness, mask)
        end

        mkpath(dirname(output_path))
        write(output_path, Rasters.GDALsource(), thickness; force=true)
    end

    return output_path
end

function _warp_raster_to_match(source_path, target_path, output_path; method="bilinear")
    ArchGDAL.read(target_path) do target
        width = ArchGDAL.width(target)
        height = ArchGDAL.height(target)
        gt = ArchGDAL.getgeotransform(target)
        xmin = gt[1]
        xmax = gt[1] + gt[2] * width
        ymax = gt[4]
        ymin = gt[4] + gt[6] * height
        target_epsg = _target_epsg_code(target, gt)

        ArchGDAL.read(source_path) do source
            source_epsg = _target_epsg_code(source, ArchGDAL.getgeotransform(source))
            options = [
                "-s_srs", "EPSG:" * source_epsg,
                "-t_srs", "EPSG:" * target_epsg,
                "-te", string(xmin), string(ymin), string(xmax), string(ymax),
                "-ts", string(width), string(height),
                "-r", method,
                "-of", "GTiff",
            ]

            warped = ArchGDAL.unsafe_gdalwarp([source], options; dest=output_path)
            ArchGDAL.destroy(warped)
        end
    end

    return output_path
end

function _last_epsg_code(wkt)
    matches = collect(eachmatch(r"""AUTHORITY\["EPSG","(\d+)"\]""", wkt))
    isempty(matches) && error("Could not find an EPSG code in target raster projection.")
    return last(matches).captures[1]
end

"""
    plot_aletsch_thickness_overview(thickness_paths, labels; output_path)

Plot Aletsch ice-thickness products with one shared color scale, white at zero,
a glacier outline, and a marker on the maximum value of each panel.
"""
function plot_aletsch_thickness_overview(
    thickness_paths,
    labels;
    output_path="output/plots/aletsch_thick_overview.png",
    colormap=cgrad([:white, :lightcyan, :deepskyblue3, :royalblue4, :midnightblue]),
    outline_path=nothing,
)
    length(thickness_paths) == length(labels) || error("thickness_paths and labels must have the same length.")
    isempty(thickness_paths) && error("No thickness rasters were provided.")

    rasters = Raster.(thickness_paths)
    maxima = _raster_maximum_location.(thickness_paths)
    color_max = maximum(m -> m.value, maxima)
    colorrange = (0, color_max)

    fig = Figure(size=(440 * length(rasters) + 110, 450))
    heatmap_plot = nothing

    mktempdir() do tmpdir
        for (i, (raster, path, label, max_point)) in enumerate(zip(rasters, thickness_paths, labels, maxima))
            ax = Axis(fig[1, i], title=label, aspect=DataAspect())
            heatmap_plot = heatmap!(ax, raster; colormap, colorrange)
            tightlimits!(ax)

            if outline_path !== nothing
                warped_outline_path = joinpath(tmpdir, "outline_$i.tif")
                _warp_raster_to_match(outline_path, path, warped_outline_path; method="near")
                outline = map(v -> (!ismissing(v) && isfinite(Float64(v)) && Float64(v) > 0) ? 1.0 : 0.0, Raster(warped_outline_path))
                contour!(ax, outline; levels=[0.5], color=:black, linewidth=1.5)
            end

            scatter!(ax, [max_point.x], [max_point.y]; color=:red, strokecolor=:white, strokewidth=1.5, markersize=11)
            text!(
                ax,
                0.03,
                0.04;
                text="max. thick. = " * string(round(Int, max_point.value)) * " m",
                space=:relative,
                fontsize=13,
                color=:black,
                align=(:left, :bottom),
            )
        end
    end

    Colorbar(
        fig[1, length(rasters) + 1],
        heatmap_plot;
        label="Ice thickness (m)",
        ticks=_colorbar_ticks(color_max),
    )
    mkpath(dirname(output_path))
    save(output_path, fig)
    return fig
end

function _colorbar_ticks(color_max)
    top = round(Int, color_max)
    base_ticks = collect(0:100:(100 * floor(Int, color_max / 100)))
    ticks = unique(vcat(base_ticks, top))
    return ticks, string.(ticks)
end

function _raster_maximum_location(path)
    ArchGDAL.read(path) do dataset
        band = ArchGDAL.getband(dataset, 1)
        array = ArchGDAL.read(band)
        gt = ArchGDAL.getgeotransform(dataset)

        best_value = -Inf
        best_index = CartesianIndex(1, 1)
        for index in CartesianIndices(array)
            value = Float64(array[index])
            if isfinite(value) && value > best_value
                best_value = value
                best_index = index
            end
        end

        best_value == -Inf && error("Raster has no finite values: " * path)
        x = gt[1] + (best_index[1] - 0.5) * gt[2] + (best_index[2] - 0.5) * gt[3]
        y = gt[4] + (best_index[1] - 0.5) * gt[5] + (best_index[2] - 0.5) * gt[6]
        return (; value=best_value, x=x, y=y)
    end
end

function _raster_maximum_location_on_positive_mask(value_path, mask_path)
    ArchGDAL.read(value_path) do value_dataset
        ArchGDAL.read(mask_path) do mask_dataset
            value_array = ArchGDAL.read(ArchGDAL.getband(value_dataset, 1))
            mask_array = ArchGDAL.read(ArchGDAL.getband(mask_dataset, 1))
            size(value_array) == size(mask_array) || error("Value and mask rasters must be on the same grid.")
            gt = ArchGDAL.getgeotransform(value_dataset)

            best_value = -Inf
            best_index = CartesianIndex(1, 1)
            for index in CartesianIndices(value_array)
                value = Float64(value_array[index])
                mask_value = Float64(mask_array[index])
                if isfinite(value) && isfinite(mask_value) && mask_value > 0 && value > best_value
                    best_value = value
                    best_index = index
                end
            end

            best_value == -Inf && error("Raster has no finite values inside positive mask: " * value_path)
            x = gt[1] + (best_index[1] - 0.5) * gt[2] + (best_index[2] - 0.5) * gt[3]
            y = gt[4] + (best_index[1] - 0.5) * gt[5] + (best_index[2] - 0.5) * gt[6]
            return (; value=best_value, x=x, y=y)
        end
    end
end

function _annotate_maximum_thickness!(ax, max_point)
    scatter!(ax, [max_point.x], [max_point.y]; color=:red, strokecolor=:white, strokewidth=1.5, markersize=11)
    text!(
        ax,
        0.03,
        0.04;
        text="max. thick. = " * string(round(Int, max_point.value)) * " m",
        space=:relative,
        fontsize=13,
        color=:black,
        align=(:left, :bottom),
    )
end

function _finite_raster_maximum(raster)
    values = Float64[]
    for value in parent(raster)
        if !ismissing(value)
            numeric_value = Float64(value)
            isfinite(numeric_value) && push!(values, numeric_value)
        end
    end
    isempty(values) && error("Raster has no finite values.")
    return maximum(values)
end

function _default_data_dir_for_region(region)
    data_root = joinpath(dirname(@__DIR__), "data")
    region == "11" && return joinpath(data_root, "Aletsch")
    region == "14" && return joinpath(data_root, "data_biafo")
    return data_root
end

"""
    extract_glacier_from_Frank_et_al(rgi_id; data_dir=nothing)

Extract Frank et al. thickness and bedrock GeoTIFFs for one RGI6 glacier ID from
`Thk.zip` and `Topg.zip`, then save them in `data_dir`.
"""
function extract_glacier_from_Frank_et_al(rgi_id; data_dir=nothing)
    id = rgi_id isa Real ? @sprintf("%.5f", rgi_id) : string(rgi_id)
    occursin(r"^\d{1,2}\.\d{5}$", id) ||
        error("Expected an RGI6 glacier ID like \"14.00005\".")

    region = lpad(first(split(id, ".")), 2, "0")
    data_dir = data_dir === nothing ? _default_data_dir_for_region(region) : data_dir
    full_id = "RGI60-$id"

    thk_zip = joinpath(data_dir, "Thk.zip")
    topg_zip = joinpath(data_dir, "Topg.zip")
    isfile(thk_zip) || error("Missing thickness archive: $thk_zip")
    isfile(topg_zip) || error("Missing bedrock archive: $topg_zip")

    outputs = Dict{Symbol,String}()

    mktempdir() do tmpdir
        outputs[:thickness] = _extract_frank_tif(
            thk_zip,
            "RGI-$region" * "_thk.zip",
            "thk_final/$full_id" * "_thk.tif",
            joinpath(data_dir, full_id * "_thickness_frank2026.tif"),
            tmpdir,
        )

        outputs[:bedrock] = _extract_frank_tif(
            topg_zip,
            "RGI-$region" * "_topg.zip",
            "topg_final/$full_id" * "_topg.tif",
            joinpath(data_dir, full_id * "_bedrock_franck2026.tif"),
            tmpdir,
        )
    end

    return outputs
end

function _extract_frank_tif(outer_zip, region_zip_name, tif_name, output_path, tmpdir)
    region_zip_path = joinpath(tmpdir, region_zip_name)

    open(region_zip_path, "w") do io
        run(pipeline(`unzip -p $outer_zip $region_zip_name`, stdout=io))
    end

    inner_files = split(read(`unzip -Z1 $region_zip_path`, String), '\n')
    tif_name in inner_files ||
        error("Could not find $tif_name inside $region_zip_name.")

    mkpath(dirname(output_path))
    open(output_path, "w") do io
        run(pipeline(`unzip -p $region_zip_path $tif_name`, stdout=io))
    end

    return output_path
end

function _target_epsg_code(target, gt)
    matches = collect(eachmatch(r"""AUTHORITY\["EPSG","(\d+)"\]""", ArchGDAL.getproj(target)))
    !isempty(matches) && return last(matches).captures[1]

    # ASCII grids often lack CRS metadata; infer from the coordinate range used here.
    return gt[1] > 1_000_000 ? "2056" : "32632"
end

"""
    plot_aletsch_bedrock_overview(bedrock_paths, labels; output_path)

Plot Aletsch bedrock/elevation products with one shared color scale and an
optional glacier outline.
"""
function plot_aletsch_bedrock_overview(
    bedrock_paths,
    labels;
    output_path="output/plots/aletsch_bedrock_overview.png",
    colormap=:oleron,
    outline_path=nothing,
)
    length(bedrock_paths) == length(labels) || error("bedrock_paths and labels must have the same length.")
    isempty(bedrock_paths) && error("No bedrock rasters were provided.")

    rasters = Raster.(bedrock_paths)
    ranges = _finite_raster_range.(rasters)
    colorrange = (minimum(first, ranges), maximum(last, ranges))

    fig = Figure(size=(480 * length(rasters) + 110, 470))
    heatmap_plot = nothing

    mktempdir() do tmpdir
        for (i, (raster, path, label)) in enumerate(zip(rasters, bedrock_paths, labels))
            ax = Axis(fig[1, i], title=label, aspect=DataAspect())
            heatmap_plot = heatmap!(ax, raster; colormap, colorrange)
            tightlimits!(ax)

            if outline_path !== nothing
                warped_outline_path = joinpath(tmpdir, "bedrock_outline_$i.tif")
                _warp_raster_to_match(outline_path, path, warped_outline_path; method="near")
                outline = map(v -> (!ismissing(v) && isfinite(Float64(v)) && Float64(v) > 0) ? 1.0 : 0.0, Raster(warped_outline_path))
                contour!(ax, outline; levels=[0.5], color=:black, linewidth=1.5)
            end
        end
    end

    Colorbar(fig[1, length(rasters) + 1], heatmap_plot; label="Bedrock elevation (m)")
    mkpath(dirname(output_path))
    save(output_path, fig)
    return fig
end

function _finite_raster_range(raster)
    min_value = Inf
    max_value = -Inf
    for value in parent(raster)
        if !ismissing(value)
            numeric_value = Float64(value)
            if isfinite(numeric_value)
                min_value = min(min_value, numeric_value)
                max_value = max(max_value, numeric_value)
            end
        end
    end
    min_value == Inf && error("Raster has no finite values.")
    return (min_value, max_value)
end

"""
    bedrock_difference_raster(grab_path, frank_path, output_path; force=false)

Save `Grab 2021 - Frank 2026` bedrock difference on an already-matched grid.
"""
function bedrock_difference_raster(grab_path, frank_path, output_path; force=false)
    if isfile(output_path) && !force
        return output_path
    end

    grab = Raster(grab_path)
    frank = Raster(frank_path)
    difference = map((a, b) -> (!ismissing(a) && !ismissing(b) && isfinite(Float64(a)) && isfinite(Float64(b))) ? Float32(Float64(a) - Float64(b)) : Float32(NaN), grab, frank)

    mkpath(dirname(output_path))
    write(output_path, Rasters.GDALsource(), difference; force=true)
    return output_path
end

"""
    plot_aletsch_bedrock_difference(difference_path; output_path, outline_path)

Plot one bedrock difference raster using a red colorscale.
"""
function plot_aletsch_bedrock_difference(
    difference_path;
    output_path="output/plots/aletsch_bedrock_difference.png",
    outline_path=nothing,
    colormap=cgrad([:white, :mistyrose, :red2, :darkred]),
)
    raster = Raster(difference_path)
    colorrange = _finite_raster_range(raster)

    fig = Figure(size=(620, 520))
    ax = Axis(fig[1, 1], title="Grab 2021 - Frank 2026", aspect=DataAspect())
    hm = heatmap!(ax, raster; colormap, colorrange)
    tightlimits!(ax)

    if outline_path !== nothing
        mktempdir() do tmpdir
            warped_outline_path = joinpath(tmpdir, "bedrock_difference_outline.tif")
            _warp_raster_to_match(outline_path, difference_path, warped_outline_path; method="near")
            outline = map(v -> (!ismissing(v) && isfinite(Float64(v)) && Float64(v) > 0) ? 1.0 : 0.0, Raster(warped_outline_path))
            contour!(ax, outline; levels=[0.5], color=:black, linewidth=1.5)
        end
    end

    Colorbar(fig[1, 2], hm; label="Bedrock difference (m)")
    mkpath(dirname(output_path))
    save(output_path, fig)
    return fig
end

"""
    absolute_difference_raster(reference_path, comparison_path, output_path; force=false)

Save `abs(reference - comparison)` on an already-matched grid.
"""
function absolute_difference_raster(reference_path, comparison_path, output_path; force=false)
    if isfile(output_path) && !force
        return output_path
    end

    reference = Raster(reference_path)
    comparison = Raster(comparison_path)
    difference = map(
        (a, b) -> (!ismissing(a) && !ismissing(b) && isfinite(Float64(a)) && isfinite(Float64(b)))  ? Float64(a) - Float64(b) : missing,
        reference,
        comparison,
    )

    mkpath(dirname(output_path))
    write(output_path, Rasters.GDALsource(), difference; force=true, missingval=-9999f0)
    return output_path
end

"""
    plot_biafo_thickness_comparison(consensus_path, frank_path, difference_path; output_path)

Plot Biafo consensus thickness, Frank thickness on the consensus grid, and the
absolute difference between them.
"""
function plot_biafo_thickness_comparison(
    consensus_path,
    frank_path,
    difference_path;
    output_path="output/plots/biafo_thickness_comparison.png",
    thickness_colormap=cgrad([:white, :lightcyan, :deepskyblue3, :royalblue4, :midnightblue]),
    difference_colormap=cgrad([:navy, :white, :darkred]),
    outline_path=nothing,
)
    consensus = Raster(consensus_path)
    frank = Raster(frank_path)
    difference = Raster(difference_path)

    thickness_max = max(_finite_raster_range(consensus)[2], _finite_raster_range(frank)[2])
    difference_range = _finite_raster_range(difference)
    difference_max = maximum(abs, difference_range)
    consensus_max_point = _raster_maximum_location_on_positive_mask(consensus_path, consensus_path)
    frank_max_point = _raster_maximum_location_on_positive_mask(frank_path, consensus_path)

    fig = Figure(size=(1380, 440))

    ax1 = Axis(fig[1, 1], title="Consensus 2019", aspect=DataAspect())
    hm1 = heatmap!(ax1, consensus; colormap=thickness_colormap, colorrange=(0, thickness_max))
    tightlimits!(ax1)

    if outline_path !== nothing
        _plot_outline!(ax1, outline_path, consensus_path)
    end

    ax2 = Axis(fig[1, 2], title="Frank 2026", aspect=DataAspect())
    heatmap!(ax2, frank; colormap=thickness_colormap, colorrange=(0, thickness_max))
    tightlimits!(ax2)

    if outline_path !== nothing
        _plot_outline!(ax2, outline_path, frank_path)
    end

    ax3 = Axis(fig[1, 3], title="Consensus 2019 - Frank 2026", aspect=DataAspect())
    hm3 = heatmap!(ax3, difference; colormap=difference_colormap, colorrange=(-difference_max, difference_max))
    tightlimits!(ax3)

    if outline_path !== nothing
        _plot_outline!(ax3, outline_path, difference_path)
    end

    Colorbar(fig[1, 4], hm1; label="Ice thickness (m)")
    Colorbar(fig[1, 5], hm3; label="Consensus - Frank (m)")

    mkpath(dirname(output_path))
    save(output_path, fig)
    return fig
end

"""
    mask_raster_by_positive_reference(input_path, reference_path, output_path; force=false)

Set cells to missing wherever `reference_path` is missing, non-finite, or <= 0.
The input and reference rasters must already be on the same grid.
"""
function mask_raster_by_positive_reference(input_path, reference_path, output_path; force=false)
    if isfile(output_path) && !force
        return output_path
    end

    input = Raster(input_path)
    reference = Raster(reference_path)
    masked = map(
        (value, mask_value) -> (!ismissing(value) && !ismissing(mask_value) && isfinite(Float64(value)) && isfinite(Float64(mask_value)) && Float64(mask_value) > 0) ? value : missing,
        input,
        reference,
    )

    mkpath(dirname(output_path))
    write(output_path, Rasters.GDALsource(), masked; force=true, missingval=-9999f0)
    return output_path
end

function _plot_outline!(ax, outline_path, target_path)
    mktempdir() do tmpdir
        warped_outline_path = joinpath(tmpdir, "outline.tif")
        _warp_raster_to_match(outline_path, target_path, warped_outline_path; method="near")
        outline = map(v -> (!ismissing(v) && isfinite(Float64(v)) && Float64(v) > 0) ? 1.0 : 0.0, Raster(warped_outline_path))
        contour!(ax, outline; levels=[0.5], color=:black, linewidth=1.5)
    end
end

"""
    plot_biafo_thickness_comparison(consensus_path, frank_path; output_path)

Plot Biafo consensus thickness, Frank thickness on the consensus grid, and the
signed difference. Cells outside the positive consensus footprint are hidden.
"""
function plot_biafo_thickness_comparison(
    consensus_path,
    frank_path;
    output_path="output/plots/biafo_thickness_comparison.png",
    thickness_colormap=cgrad([:white, :lightcyan, :deepskyblue3, :royalblue4, :midnightblue]),
    difference_colormap=cgrad([:navy, :white, :darkred]),
)
    raw_consensus = Raster(consensus_path)
    raw_frank = Raster(frank_path)

    consensus = map(v -> (!ismissing(v) && isfinite(Float64(v)) && Float64(v) > 0) ? v : missing, raw_consensus)
    frank = map(
        (value, mask_value) -> (!ismissing(value) && !ismissing(mask_value) && isfinite(Float64(value)) && isfinite(Float64(mask_value)) && Float64(mask_value) > 0) ? value : missing,
        raw_frank,
        raw_consensus,
    )
    difference = map(
        (a, b) -> (!ismissing(a) && !ismissing(b) && isfinite(Float64(a)) && isfinite(Float64(b))) ? Float64(a) - Float64(b) : missing,
        consensus,
        frank,
    )
    outline = map(v -> (!ismissing(v) && isfinite(Float64(v)) && Float64(v) > 0) ? 1.0 : 0.0, raw_consensus)

    thickness_max = max(_finite_raster_range(consensus)[2], _finite_raster_range(frank)[2])
    difference_range = _finite_raster_range(difference)
    difference_max = maximum(abs, difference_range)
    consensus_max_point = _raster_maximum_location_on_positive_mask(consensus_path, consensus_path)
    frank_max_point = _raster_maximum_location_on_positive_mask(frank_path, consensus_path)

    fig = Figure(size=(1380, 440))

    ax1 = Axis(fig[1, 1], title="Consensus 2019", aspect=DataAspect())
    hm1 = heatmap!(ax1, consensus; colormap=thickness_colormap, colorrange=(0, thickness_max))
    contour!(ax1, outline; levels=[0.5], color=:black, linewidth=1.5)
    tightlimits!(ax1)
    _annotate_maximum_thickness!(ax1, consensus_max_point)

    ax2 = Axis(fig[1, 2], title="Frank 2026", aspect=DataAspect())
    heatmap!(ax2, frank; colormap=thickness_colormap, colorrange=(0, thickness_max))
    contour!(ax2, outline; levels=[0.5], color=:black, linewidth=1.5)
    tightlimits!(ax2)
    _annotate_maximum_thickness!(ax2, frank_max_point)

    ax3 = Axis(fig[1, 3], title="Consensus 2019 - Frank 2026", aspect=DataAspect())
    hm3 = heatmap!(ax3, difference; colormap=difference_colormap, colorrange=(-difference_max, difference_max))
    contour!(ax3, outline; levels=[0.5], color=:black, linewidth=1.5)
    tightlimits!(ax3)

    Colorbar(fig[1, 4], hm1; label="Ice thickness (m)")
    Colorbar(fig[1, 5], hm3; label="Consensus - Frank (m)")

    mkpath(dirname(output_path))
    save(output_path, fig)
    return fig
end

"""
    plot_biafo_thickness_comparison_with_milan(consensus_path, milan_path, frank_path; output_path)

Plot Biafo consensus, Milan 2022, Frank 2026, and the signed Consensus - Frank difference.
Cells outside the positive consensus footprint are hidden.
"""
function plot_biafo_thickness_comparison_with_milan(
    consensus_path,
    milan_path,
    frank_path;
    output_path="output/plots/biafo_thickness_comparison.png",
    thickness_colormap=cgrad([:white, :lightcyan, :deepskyblue3, :royalblue4, :midnightblue]),
    difference_colormap=cgrad([:navy, :white, :darkred]),
)
    raw_consensus = Raster(consensus_path)
    raw_milan = Raster(milan_path)
    raw_frank = Raster(frank_path)

    consensus = map(v -> (!ismissing(v) && isfinite(Float64(v)) && Float64(v) > 0) ? v : missing, raw_consensus)
    milan = map(
        (value, mask_value) -> (!ismissing(value) && !ismissing(mask_value) && isfinite(Float64(value)) && isfinite(Float64(mask_value)) && Float64(mask_value) > 0) ? value : missing,
        raw_milan,
        raw_consensus,
    )
    frank = map(
        (value, mask_value) -> (!ismissing(value) && !ismissing(mask_value) && isfinite(Float64(value)) && isfinite(Float64(mask_value)) && Float64(mask_value) > 0) ? value : missing,
        raw_frank,
        raw_consensus,
    )
    difference = map(
        (a, b) -> (!ismissing(a) && !ismissing(b) && isfinite(Float64(a)) && isfinite(Float64(b))) ? Float64(a) - Float64(b) : missing,
        consensus,
        frank,
    )
    outline = map(v -> (!ismissing(v) && isfinite(Float64(v)) && Float64(v) > 0) ? 1.0 : 0.0, raw_consensus)

    thickness_max = maximum([_finite_raster_range(consensus)[2], _finite_raster_range(milan)[2], _finite_raster_range(frank)[2]])
    difference_range = _finite_raster_range(difference)
    difference_max = maximum(abs, difference_range)
    max_points = [
        _raster_maximum_location_on_positive_mask(consensus_path, consensus_path),
        _raster_maximum_location_on_positive_mask(milan_path, consensus_path),
        _raster_maximum_location_on_positive_mask(frank_path, consensus_path),
    ]

    fig = Figure(size=(1780, 440))

    thickness_rasters = [consensus, milan, frank]
    titles = ["Consensus 2019", "Milan 2022", "Frank 2026"]
    heatmap_plot = nothing
    for i in eachindex(thickness_rasters)
        ax = Axis(fig[1, i], title=titles[i], aspect=DataAspect())
        heatmap_plot = heatmap!(ax, thickness_rasters[i]; colormap=thickness_colormap, colorrange=(0, thickness_max))
        contour!(ax, outline; levels=[0.5], color=:black, linewidth=1.5)
        tightlimits!(ax)
        _annotate_maximum_thickness!(ax, max_points[i])
    end

    ax4 = Axis(fig[1, 4], title="Consensus 2019 - Frank 2026", aspect=DataAspect())
    hm4 = heatmap!(ax4, difference; colormap=difference_colormap, colorrange=(-difference_max, difference_max))
    contour!(ax4, outline; levels=[0.5], color=:black, linewidth=1.5)
    tightlimits!(ax4)

    Colorbar(fig[1, 5], heatmap_plot; label="Ice thickness (m)")
    Colorbar(fig[1, 6], hm4; label="Consensus - Frank (m)")

    mkpath(dirname(output_path))
    save(output_path, fig)
    return fig
end
