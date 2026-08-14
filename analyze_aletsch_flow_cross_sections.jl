include("plot_aletsch_gpr_profiles.jl")

const FLOW_GROUPS = [[4, 5], [7], [3, 9], [1, 13]]
const GRAB_PROFILE_PATHS = [
    joinpath(@__DIR__, "output", "profiles", "Profile_1_B36-26_GlacierBed_Grab2021.csv"),
    joinpath(@__DIR__, "output", "profiles", "Profile_2_B36-26_GlacierBed_Grab2021.csv"),
]

function read_grab_profile(path)
    rows = readlines(path)[2:end]
    x, y, bed, surface = Float64[], Float64[], Float64[], Float64[]
    for row in rows
        fields = split(row, ',')
        length(fields) >= 3 || continue
        parsed_x, parsed_y = tryparse(Float64, fields[1]), tryparse(Float64, fields[2])
        (parsed_x === nothing || parsed_y === nothing) && continue
        parsed_bed = tryparse(Float64, fields[3])
        parsed_surface = length(fields) >= 6 ? tryparse(Float64, fields[6]) : nothing
        push!(x, parsed_x)
        push!(y, parsed_y)
        push!(bed, parsed_bed === nothing ? NaN : parsed_bed)
        push!(surface, parsed_surface === nothing ? NaN : parsed_surface)
    end
    if x[end] < x[1]
        reverse!(x); reverse!(y); reverse!(bed); reverse!(surface)
    end
    distance = zeros(length(x))
    for i in 2:length(x)
        distance[i] = distance[i-1] + hypot(x[i]-x[i-1], y[i]-y[i-1])
    end
    (; x, y, distance, bed, surface)
end

function interp_linear(x, y, query)
    good = findall(isfinite, y)
    length(good) >= 2 || return fill(NaN, length(query))
    xx, yy = x[good], y[good]
    order = sortperm(xx); xx, yy = xx[order], yy[order]
    out = fill(NaN, length(query))
    for (k, q) in enumerate(query)
        (xx[1] <= q <= xx[end]) || continue
        hi = searchsortedfirst(xx, q)
        hi == 1 && (out[k] = yy[1]; continue)
        hi > length(xx) && (out[k] = yy[end]; continue)
        lo = hi - 1
        t = (q - xx[lo]) / (xx[hi] - xx[lo])
        out[k] = (1-t)*yy[lo] + t*yy[hi]
    end
    out
end

function extend_profile_to_surface(distance, bed, surface)
    order = sortperm(distance)
    d, b, s = distance[order], bed[order], surface[order]
    valid = findall(isfinite.(b) .& (b .> 0))
    isempty(valid) && error("Profile has no valid bedrock picks")
    complete = fill(NaN, length(b)); measured = fill(NaN, length(b))
    complete[valid] .= b[valid]; measured[valid] .= b[valid]
    for (left, right) in zip(valid[1:end-1], valid[2:end])
        right > left + 1 && (complete[left:right] .= range(b[left], b[right], length=right-left+1))
    end
    center_position = argmin(b[valid])
    thickness = s .- b
    left_indices, right_indices = reverse(valid[1:center_position]), valid[center_position:end]
    left_cross = findfirst(i -> isfinite(thickness[i]) && thickness[i] <= 0, left_indices)
    right_cross = findfirst(i -> isfinite(thickness[i]) && thickness[i] <= 0, right_indices)
    if left_cross === nothing
        first_pick = first(valid)
        complete[1:first_pick] .= range(s[1], b[first_pick], length=first_pick)
    else
        outside = left_indices[left_cross]
        complete[1:outside-1] .= NaN; measured[1:outside-1] .= NaN
        complete[outside] = s[outside]; measured[outside] = s[outside]
    end
    if right_cross === nothing
        last_pick = last(valid)
        complete[last_pick:end] .= range(b[last_pick], s[end], length=length(b)-last_pick+1)
    else
        outside = right_indices[right_cross]
        complete[outside+1:end] .= NaN; measured[outside+1:end] .= NaN
        complete[outside] = s[outside]; measured[outside] = s[outside]
    end
    (; distance=d, complete, measured, surface=s)
end

function average_group(group, profiles, surface_raster; spacing=2.0)
    distances = common_profile_distances(profiles[group])
    sections = []
    for (j, profile_index) in enumerate(group)
        profile = profiles[profile_index]
        surface = [sample_bilinear(surface_raster, x, y) for (x,y) in zip(profile.x, profile.y)]
        push!(sections, extend_profile_to_surface(distances[j], profile.bed, surface))
    end
    xmin = minimum(minimum(s.distance) for s in sections)
    xmax = maximum(maximum(s.distance) for s in sections)
    grid = collect(xmin:spacing:xmax); grid[end] < xmax && push!(grid, xmax)
    complete = [interp_linear(s.distance, s.complete, grid) for s in sections]
    measured = [interp_linear(s.distance, s.measured, grid) for s in sections]
    surfaces = [interp_linear(s.distance, s.surface, grid) for s in sections]
    fit_member = findfirst(j -> PROFILE_SPECS[group[j]].year == 2011, eachindex(group))
    fit_member === nothing && error("Every fitted group must contain a 2011 profile")
    bed, surf, measured_bed = complete[fit_member], surfaces[fit_member], measured[fit_member]
    extended_bed = [isfinite(bed[k]) && !isfinite(measured_bed[k]) ? bed[k] : NaN for k in eachindex(grid)]
    (; distance=grid, bed, surface=surf, measured=measured_bed, extended=extended_bed, individual=measured)
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
                    simplex[i] = simplex[1] + 0.5(simplex[i]-simplex[1]); values[i] = objective(simplex[i])
                end
            end
        end
    end
    simplex[argmin(values)]
end

function fit_power_shape(section)
    good = findall(isfinite, section.bed); x, z = section.distance[good], section.bed[good]
    imin = argmin(z); width, relief = maximum(x)-minimum(x), maximum(z)-minimum(z)
    start = [x[imin], z[imin], log(max(relief/(0.5width)^2,1e-6)), log(2.0)]
    objective(p) = begin
        x0,z0,loga,logb=p; a,b=exp(loga),exp(logb)
        (0.35 <= b <= 6 && minimum(x)-0.2width <= x0 <= maximum(x)+0.2width) || return 1e30
        sum(abs2, z .- (z0 .+ a .* abs.(x.-x0).^b))/length(x)
    end
    p = nelder_mead(objective,start,[0.05width,0.05max(relief,1),0.3,0.15])
    x0,z0,loga,logb=p; a,b=exp(loga),exp(logb)
    (; a,b,x0,z0,rmse=sqrt(objective(p)),fitted=z0 .+ a .* abs.(section.distance.-x0).^b)
end

function fit_label_position(section)
    finite_bed = filter(isfinite, section.bed)
    finite_surface = filter(isfinite, section.surface)
    x = (minimum(section.distance) + maximum(section.distance)) / 2
    y = (minimum(finite_bed) + maximum(finite_surface)) / 2
    (; x, y)
end

function label_curve!(ax, distance, values, label, color; fraction=0.7, offset=(8,8))
    finite = findall(isfinite, values)
    isempty(finite) && return
    target_x = minimum(distance[finite]) + fraction * (maximum(distance[finite]) - minimum(distance[finite]))
    index = finite[argmin(abs.(distance[finite] .- target_x))]
    text!(ax,distance[index],values[index];text=label,color,fontsize=14,font=:bold,
          align=(:left,:bottom),offset)
end

function main_flow_analysis()
    mkpath(OUTPUT_DIR)
    profiles = [read_profile(profile_path(s)) for s in PROFILE_SPECS]
    surface = read_raster(joinpath(ALETSCH,"alti3D_2m_Aletsch_2023.tif"))
    sections = [average_group(g,profiles,surface) for g in FLOW_GROUPS]
    fits = fit_power_shape.(sections)
    grab_sections = read_grab_profile.(GRAB_PROFILE_PATHS)
    grab_fits = fit_power_shape.(grab_sections)
    fig = Figure(size=(2100,1450))
    for (i,(group,section,fit)) in enumerate(zip(FLOW_GROUPS,sections,fits))
        fit_index = first(filter(j -> PROFILE_SPECS[j].year == 2011, group))
        label=profile_label(PROFILE_SPECS[fit_index])
        ax=Axis(fig[cld(i,2),mod1(i,2)],title="$label   RMSE = $(round(fit.rmse,digits=2)) m",
                xlabel="West-to-east distance (m)",ylabel="Elevation (m a.s.l.)")
        lines!(ax,section.distance,section.surface;color=:dodgerblue3,linewidth=2,label="2023 surface")
        lines!(ax,section.distance,section.measured;color=:black,linewidth=2.8,label=length(group) > 1 ? "2011 bedrock used for fit" : "GPR bedrock used for fit")
        lines!(ax,section.distance,section.extended;color=:red3,linewidth=3,label="Interpolated to surface")
        lines!(ax,section.distance,fit.fitted;color=:forestgreen,linestyle=:dash,linewidth=2.5,label="power-law fit")
        fit_label = fit_label_position(section)
        text!(ax,fit_label.x,fit_label.y;text="b = $(round(fit.b,digits=2))",
              color=:forestgreen,fontsize=20,font=:bold,align=(:center,:center))
        if i == 1
            label_curve!(ax,section.distance,section.surface,"2023 glacier surface",:dodgerblue3;fraction=0.72)
            label_curve!(ax,section.distance,section.measured,"2011 GPR (used for fit)",:black;fraction=0.62)
            label_curve!(ax,section.distance,fit.fitted,"power-law fit",:forestgreen;fraction=0.48)
        end
    end
    for (i,(section,fit)) in enumerate(zip(grab_sections,grab_fits))
        ax=Axis(fig[3,i],title="Gravimetry Profile $i   RMSE = $(round(fit.rmse,digits=2)) m",
                xlabel="West-to-east distance (m)",ylabel="Elevation (m a.s.l.)")
        lines!(ax,section.distance,section.surface;color=:dodgerblue3,linewidth=2)
        lines!(ax,section.distance,section.bed;color=:deeppink2,linewidth=3)
        lines!(ax,section.distance,fit.fitted;color=:forestgreen,linestyle=:dash,linewidth=2.5)
        fit_label = fit_label_position(section)
        text!(ax,fit_label.x,fit_label.y;text="b = $(round(fit.b,digits=2))",
              color=:forestgreen,fontsize=20,font=:bold,align=(:center,:center))
        text!(ax,0.5,0.04;text="interpolated bedrock from Grab et al. 2021",space=:relative,
              color=:deeppink2,fontsize=14,font=:bold,align=(:center,:bottom))
    end
    fit_profile_indices = [first(filter(i -> PROFILE_SPECS[i].year == 2011, group)) for group in FLOW_GROUPS]
    gravity = read_gravimetry(joinpath(ALETSCH,"profiles_gravi_aletsch","gps_points_20260513.txt"))
    map_ax = plot_location_map!(fig[1:3,3], PROFILE_SPECS[fit_profile_indices],
                                profiles[fit_profile_indices], gravity)
    for (i,section) in enumerate(grab_sections)
        lines!(map_ax,section.x,section.y;color=:deeppink2,linewidth=4)
        label_index = cld(length(section.x),2)
        text!(map_ax,section.x[label_index],section.y[label_index];text="Profile $i",
              color=:deeppink4,fontsize=14,font=:bold,align=(:center,:bottom),offset=(0,10))
    end
    Label(fig[0,1:3],"Aletsch flow-transverse cross-sections: z = z₀ + a|x-x₀|ᵇ",fontsize=22)
    output=joinpath(OUTPUT_DIR,"aletsch_flow_cross_section_powerlaw_fits.png")
    save(output,fig;px_per_unit=1.5)
    csv=joinpath(OUTPUT_DIR,"aletsch_flow_cross_section_powerlaw_fits.csv")
    open(csv,"w") do io
        println(io,"panel_profiles,fit_profile,a,b,x0_m,z0_m,rmse_m")
        for (group,fit) in zip(FLOW_GROUPS,fits)
            fit_index = first(filter(i -> PROFILE_SPECS[i].year == 2011, group))
            println(io, join((join(profile_label.(PROFILE_SPECS[group]), "+"), profile_label(PROFILE_SPECS[fit_index]), fit.a, fit.b, fit.x0, fit.z0, fit.rmse), ","))
        end
        for (i,fit) in enumerate(grab_fits)
            println(io, join(("Gravimetry Profile $i", "Grab et al. (2021)", fit.a, fit.b, fit.x0, fit.z0, fit.rmse), ","))
        end
    end
    for (group,fit) in zip(FLOW_GROUPS,fits)
        println(join(profile_label.(PROFILE_SPECS[group])," + "),": b=",round(fit.b,digits=4),", RMSE=",round(fit.rmse,digits=2)," m")
    end
    for (i,fit) in enumerate(grab_fits)
        println("Gravimetry Profile $i (Grab et al. 2021): a=",fit.a,", b=",round(fit.b,digits=4),
                ", x0=",round(fit.x0,digits=2)," m, z0=",round(fit.z0,digits=2),
                " m, RMSE=",round(fit.rmse,digits=2)," m")
    end
    println("Saved ",output," and ",csv)
end

main_flow_analysis()
