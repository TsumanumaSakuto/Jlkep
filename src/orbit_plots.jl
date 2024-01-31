module orbit_plots

using GLMakie
include("core_module.jl")

function plot_planet!(planet,ax ,t₀ = 0, tf = nothing, num_grid = 60)
	if tf === nothing
		x₀, _ = SPICE.spkezr(planet, t₀, "ECLIPJ2000", "NONE", "SUN")
		per = ic2par(x₀[1:3], x₀[4:6], MU_SUN)
		tf = t₀ + calc_period(per.a, MU_SUN)
	end
	# Pre-calculation
	et_all = LinRange(t₀, tf, num_grid)
	x_planet = zeros(num_grid, 6)

	# Iteration
	for i ∈ 1:num_grid
		x_planet[i, :], _ = SPICE.spkezr(planet, et_all[i], "ECLIPJ2000", "NONE", "SUN")
	end
	lines!(ax, x_planet[:, 1], x_planet[:, 2], x_planet[:, 3])

end

end
