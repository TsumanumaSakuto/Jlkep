using Test
using Revise
import Jlkep as jk
using GLMakie

# @testset "Jlkep.jl" begin
# r = [0, 6599.206,0]
# v = [9.045803, 0, 0]
# rf, vf = jk.propagate_lagrangian(r, v, 100, jk.MU_EARTH)
# println(rf, ",", vf)
# end
fig = Figure()

ax = Axis3(fig[1, 1], aspect = :data)

jk.orbit_plots.plot_planet!("Io",ax,2458258*365.25)
jk.orbit_plots.plot_planet!("Europa",ax,2458258*365.25)
display(fig)

# jk.orbit_plots.plot_planet("EARTH BARYCENTER")


