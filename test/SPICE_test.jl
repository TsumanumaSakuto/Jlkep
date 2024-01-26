using SPICE: SPICE
using GLMakie


genker_naif0012 = raw"C:\Users\Sakuto\AppData\Local\Temp\jl_8jAxylw7m5"
genker_gm_de431 = raw"C:\Users\Sakuto\AppData\Local\Temp\jl_GRCB0zaPmh"
genker_de440 = raw"C:\Users\Sakuto\AppData\Local\Temp\jl_e1o1WtST6x"

# Load generic kernels
SPICE.furnsh(genker_naif0012) # Leap seconds kernel
SPICE.furnsh(genker_gm_de431) # Gravity Constant
SPICE.furnsh(genker_de440) # Planetary ephemeris kernel

# Gravity Constant
μ = SPICE.bodvrd("SUN", "GM", 1)[1]

# Parameter Setting
et0 = SPICE.str2et("2021/04/06 13:52:32 UTC") # Initial Epoch
etf = SPICE.str2et("2022/04/06 13:52:32 UTC") # Final Epoch
num_grid = 10000 # Number of Grid

# Pre-calculation
et_all = LinRange(et0, etf, num_grid)
x_earth = zeros(num_grid, 6)
x_mars = zeros(num_grid, 6)

# Iteration
for i ∈ 1:num_grid
	x_earth[i, :], _ = SPICE.spkez(399, et_all[i], "ECLIPJ2000", "NONE", 10)
	x_mars[i, :], _ = SPICE.spkez(4, et_all[i], "ECLIPJ2000", "NONE", 10)
end


fig = Figure()
ax = Axis3(fig[1, 1], aspect = :data)
lines!(ax, x_earth[:, 1], x_earth[:, 2], x_earth[:, 3])
lines!(ax, x_mars[:, 1], x_mars[:, 2], x_mars[:, 3])
# lines!(ax,r2)
display(fig)
