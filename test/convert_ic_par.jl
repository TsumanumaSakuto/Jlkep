using Revise
using Jlkep
using Test

r = [0, -3299.603, -6599.206]
v = [9.645803, 0, 0]

println(Jlkep.ic2par(r, v, Jlkep.MU_EARTH))


par = Jlkep.Par(
	149598261.1504, # a, km
	0.01671123, # e, -
	-2.6720990848033185e-07, # i, rad
	0.0, # Ω, rad
	1.796601474049171, # ω, rad
	0.0, # E, rad
)

println(Jlkep.par2ic(par, Jlkep.MU_SUN))
