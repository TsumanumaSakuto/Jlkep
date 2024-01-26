using Downloads

# Download generic kernels
println("genker_naif0012 = raw\"", Downloads.download("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"),"\"")
println("genker_gm_de431 = raw\"", Downloads.download("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc"),"\"")
println("genker_de440 = raw\"", Downloads.download("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets//de440.bsp"),"\"")


#Result
# genker_naif0012 = raw"C:\Users\Sakuto\AppData\Local\Temp\jl_8jAxylw7m5"
# genker_gm_de431 = raw"C:\Users\Sakuto\AppData\Local\Temp\jl_GRCB0zaPmh"
# genker_de440 = raw"C:\Users\Sakuto\AppData\Local\Temp\jl_e1o1WtST6x"