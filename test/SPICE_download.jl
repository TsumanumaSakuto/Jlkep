using Downloads

# Download generic kernels
println("genker_naif0012 = raw\"", Downloads.download("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/lsk/naif0012.tls"),"\"")
println("genker_gm_de431 = raw\"", Downloads.download("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/pck/gm_de431.tpc"),"\"")
println("genker_de440 = raw\"", Downloads.download("https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets//de440.bsp"),"\"")


#Result
# genker_naif0012 = "C:\Users\Sakuto\AppData\Local\Temp\jl_7SBgzQPU90"
# genker_gm_de431 = "C:\Users\Sakuto\AppData\Local\Temp\jl_nq1AlaR6hU"
# genker_de440 = "C:\Users\Sakuto\AppData\Local\Temp\jl_HcZqPCVML0"