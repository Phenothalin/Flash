using Clapeyron

# model = SRK(["NITROGEN","CARBON DIOXIDE","METHANE",   "ETHANE",   "PROPANE",
#   "ISOBUTANE", "n-BUTANE", "ISOPENTANE",
#   "n-PENTANE", "n-HEXANE", "n-HEPTANE"])
model = SRK(["ETHANE",   "PROPANE"])
feed = [0.5,0.5]

compress_factor1 = compressibility_factor(model , 101325 , 298 , feed , ; phase=:v, threaded=true, vol0=nothing)
compress_factor2 = compressibility_factor(model , 101325 , 298 , feed , ; phase=:l, threaded=true, vol0=nothing)

println("vapor compress_factor= ",compress_factor1)
println("liquid compress_factor= ",compress_factor2)