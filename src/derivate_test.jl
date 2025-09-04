using Clapeyron

#∂μ_∂x, ∂μ_∂n
function chemical_potential_derivatives(model::EoSModel, p, T, n; phase=:unknown)
    N = sum(n)  
    z = n ./ N   
    
    # ∂μ_∂x
    μ(z_vec) = Clapeyron.chemical_potential(model, p, T, z_vec; phase=phase)
    J = ForwardDiff.jacobian(μ, z)
    ∂μ_∂x = J

    # ∂μ_∂n
    ∂μ_∂n = (J .- sum(J .* z', dims=2)) ./ N
    
    return ∂μ_∂x, ∂μ_∂n
end

#∂ϕ_∂x, ∂ϕ_∂n
function fugacity_cofficient_derivatives(model::EoSModel, p, T, n; phase=:unknown)
    N = sum(n)  
    z = n ./ N   
    
    # ∂ϕ_∂x
    ϕ(z_vec) = Clapeyron.fugacity_coefficient(model, p, T, z_vec; phase=phase)
    J = ForwardDiff.jacobian(ϕ, z)
    ∂ϕ_∂x = J

    # ∂ϕ_∂n
    ∂ϕ_∂n = (J .- sum(J .* z', dims=2)) ./ N
    
    return ∂ϕ_∂x, ∂ϕ_∂n
end

model1 = SRK(["ethane", "propane"])
fugacity_coefficient1 = fugacity_coefficient(model1, 1832580, 295, [0.62227058444458, 0.37772941555542]; phase=:v, threaded=true, vol0=nothing)
fugacity_coefficient2 = fugacity_coefficient(model1, 1832580, 295, [0.369828657899892, 0.630171342100108]; phase=:l, threaded=true, vol0=nothing)
println("Fugacity Coefficient (Vapor Phase): ", fugacity_coefficient1)
println("Fugacity Coefficient (Liquid Phase): ", fugacity_coefficient2)

derivate_phicoeff_x1 = fugacity_cofficient_derivatives(model1,1832580 ,295,[0.62227058444458, 0.37772941555542];phase = :v)
derivate_phicoeff_x2 = fugacity_cofficient_derivatives(model1,1832580 ,295,[0.369828657899892, 0.630171342100108];phase = :l)
println("dϕdx (Vapor Phase): ", derivate_phicoeff_x1)
println("dϕdx (Liquid Phase): ", derivate_phicoeff_x2)
# println("dϕdn (Vapor Phase): ", derivate_phicoeff_x1[2])
# println("dϕdn (Liquid Phase): ", derivate_phicoeff_x2[2])

derivate_phicoeff_n1 = Clapeyron.∂lnϕ∂n∂P∂T(model1,1832580 ,295,[0.62227058444458, 0.37772941555542];phase = :v)
derivate_phicoeff_n2 = Clapeyron.∂lnϕ∂n∂P∂T(model1,1832580 ,295,[0.369828657899892, 0.630171342100108];phase = :l)
println("dlnϕdn (Vapor Phase): ", derivate_phicoeff_n1[2])
println("dlnϕdn (Liquid Phase): ", derivate_phicoeff_n2[2])