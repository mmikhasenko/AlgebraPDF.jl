
standardgauss(x,σ) = exp(-x^2/(2*σ^2)) / sqrt(2π) / σ
standardBW(x,m,Γ) = abs2(m*Γ/(m^2-x^2-1im*m*Γ))
# 
function aGauss(p, lims)
    μ, σ = keys(p)
    return pdf((x;p)->standardgauss.(x .- getproperty(p,μ), getproperty(p,σ)), p, lims)
end
function aBreitWigner(p, lims)
    m, Γ = keys(p)
    return pdf((x;p)->standardBW.(x,getproperty(p,m), getproperty(p,Γ)), p, lims)
end
function aExp(p, lims)
    α, = keys(p)
    return pdf((x;p)->exp.(x .*getproperty(p,α)), p, lims)
end