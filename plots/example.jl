using AlgebraPDF
using DelimitedFiles

using Plots
using LaTeXStrings
theme(:wong2, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    xlim=(:auto,:auto), ylim=(:auto,:auto),
    xlab=L"m(\Xi_c^+K^-)\,\,[\mathrm{GeV}]", ylab=L"\mathrm{Candidates}")
# 
# 
const xth = 2.96
const fitrange = (2.960, 3.180)

#########################################################################

data = readdlm(joinpath("plots", "data.txt"))[:,1]
Nev = size(data,1)

# amplitude for the signal
Φ2(x) = sqrt(x-xth)
Γ(x,m,Γ₀) = Γ₀*Φ2(x)/Φ2(m)
breitwigner(x,m,Γ₀) = 1/(m^2-x^2-1im*m*Γ(x,m,Γ₀)) # complex function

# create
bw = FunctionWithParameters((x;p)->breitwigner(x,p.m,p.Γ), (m=3.05,Γ=0.02))
# check
bw(3.02) # works on a number
bw([3.05,3.065]) # works on a vector
pars(bw) # gives namedtuple
bw(3.02; p=(m=3.05,Γ=0.02)) # works parameters
bw(3.02, [3.05, 0.02]) # 

# operators like `abs2`, `log` can be applied
bwsq = abs2(bw) # real function
plot(bwsq, fitrange...) # can plot, range is needed

# making PDF from the function just by normalizing
bwsq_norm = Normalized(bwsq, fitrange) # needs limits (ranges, support)
plot(bwsq_norm) # can plot


# define a type SimpleBW and method `func` for dispatch
struct SimpleBW{P} <: AbstractFunctionWithParameters
    p::P
end
import AlgebraPDF:func
function func(bw::SimpleBW, x::NumberOrTuple; p=pars(bw))
    m,Γ = (getproperty(p,s) for s in keys(bw.p))
    breitwigner(x, m, Γ)
end

# create four instances of the same PDF with different parameters
signalamps = []
push!(signalamps, SimpleBW(( m1=3.00, Γ1=6.5e-3 )))
push!(signalamps, SimpleBW(( m2=3.05, Γ2=2.3e-3 )))
push!(signalamps, SimpleBW(( m3=3.06, Γ3=4.0e-3 )))
push!(signalamps, SimpleBW(( m4=3.09, Γ4=9.9e-3 )))
@time plot(); plot!.(abs2.(signalamps), fitrange...); plot!()


# |A|^2 * phase_space
phasespace = FunctionWithParameters((x;p)->Φ2(x), ∅)
signalpdfs = [Normalized(abs2(A)*phasespace, fitrange) for A in signalamps]
@time plot(); plot!.(signalpdfs, fill=0, α=0.3); plot!()


#  - BW S-wave
@makefuntype SimpleBWg(x;p) =
    1/(p.m0^2 - x^2 - 1im*p.g^2*Φ2(x))
#
signalpdf0 = Normalized(
        abs2(SimpleBWg((m0=2.95, g=0.01)))*phasespace,
        fitrange)
plot!(signalpdf0, fill=0, α=0.3)
# 
backgrpdf = Normalized(phasespace, fitrange)
plot(backgrpdf)

#########################################################################
# full model
model0 = 
    signalpdf0 * (f0=0.05Nev,) +
    signalpdfs[1] * (f1=0.25Nev,) + signalpdfs[2] * (f2=0.15Nev,) + 
    signalpdfs[3] * (f3=0.15Nev,) + signalpdfs[4] * (f4=0.55Nev,) +
    backgrpdf * Ext(fb=0.3Nev,) # Ext - to be able to fix

model1 = fixpar(model0, :fb, 17.2)

@time fit_summary = fit(
    Extended(NegativeLogLikelihood(model1, data)))
# 
let bins = range(fitrange..., length=80)
    Ns = scaletobinneddata(bins)
    stephist(data; bins, lc=:black, lab="Data")
    for i in 1:5
        plot!(fit_summary.best_model[i], Ns, 300,
            fillto=0.0, α=0.4, lab="Signal$(i)")
    end
    plot!(fit_summary.best_model, Ns, 300, lc=:red, lab="Total fit")
    # 
    m = fit_summary.measurements
    for (i,(k,v)) in enumerate(zip(keys(m), m))
        annotate!(xth+0.01,25-1.2i,text("$k=$(string(v))",:left, 6))
    end
    plot!()
end
savefig(joinpath("plots", "sixcompfit.png"))
