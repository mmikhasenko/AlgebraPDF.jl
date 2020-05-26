module AlgebraPDF
#
using Parameters
using QuadGK
using Optim
#
import Base.+
import Base.*

export fit_llh, llh
include("fit.jl")

export pdf, npars, p2v, v2p
export integral
include("pdf.jl")

export generate
include("generation.jl")

# export pol0, pol1, pol2, pol3, pol4
# export gaus, BW
# include("basic_lineshapes.jl")

end # module
