module AlgebraPDF
#
using Parameters
using QuadGK
using Optim
#
import Base.+
import Base.*

export fit_llh
include("fit.jl")

export pdf, npars, p2v, v2p
include("pdf.jl")

end # module
