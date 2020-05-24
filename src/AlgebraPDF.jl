module AlgebraPDF
#
using Parameters
using QuadGK
using Optim
# 
import Base: +, *

#
export fit_llh
import("fit.jl")

export pdf, npars, p2v, v2p
#
export +, *
import("pdf.jl")

end # module
