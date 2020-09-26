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
export fix_parameters
include("pdf.jl")

export generate
include("generation.jl")

export conv_with_gauss, conv_with_gauss_sampling
include("convolution.jl")

end # module
