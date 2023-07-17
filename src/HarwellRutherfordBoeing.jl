module HarwellRutherfordBoeing

using Printf
using SparseArrays

export HarwellBoeingMatrix, RutherfordBoeingData
export write_rb, write_rb_aux

include("hrb_utils.jl")
include("harwell_boeing.jl")
include("rutherford_boeing.jl")
include("write_rb.jl")

end
