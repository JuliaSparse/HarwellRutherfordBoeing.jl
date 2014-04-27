# A [Julia](http://julialang.org) Reader for the [Harwell-Boeing Format](http://math.nist.gov/MatrixMarket/data/Harwell-Boeing)

## Obtaining the Harwell-Boeing Collection

Retrieve the systems from

[ftp://ftp.cerfacs.fr/pub/algo/matrices/harwell_boeing](ftp://ftp.cerfacs.fr/pub/algo/matrices/harwell_boeing)

Build `hsplit.c` using `cc -o hsplit hsplit.c`. This tool may be used to split a data file into its constituents. This module will only read one set per file.

## Example

````JULIA
julia> import("rb.jl")
julia> M = HarwellBoeingMatrix("well1850.rra")
Harwell-Boeing matrix WELL1850 of type RRA
1850 rows, 712 cols, 8758 nonzeros
1 right-hand sides, 0 guesses, 0 solutions

julia > M.matrix
1850x712 sparse matrix with 8758 Float64 entries:
[1   ,    1]  =  0.027735
[3   ,    1]  =  0.027735  # etc...

julia> M.rhs'
1x1850 Array{Float64,2}:
6.40676  0.58834  6.40279  0.595772  …  -3.30846  -2.91383  -2.91705
````

## Testing

The script `test_rb.jl` may be used to scan through the entire Harwell-Boeing collection and read each matrix.

[![GPLv3](http://www.gnu.org/graphics/gplv3-88x31.png)](http://www.gnu.org/licenses/gpl.html "GPLv3")
