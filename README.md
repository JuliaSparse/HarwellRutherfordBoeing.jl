# A [Julia](http://julialang.org) Reader for the [Harwell-Boeing](http://math.nist.gov/MatrixMarket/formats.html#hb) and [Rutherford-Boeing](https://www.cise.ufl.edu/research/sparse/matrices/DOC/rb.pdf) Formats

[![Build Status](https://travis-ci.org/JuliaSparse/HarwellRutherfordBoeing.jl.svg?branch=ci)](https://travis-ci.org/JuliaSparse/HarwellRutherfordBoeing.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/3qrjx53tfff2hnrl?svg=true)](https://ci.appveyor.com/project/dpo/harwellrutherfordboeing-jl)
[![Coverage Status](https://coveralls.io/repos/JuliaSparse/HarwellRutherfordBoeing.jl/badge.svg?branch=master)](https://coveralls.io/r/JuliaSparse/HarwellRutherfordBoeing.jl?branch=master)


## Installing

````JULIA
julia> Pkg.add("HarwellRutherfordBoeing")
julia> Pkg.test("HarwellRutherfordBoeing")
````

## Obtaining the Harwell-Boeing Collection

Retrieve the systems from

[ftp://ftp.cerfacs.fr/pub/algo/matrices/harwell_boeing](ftp://ftp.cerfacs.fr/pub/algo/matrices/harwell_boeing)

Build `hsplit.c` using `cc -o hsplit hsplit.c`. This tool may be used to split a data file into its constituents. This module will only read one set per file.

## Obtaining Matrices and Supplementary Data in Rutherford-Boeing Format

The best source is the [University of Florida Sparse Matrix Collection](http://www.cise.ufl.edu/research/sparse/matrices).

## Example

````JULIA
julia> using HarwellRutherfordBoeing
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

julia> rb = RutherfordBoeingData("aa3.rb")
Rutherford-Boeing data 1681 of type pra
825 rows, 8627 cols, 70806 nonzeros

julia> using PyPlot

julia> spy(rb.data, markersize=2)
````

This content is released under the [MIT](http://opensource.org/licenses/MIT) License.
<a rel="license" href="http://opensource.org/licenses/MIT">
<img alt="MIT license" height="40" src="http://upload.wikimedia.org/wikipedia/commons/c/c3/License_icon-mit.svg" /></a>
