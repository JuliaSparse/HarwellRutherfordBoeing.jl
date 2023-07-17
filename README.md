# A [Julia](http://julialang.org) reader for the [Harwell-Boeing](http://math.nist.gov/MatrixMarket/formats.html#hb) and [Rutherford-Boeing](https://www.cise.ufl.edu/research/sparse/matrices/DOC/rb.pdf) formats

[![build-gh][build-gh-img]][build-gh-url] [![build-cirrus][build-cirrus-img]][build-cirrus-url]

[build-gh-img]: https://github.com/JuliaSmoothOptimizers/Krylov.jl/workflows/CI/badge.svg?branch=main
[build-gh-url]: https://github.com/JuliaSmoothOptimizers/Krylov.jl/actions
[build-cirrus-img]: https://img.shields.io/cirrus/github/JuliaSmoothOptimizers/Krylov.jl?logo=Cirrus%20CI
[build-cirrus-url]: https://cirrus-ci.com/github/JuliaSmoothOptimizers/Krylov.jl

## Installing

```julia
julia> ]
pkg> add HarwellRutherfordBoeing
pkg> test HarwellRutherfordBoeing
```

## Obtaining the Harwell-Boeing Collection

Retrieve the systems from the [Harwell-Boeing Collection](https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/).

<!-- Build `hsplit.c` using `cc -o hsplit hsplit.c`. This tool may be used to split a data file into its constituents. This module will only read one set per file. -->

## Obtaining Matrices and Supplementary Data in Rutherford-Boeing Format

The best source is the [University of Florida Sparse Matrix Collection](https://sparse.tamu.edu/).

## Example

```julia
using HarwellRutherfordBoeing
M = HarwellBoeingMatrix("well1850.rra")
```
```julia
Harwell-Boeing matrix WELL1850 of type RRA
1850 rows, 712 cols, 8758 nonzeros
1 right-hand sides, 0 guesses, 0 solutions
```
```julia
M.matrix
```
```julia
1850×712 SparseMatrixCSC{Float64, Int64} with 8758 stored entries:
⎡⢧⠀⠀⠀⠀⢹⠀⠀⠘⡆⠀⢳⠀⠀⡇⎤
⎢⣼⡀⠀⠀⠀⢸⡀⠀⠀⡇⠀⣸⡀⠀⡇⎥
⎢⠈⣧⠀⠀⠀⠀⡇⠀⠀⢹⠀⠈⡇⠀⡇⎥
⎢⠀⠻⡄⠀⠀⠀⢧⠀⠀⢸⠀⠀⢧⠀⡇⎥
⎢⠈⠳⣧⠀⠀⠀⢸⡀⠀⠘⡆⠈⢻⡀⡇⎥
⎢⠀⠀⢸⡄⠀⠀⠀⡇⠀⠀⡇⠀⠈⡇⡇⎥
⎢⠀⠀⠠⣷⠀⠀⠀⢹⠀⠀⢹⠀⠀⣹⡇⎥
⎢⠀⠀⠰⢾⡀⠀⠀⠸⡄⠀⢸⠀⠀⠾⡇⎥
⎢⠀⠀⢀⣨⢧⠀⠀⠀⣇⠀⠘⡆⠀⣈⣇⎥
⎢⠀⠀⢰⠙⢒⠀⠀⠀⢸⠀⠀⡇⠀⡟⣻⎥
⎢⠀⢤⣾⡃⠘⠀⠀⠀⢸⠀⠀⡇⢴⣗⢸⎥
⎢⠀⡔⡄⠃⢸⠀⠀⠀⢸⠀⠀⡇⣴⠛⢸⎥
⎢⠀⠡⡇⠀⢸⠀⠀⠀⢸⠀⠀⡇⢹⠀⢸⎥
⎢⡴⣌⣁⠀⠘⠀⠀⠀⢸⠀⠀⡷⣌⡀⢻⎥
⎢⠀⠀⢈⣷⢘⡆⠀⠀⢸⠀⠀⡇⠀⣿⢸⎥
⎢⠀⢀⣀⢿⢮⡇⠀⠀⢸⠀⠀⡇⢀⣻⢾⎥
⎢⣭⣷⠀⠀⢸⡆⠀⠀⢸⠀⠀⣯⡇⠀⣿⎥
⎢⠰⣼⠆⠀⠘⠇⠀⠀⢸⠀⠀⡷⡴⠀⢸⎥
⎢⠉⠙⠓⢲⣤⠀⠀⠀⢸⠀⠀⡏⠛⢲⣼⎥
⎣⢤⠤⡴⠛⢉⠀⠀⠀⢸⠀⠀⣧⢤⠞⢹⎦

```
```julia
M.rhs'
```
```julia
1×1850 adjoint(::Matrix{Float64}) with eltype Float64:
 6.40676  0.58834  6.40279  0.595772  …  -3.30846  -2.91383  -2.91705
```
```julia
rb = RutherfordBoeingData("aa3.rb")
```
```julia
Rutherford-Boeing data 1681 of type pra
825 rows, 8627 cols, 70806 nonzeros
```
```julia
rb.data
```
```julia
825×8627 SparseMatrixCSC{Float64, Int64} with 70806 stored entries:
⎡⣿⣾⣿⣿⣿⣿⣻⣾⣿⣾⣿⣷⣿⣿⣿⣿⣿⣾⣿⣿⣿⣿⣿⣿⣷⣾⣿⣾⣿⣾⣟⣿⣿⣻⣿⣿⣿⣿⣷⡿⎤
⎣⢸⣿⠿⣿⢿⢿⣿⡿⣿⣿⣿⣿⣿⣽⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣿⣾⣿⣿⢿⣿⣿⡿⣿⣿⣿⣿⣿⣿⣾⣿⎦
```

This content is released under the [MIT](http://opensource.org/licenses/MIT) License.
<a rel="license" href="http://opensource.org/licenses/MIT">
<img alt="MIT license" height="40" src="http://upload.wikimedia.org/wikipedia/commons/c/c3/License_icon-mit.svg" /></a>
