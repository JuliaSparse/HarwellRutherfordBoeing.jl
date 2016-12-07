# Writing data to file involves writing full lines and partial lines.
# Because @sprintf will not accept computed formats, we predefine
# all the variants of @sprintf that we may need, both for float
# and integer data. That's a truckload of variants, but the alternative
# is very slow.


global const real_fortfmt = ["(8E10.1E3)", "(7E11.2E3)", "(6E12.3E3)", "(6E13.4E3)",
                             "(5E14.5E3)", "(5E15.6E3)", "(5E16.7E3)", "(4E17.8E3)",
                             "(4E18.9E3)", "(4E19.10E3)", "(4E20.11E3)", "(3E21.12E3)",
                             "(3E22.13E3)", "(3E23.14E3)", "(3E24.15E3)", "(3E25.16E3)"]
global const real_fmts = ["%10.1e", "%11.2e", "%12.3e", "%13.4e",
                          "%14.5e", "%15.6e", "%16.7e", "%17.8e",
                          "%18.9e", "%19.10e", "%20.11e", "%21.12e",
                          "%22.13e", "%23.14e", "%24.15e", "%25.16e"]
global const real_lens = 10:25
global const real_nitems = [8, 7, 6, 6, 5, 5, 5, 4, 4, 4, 4, 3, 3, 3, 3, 3]

# define several variants of @sprintf per format
for (fmt, len, nitem) in zip(real_fmts, real_lens, real_nitems)
  for k = 1 : nitem  # need to be able to print 1, 2, ..., nitems items for each format

    fmt1 = fmt^k * "\n"
    @eval begin
      function $(Symbol("sprintf_real_$(k)_$len"))(data)
        @sprintf($fmt1, data...)
      end
    end

  end
end

global const int_lens = vcat(collect(2:11), [20])
global const int_fmts = ["%$(k)d" for k in int_lens]
global const int_nitems = [40, 26, 20, 16, 13, 11, 10, 8, 8, 7, 4]

# define several variants of @sprintf per format
for (fmt, len, nitem) in zip(int_fmts, int_lens, int_nitems)
  for k = 1 : nitem  # need to be able to print 1, 2, ..., nitems items for each format

    fmt1 = fmt^k * "\n"
    @eval begin
      function $(Symbol("sprintf_int_$(k)_$len"))(data)
        @sprintf($fmt1, data...)
      end
    end

  end
end


"Return appropriate format for integer arrays."
function get_int_fmt(n :: Int)
  ndigits = max(2, Int(1 + floor(log10(n))))
  if ndigits <= 10
    return (ndigits,
            @sprintf("%dI%d", int_nitems[ndigits-1], ndigits),
            int_nitems[ndigits-1])
  end
  return (int_lens[11],
          @sprintf("%dI%d", int_nitems[11], int_lens[11]),
          int_nitems[11])
end

"Return appropriate format for real array (1 <= p <= 16)."
function get_real_fmt(p :: Int)
  return (real_lens[p], real_fortfmt[p], real_nitems[p])
end


function fortran_write_line{T <: AbstractFloat}(io :: IOStream, data :: Vector{T}, sprintf_variant :: Function)
  length(data) > 0 || return
  print(io, replace(sprintf_variant(data), r"[eE]", "D"))
end


function fortran_write_line{T <: Integer}(io :: IOStream, data :: Vector{T}, sprintf_variant :: Function)
  length(data) > 0 || return
  print(io, sprintf_variant(data))
end


function fortran_write_array{T <: Integer}(io :: IOStream, data :: Vector{T}, chunk_size :: Int, ndigits :: Int)
  nelts = length(data)
  nlines, leftover = fldmod(nelts, chunk_size)
  sprintf_full_line = eval(Symbol("sprintf_int_$(chunk_size)_$ndigits"))
  for k = 0:nlines-1
    fortran_write_line(io, data[k * chunk_size + 1 : (k + 1) * chunk_size], sprintf_full_line)
  end
  if leftover > 0
    sprintf_partial_line = eval(Symbol("sprintf_int_$(leftover)_$ndigits"))
    fortran_write_line(io, data[nlines * chunk_size + 1 : nelts], sprintf_partial_line)
  end
end


function fortran_write_array{T <: AbstractFloat}(io :: IOStream, data :: Vector{T}, chunk_size :: Int, width :: Int)
  nelts = length(data)
  nlines, leftover = fldmod(nelts, chunk_size)
  sprintf_full_line = eval(Symbol("sprintf_real_$(chunk_size)_$width"))
  for k = 0:nlines-1
    fortran_write_line(io, data[k * chunk_size + 1 : (k + 1) * chunk_size], sprintf_full_line)
  end
  if leftover > 0
    sprintf_partial_line = eval(Symbol("sprintf_real_$(leftover)_$width"))
    fortran_write_line(io, data[nlines * chunk_size + 1 : nelts], sprintf_partial_line)
  end
end


global const data_types = ["ord", "rhs", "sln", "est", "evl", "svl", "evc",
                           "svc", "sbv", "sbm", "sbp", "ipt", "icv", "lvl",
                           "geo", "avl"]

global const organizations = ['s', 'd', 'e']
global const positions = ['r', 'l', 's']

"""Write auxilliary data to file in Rutherford-Boeing format.

    write_rb_aux(filename, vecs)

# Arguments
* `filename::String`: the name of the output file
* `vecs::SparseMatrixCSC`: the data to write to file in sparse format

# Keyword arguments
* `title::String`: a title to be recorded in the output file (72 chars max)
* `key::String`: a shorthand title (8 chars max)
* `caseid::String`: case identifier to distinguish multiple sets of data (8 chars max)
* `dattyp::String`: auxilliary data type (rhs, sln, etc.)
* `positn::Char`: `l` (left), `r` (right) or `s` (symmetric) depending on the type of data
* `orgniz::Char`: `d` (dense), `s` (sparse) or `e` (elemental)
* `nauxd::Int`: supplementary integer parameters
* `precision::Int`: number of significant digits to write for real data
"""
function write_rb_aux{Ti <: Integer, Tr <: Real}(
  filename :: String, vecs :: SparseMatrixCSC{Tr,Ti};
  title :: String="Generic", key :: String="Generic", caseid :: String="Generic",
  dattyp :: String="rhs", positn :: Char='r',
  orgniz :: Char='d', nauxd :: Int=0, precision :: Int=16)

  dattyp in data_types || error("unknown data type: $dattyp")
  positn in positions || error("unknown position: $positn")
  if dattyp == "rhs"
    orgniz in organizations || error("unknown organization: $orgniz")
  end

  nrow, nvecs = size(vecs)
  if dattyp in ["evl", "svl", "lvl", "sbp"]
    nvecs == 1 || error("data type $dattyp requires a single vector")
  end
  if dattyp in ["evl", "svl", "lvl", "sbv", "sbm", "sbp", "avl"]
    positn = ' '
  end
  if dattyp != "rhs"
    orgniz = ' '
  end

  numerf = dattyp == "ord" ? 'i' : (dattyp in ["ipt", "icv"] ? 'p' : 'r')
  if orgniz != 'e'
    nauxd = 0
    if orgniz == 's' || dattyp in ["icv", "ipt"]
      nauxd = nvecs
    elseif orgniz == 'd'
      nauxd = prod(size(vecs))
    end
  end

  fm1 = fm2 = fm3 = " "

  if dattyp in ["ipt", "icv"]
    ndigits1, fm1, nitems1 = get_int_fmt(nauxd + 1)
    ndigits2, fm2, nitems2 = get_int_fmt(nrow)
  elseif dattyp == "ord"
    ndigits1, fm1, nitems1 = get_int_fmt(nrow)
  else
    precision = min(max(1, precision), 16)
    if dattyp == "rhs" && orgniz == 's'
      ndigits1, fm1, nitems1 = get_int_fmt(nauxd + 1)
      ndigits2, fm2, nitems2 = get_int_fmt(nrow)
      ndigits3, fm3, nitems3 = get_real_fmt(precision)
    else
      ndigits1, fm1, nitems1 = get_real_fmt(precision)
    end
  end

  fp = open(filename * ".rb", "w")

  # write header
  lt = min(length(title), 72)
  lk = min(length(key), 8)
  lc = min(length(caseid), 8)
  print(fp, @sprintf("%-72s%-8s\n", title[1:lt], key[1:lk]))
  print(fp, @sprintf("%-3s%c%c %-8s %c %-13d %-13d %-13d\n",
                     dattyp, positn, orgniz, caseid[1:lc], numerf, nrow, nvecs, nauxd))
  print(fp, @sprintf("%-20s%-20s%-20s\n", fm1, fm2, fm3))

  # write integer data
  if (dattyp == "rhs" && orgniz == 's') || dattyp in ["ipt", "icv"]
    fortran_write_array(fp, vecs.colptr, nitems1, ndigits1)
    fortran_write_array(fp, vecs.rowval, nitems2, ndigits2)
  elseif dattyp == "ord"
    fortran_write_array(fp, vecs.rowval, nitems1, ndigits1)
  end

  # write entries
  if dattyp == "rhs" && orgniz == 's'
    fortran_write_array(fp, vecs.nzval, nitems3, ndigits3)
  elseif !(dattyp in ["ipt", "icv"])
    fortran_write_array(fp, vecs.nzval, nitems1, ndigits1)
  end
  close(fp)
end


"""Write sparse matrix to file in Rutherford-Boeing format.

    write_rb(filename, matrix)

# Arguments
* `filename::String`: the name of the output file
* `matrix::SparseMatrixCSC`: the matrix to write to file in sparse format

# Keyword arguments
* `title::String`: a title to be recorded in the output file (72 chars max)
* `key::String`: a shorthand title (8 chars max)
* `symmetric::Bool`: only store the lower triangle of the matrix
* `skew::Bool`: only store the strict lower triangle of the matrix
* `pattern_only::Bool`: only store the sparsity pattern
* `precision::Int`: number of significant digits to write for real data
"""
function write_rb{Ti <: Integer, Tr <: Real}(
  filename :: String, matrix :: SparseMatrixCSC{Tr,Ti};
  title :: String="Generic", key :: String="Generic",
  symmetric :: Bool=false, skew :: Bool=false,
  pattern_only :: Bool=false, precision :: Int=16)

  nrow, ncol = size(matrix)
  ne = nnz(matrix)

  mxtype0 = pattern_only ? "p" : "r"  # only real matrices for now
  mxtype1 = nrow == ncol ? (symmetric ? "s" : "u") : "r"
  mxtype2 = "a"  # only assembled matrices for now
  mxtype = mxtype0 * mxtype1 * "a"

  ndigits_ptr, fmt_ptr, nitems_ptr = get_int_fmt(ne + 1)
  nlines_ptr = Int(floor(ncol / nitems_ptr)) + 1

  ndigits_ind, fmt_ind, nitems_ind = get_int_fmt(nrow)
  nlines_ind = Int(floor((ne - 1) / nitems_ind)) + 1

  if pattern_only
    nlines_val = 0
    fmt_val = " "
  else
    ndigits_val, fmt_val, nitems_val = get_real_fmt(precision)
    nlines_val = Int(floor((ne - 1) / nitems_val)) + 1
  end

  nlines = nlines_ptr + nlines_ind + nlines_val
  neltvl = 0

  fp = open(filename * ".rb", "w")

  lt = min(length(title), 72)
  lk = min(length(key), 8)
  print(fp, @sprintf("%-72s%-8s\n", title[1:lt], key[1:lk]))
  print(fp, @sprintf("%14d %13d %13d %13d\n", nlines, nlines_ptr, nlines_ind, nlines_val))
  print(fp, @sprintf("%3s            %13d %13d %13d %13d\n", mxtype, nrow, ncol, ne, neltvl))
  print(fp, @sprintf("%-16s%-16s%-20s\n", fmt_ptr, fmt_ind, fmt_val))

  # ensure row indices are sorted in each column
  sortsparse!(matrix.colptr, matrix.rowval, matrix.nzval)

  fortran_write_array(fp, matrix.colptr, nitems_ptr, ndigits_ptr)
  fortran_write_array(fp, matrix.rowval, nitems_ind, ndigits_ind)
  pattern_only || fortran_write_array(fp, matrix.nzval, nitems_val, ndigits_val)

  close(fp)
end
