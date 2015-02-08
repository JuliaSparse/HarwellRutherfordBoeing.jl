module HarwellBoeing

export HarwellBoeingMatrix


type HBMeta
  # Metadata attached to a Harwell-Boeing matrix.
  title :: String
  key   :: String

  totcrd :: Int
  ptrcrd :: Int
  indcrd :: Int
  valcrd :: Int
  rhscrd :: Int

  mxtype :: String
  nrow :: Int
  ncol :: Int
  nnzero :: Int
  neltvl :: Int

  hermitian :: Bool
  assembled :: Bool

  ptrfmt :: String
  indfmt :: String
  valfmt :: String
  rhsfmt :: String

  nrhs   :: Int
  rhstyp :: String
  nrhsix :: Int
end


RHSType = Union(Array{Float64,2}, Array{Complex64,2}, SparseMatrixCSC)


type HarwellBoeingMatrix

  meta   :: HBMeta
  matrix :: SparseMatrixCSC
  rhs    :: RHSType  # Right-hand side, if any.
  guess  :: RHSType  # Initial guesses, if any.
  sol    :: RHSType  # Solutions, if any.

  function HarwellBoeingMatrix(file_name :: String)
    hb = open(file_name)

    # Read header.
    seekstart(hb)
    line  = readline(hb)
    title = strip(line[1:72])      # A72
    key   = strip(line[72:end-1])  # A8
    totcrd, ptrcrd, indcrd, valcrd, rhscrd = map(int, split(chomp(readline(hb))))

    line = readline(hb)
    mxtype = line[1:3]   # A3
    # Skip 11 blanks.
    nrow, ncol, nnzero, neltvl = map(int, [line[14 + 14*(i-1)+1 : 14 + 14*i] for i = 1 : 4])

    pattern_only = (mxtype[1] == 'P')
    is_complex = (mxtype[1] == 'C')
    hermitian = (mxtype[2] == 'S')
    assembled = (mxtype[3] == 'A')
    data_type = is_complex ? Complex64 : Float64

    line = readline(hb)
    ptrfmt = strip(line[1:16])  # 2A16
    indfmt = strip(line[17:32])

    valfmt = pattern_only ? "" : strip(line[33:52])  # 2A20
    rhsfmt = (pattern_only | rhscrd == 0) ? "" : strip(line[53:72])

    # Read right-hand side info, if present.
    if rhscrd > 0
      line = readline(hb)
      rhstyp = line[1:3]  # A3
      nrhs, nrhsix = map(int, split(chomp(line[15:end])))

    else
      rhstyp = "   "
      nrhs = 0
      nrhsix = 0
    end

    # Read integer data.
    ip = read_array(hb, ncol+1, ptrfmt)
    ind = read_array(hb, nnzero, indfmt)

    # Read values, if supplied.
    if pattern_only
      vals = ones(Float64, nnzero)  # Just so we can build a SparseMatrixCSC.
    else
      vallen = assembled ? nnzero : neltvl
      vals = read_array(hb, vallen, valfmt, is_complex=is_complex)
    end
    matrix = SparseMatrixCSC(nrow, ncol, ip, ind, vals)

    # Read right-hand sides, if any.
    if (nrhs > 0) & (rhstyp[1] != 'F') & assembled  # Sparse rhs.
        rhsptr  = read_array(hb, nrhs+1, ptrfmt)
        rhsind  = read_array(hb, nrhsix, indfmt)
        rhsvals = read_array(hb, nrhsix, rhsfmt, is_complex=is_complex)
        rhs = SparseMatrixCSC(nrow, nrhs, rhsptr, rhsind, rhsvals)

    else
      lenrhs = rhstyp[1] == 'F' ? nrow : nnzero
      rhs    = zeros(data_type, (nrow, nrhs))

      for i = 1 : nrhs
        rhs[:,i] = read_array(hb, lenrhs, rhsfmt, is_complex=is_complex)
      end
    end

    guess  = zeros(data_type, (nrow, rhstyp[2] == 'G' ? nrhs : 0))
    sol    = zeros(data_type, (nrow, rhstyp[3] == 'X' ? nrhs : 0))

    # Read initial guesses, if any.
    if rhstyp[2] == 'G'
      for i = 1 : nrhs
        guess[:,i] = read_array(hb, nrow, rhsfmt, is_complex=is_complex)
      end
    end

    # Read solution vectors, if any.
    if rhstyp[3] == 'X'
      for i = 1 : nrhs
        sol[:,i] = read_array(hb, nrow, rhsfmt, is_complex=is_complex)
      end
    end

    close(hb)

    meta = HBMeta(title, key, totcrd, ptrcrd, indcrd, valcrd, rhscrd,
                  mxtype, nrow, ncol, nnzero, neltvl, hermitian, assembled,
                  ptrfmt, indfmt, valfmt, rhsfmt, nrhs, rhstyp, nrhsix)

    new(meta, matrix, rhs, guess, sol)
  end
end


# Helper functions.

function decode_int_fmt(fmt :: String)
  if fmt[1] == '('
    fmt = uppercase(fmt[2:end-1])
  end
  return map(int, split(fmt, 'I'))
end


function decode_real_fmt(fmt :: String)
  fmt = join(split(fmt))  # Remove all white spaces.
  if fmt[1] == '('
    fmt = uppercase(fmt[2:end-1])
  end
  scale = 0
  if (',' in fmt)  # Process scale factor, e.g., 1P,5D16.9
    scale, fmt = split(fmt, ',')
    scale, _ = split(scale, 'P')
  elseif ('P' in fmt)
    scale, fmt = split(fmt, 'P')
  end
  scale = int(scale)

  fmt1 = split(fmt, '.')[1]
  if search(fmt1, 'E') > 0
    (npl, len) = map(int, split(fmt1, 'E'))
  elseif search(fmt1, 'D') > 0
    (npl, len) = map(int, split(fmt1, 'D'))
  elseif search(fmt1, 'F') > 0
    (npl, len) = map(int, split(fmt1, 'F'))
  else
    error("Malformed real format")
  end

  return (npl, len, scale)
end


function standardize_real(number_as_str :: String)
  s = join(split(number_as_str))  # for numbers in the form "0.24555165E 00".
  # change ".16000000+006" to ".16000000e+006". The first char could be +/-.
  if (search(s, 'E') + search(s, 'D') + search(s, 'e') + search(s, 'd')) == 0
    if search(s[2:end], '+') > 0
      s = s[1:1] * join(split(s[2:end], '+'), "e+")
    elseif search(s[2:end], '-') > 0
      s = s[1:1] * join(split(s[2:end], '-'), "e-")
    end
  end
  return s
end


function read_array(io :: IO, n :: Int, fmt :: String; is_complex :: Bool = false)
  if 'I' in fmt
    scale = 0
    (npl, len) = decode_int_fmt(fmt)
    conv = int
    typ = Int64
  else
    (npl, len, scale) = decode_real_fmt(fmt)
    conv = float
    typ = Float64
    if is_complex
      n *= 2
    end
  end

  x = zeros(typ, n)
  for j = 1 : div(n, npl)
    if typ == Float64
      line = join(split(uppercase(readline(io)), 'D'), 'e')
    else
      line = readline(io)
    end
    chunk = [line[len*(i-1)+1:len*i] for i = 1 : npl]
    if typ == Float64
      chunk = map(standardize_real, chunk)
    end
    x[npl * (j-1) + 1 : npl * j] = map(conv, chunk)
  end
  rem = mod(n, npl)
  if rem > 0
    if typ == Float64
      line = join(split(uppercase(readline(io)), 'D'), 'e')
    else
      line = readline(io)
    end
    chunk = [line[len*(i-1)+1:len*i] for i = 1 : rem]
    if typ == Float64
      chunk = map(standardize_real, chunk)
    end
    x[end-rem+1 : end] = map(conv, chunk)
  end
  if scale != 0
    x /= 10^scale
  end
  return is_complex ? [Complex64(x[i], x[i+1]) for i = 1 : 2 : n-1] : x
end


# Displaying Harwell-Boeing matrix instances.

import Base.show, Base.print
function show(io :: IO, hb :: HarwellBoeingMatrix)
  s  = @sprintf("Harwell-Boeing matrix %s of type %s\n", hb.meta.key, hb.meta.mxtype)
  s *= @sprintf("%d rows, %d cols, %d nonzeros\n", hb.meta.nrow, hb.meta.ncol, hb.meta.nnzero)
  s *= @sprintf("%d right-hand sides, %d guesses, %d solutions\n",
                size(hb.rhs,2), size(hb.guess,2), size(hb.sol,2))
  print(io, s)
end

function print(io :: IO, hb :: HarwellBoeingMatrix)
  @printf("Harwell-Boeing matrix %s\n", hb.meta.key)
  @printf("%s\n", hb.meta.title)
  @printf("%d rows, %d cols, %d nonzeros\n", hb.meta.nrow, hb.meta.ncol, hb.meta.nnzero)
  if hb.meta.mxtype[1] == 'P'
    dtype = "pattern only"
  elseif hb.meta.mxtype[1] == 'R'
    dtype = "real"
  else
    dtype = "complex"
  end
  herm = hb.meta.hermitian ? "hermitian" : "non-hermitian"
  assm = hb.meta.assembled ? "assembled" : "elemental"
  @printf("(%s, %s, %s)\n", dtype, herm, assm)
  @printf("%d right-hand sides, %d guesses, %d solutions\n",
           size(hb.rhs,2), size(hb.guess,2), size(hb.sol,2))
  if size(hb.rhs,2) > 0
    @printf("Right-hand side(s):\n");
    display(hb.rhs'); @printf("\n")
  end
  if size(hb.guess,2) > 0
    @printf("Initial guess(es):\n");
    display(hb.guess'); @printf("\n")
  end
  if size(hb.sol,2) > 0
    @printf("Solution(s):\n");
    display(hb.sol'); @printf("\n")
  end
end

end
