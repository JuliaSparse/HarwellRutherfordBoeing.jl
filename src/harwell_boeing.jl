mutable struct HBMeta
  # Metadata attached to a Harwell-Boeing matrix.
  title :: AbstractString
  key   :: AbstractString

  totcrd :: Int
  ptrcrd :: Int
  indcrd :: Int
  valcrd :: Int
  rhscrd :: Int

  mxtype :: AbstractString
  nrow :: Int
  ncol :: Int
  nnzero :: Int
  neltvl :: Int

  hermitian :: Bool
  assembled :: Bool

  ptrfmt :: AbstractString
  indfmt :: AbstractString
  valfmt :: AbstractString
  rhsfmt :: AbstractString

  nrhs   :: Int
  rhstyp :: AbstractString
  nrhsix :: Int
end


const RHSType = Union{Array{Float64,2}, Array{ComplexF32,2}, SparseMatrixCSC}


mutable struct HarwellBoeingMatrix

  meta   :: HBMeta
  matrix :: SparseMatrixCSC
  rhs    :: RHSType  # Right-hand side, if any.
  guess  :: RHSType  # Initial guesses, if any.
  sol    :: RHSType  # Solutions, if any.

  function HarwellBoeingMatrix(file_name :: AbstractString)
    hb = open(file_name)

    # Read header.
    seekstart(hb)
    line  = readline(hb)
    title = strip(line[1:72])      # A72
    key   = strip(line[73:end-1])  # A8
    totcrd, ptrcrd, indcrd, valcrd, rhscrd = map(s -> parse(Int, s),
                                                 split(chomp(readline(hb))))

    line = readline(hb)
    mxtype = line[1:3]   # A3
    # Skip 11 blanks.
    nrow, ncol, nnzero, neltvl = map(s -> parse(Int, s),
                                     [line[14 + 14*(i-1)+1 : 14 + 14*i] for i = 1 : 4])

    pattern_only = (mxtype[1] == 'P')
    is_complex = (mxtype[1] == 'C')
    hermitian = (mxtype[2] == 'S')
    assembled = (mxtype[3] == 'A')
    data_type = is_complex ? ComplexF32 : Float64

    line = readline(hb)
    ptrfmt = strip(line[1:16])  # 2A16
    indfmt = strip(line[17:32])

    valfmt = pattern_only ? "" : strip(line[33:52])  # 2A20
    rhsfmt = (pattern_only | rhscrd == 0) ? "" : strip(line[53:72])

    # Read right-hand side info, if present.
    if rhscrd > 0
      line = readline(hb)
      rhstyp = line[1:3]  # A3
      nrhs, nrhsix = map(s -> parse(Int, s), split(chomp(line[15:end])))

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
    # Ensure row indices are sorted in each column.
    sortsparse!(ip, ind, vals)
    matrix = SparseMatrixCSC(nrow, ncol, ip, ind, vals)

    # Read right-hand sides, if any.
    if (nrhs > 0) & (rhstyp[1] != 'F') & assembled  # Sparse rhs.
        rhsptr  = read_array(hb, nrhs+1, ptrfmt)
        rhsind  = read_array(hb, nrhsix, indfmt)
        rhsvals = read_array(hb, nrhsix, rhsfmt, is_complex=is_complex)
        # Ensure row indices are sorted in each column.
        sortsparse!(rhsptr, rhsind, rhsvals)
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
