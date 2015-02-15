include("hrb_utils.jl")

type RBMeta
  # Metadata attached to Rutherford-Boeing data.
  title :: String
  key   :: String

  mxtype :: String
  nrow :: Int
  ncol :: Int
  nnzero :: Int
  neltvl :: Int

  hermitian :: Bool
  assembled :: Bool
  pattern_only :: Bool

  ptrfmt :: String
  indfmt :: String
  valfmt :: String

  dattyp :: String
  positn :: Char
  orgniz :: Char
  caseid :: String
  numerf :: Char
  auxfm1 :: String
  auxfm2 :: String
  auxfm3 :: String
end


RBDataType = Union(Array{Int64,2}, Array{Float64,2}, Array{Complex128,2}, SparseMatrixCSC)


type RutherfordBoeingData

  meta :: RBMeta
  data :: RBDataType

  function RutherfordBoeingData(file_name :: String)
    data_types = ['r' => Float64, 'c' => Complex128, 'i' => Int64]
    rb = open(file_name)

    # Read header.
    seekstart(rb)
    line  = readline(rb)
    title = strip(line[1:72])      # A72
    key   = strip(line[73:end-1])  # A8
    buffer1 = readline(rb)         # A80
    buffer2 = readline(rb)         # A80

    if buffer2[3] in ['a', 'e']
      # Read a matrix.
      mxtype = buffer2[1:3]
      nrow, ncol, nnzero, neltvl = map(int, split(chomp(buffer2[4:end])))

      pattern_only = (mxtype[1] in ['p', 'q'])
      if pattern_only
        ptrfmt, indfmt = split(strip(readline(rb)))
      else
        data_type = data_types[mxtype[1]]
        ptrfmt, indfmt, valfmt = split(strip(readline(rb)))
      end

      # Read integer data.
      np1 = (mxtype[2:3] == "re") ? (2*ncol + 1) : (ncol + 1)
      ip  = read_array(rb, np1, ptrfmt)
      ind = read_array(rb, nnzero, indfmt)

      nreal = neltvl > 0 ? neltvl : nnzero

      if pattern_only
        vals = ones(Float64, nreal)  # To be able to build a SparseMatrixCSC
      else
        vals = read_array(rb, nreal, valfmt)
      end
      data = SparseMatrixCSC(nrow, ncol, ip, ind, vals)

      meta = RBMeta(title, key, mxtype, nrow, ncol, nnzero, neltvl,
                    mxtype[2] == 's', mxtype[3] == 'a', pattern_only,
                    ptrfmt, indfmt, pattern_only ? "" : valfmt,
                    "", '\0', '\0', "", '\0', "", "", "")

    else

      # Read supplementary data.
      types_with_values = ["rhs", "sln", "est", "ipt", "icv"]

      dattyp = buffer1[1:3]
      positn = buffer1[4]
      orgniz = buffer1[5]
      caseid = buffer1[7:14]
      numerf = buffer1[16]
      intvals = map(int, split(chomp(buffer1[17:end])))
      auxfmts = split(strip(buffer2))
      auxfm1 = auxfm2 = auxfm3 = ""
      if dattyp in types_with_values
        (nrow, ncol, nnzero) = intvals
      else
        (nrow, ncol) = intvals
      end
      if dattyp == "ord"
        nnzero = nrow * ncol
      end

      if (dattyp == "rhs" && orgniz == 's') || (dattyp in ["ord", "ipt", "icv"])
        # Read indices.
        auxfm = auxfmts[1]
        if (dattyp == "rhs" && orgniz == 's')
          auxfm1 = auxfmts[1]
          ip = read_array(rb, ncol+1, auxfm1)
          auxfm = auxfmts[2]
        end
        ind = read_array(rb, nnzero, auxfm)
      end

      if dattyp in types_with_values
        # Read values.
        if dattyp != "rhs"
          nnzero = nrow * ncol
        end
        if dattyp == "rhs" && orgniz == 's'
          auxfm = auxfmts[3]
        else
          auxfm = auxfmts[1]
        end

        val = read_array(rb, nnzero, auxfm)

        if (dattyp == "rhs" && orgniz == 's') || (dattyp in ["ipt", "icv"])
          data = SparseMatrixCSC(nrow, ncol, ip, ind, val)
        else
          data = reshape(val, (nrow, ncol))
        end
      else
        data = reshape(ind, (nrow, ncol))
      end

      meta = RBMeta(title, key, "", nrow, ncol, nnzero, 0,
                    false, orgniz != 'e', false,
                    "", "", "",
                    dattyp, positn, orgniz, caseid, numerf,
                    auxfm1, auxfm2, auxfm3)
    end

    close(rb)

    new(meta, data)
  end
end


# Displaying Rutherford-Boeing data instances.

import Base.show, Base.print
function show(io :: IO, rb :: RutherfordBoeingData)
  typ = (rb.meta.mxtype == "") ? rb.meta.dattyp : rb.meta.mxtype
  s  = @sprintf("Rutherford-Boeing data %s of type %s\n", rb.meta.key, typ)
  s *= @sprintf("%d rows, %d cols, %d nonzeros\n", rb.meta.nrow, rb.meta.ncol, rb.meta.nnzero)
  print(io, s)
end

function print(io :: IO, rb :: RutherfordBoeingData)
  @printf("Rutherford-Boeing data %s\n", rb.meta.key)
  @printf("%s\n", rb.meta.title)
  @printf("%d rows, %d cols, %d nonzeros\n", rb.meta.nrow, rb.meta.ncol, rb.meta.nnzero)
  if rb.meta.mxtype != ""
    if rb.meta.mxtype[1] == 'P'
      dtype = "pattern only"
    elseif rb.meta.mxtype[1] == 'R'
      dtype = "real"
    else
      dtype = "complex"
    end
    herm = rb.meta.hermitian ? "hermitian" : "non-hermitian"
    assm = rb.meta.assembled ? "assembled" : "elemental"
    @printf("(%s, %s, %s)\n", dtype, herm, assm)
  end
  display(rb.data)
end
