# Helper functions to read Harwell-Boeing and Rutherford-Boeing data.

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


function read_array(io :: IO, n :: Int, fmt :: String; is_complex :: Bool=false)
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
    x /= 10.0^scale
  end
  return is_complex ? [Complex64(x[i], x[i+1]) for i = 1 : 2 : n-1] : x
end
