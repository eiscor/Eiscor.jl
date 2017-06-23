module RootFreeUnitary

## root free qr function
#function rfuqr(h,its)
#
#  # fetch type
#  T = typeof(h[1])
#
#  # fetch dimension
#  l = round(Int,(length(h)+3)/5)
#  
#  # check length
#  @assert 5*l-3 == length(h)  
#  @assert l-1 == length(its)
#
#  # check unitarity
#  tol = eps(T)
#  @assert norm(h[1:5:end].^2+h[2:5:nd].^2-1,Inf) < 10*tol
#  @assert norm(h[3:5:end].^2+h[4:5:nd].^2+h[5:5:end].^2-1,Inf) < 10*tol
#
#  # square subdiagonals of unitary matrix
#  h[5:5:end] .*= h[5:5:end]
#
#  # allocate memory
#  str = zero(Int)
#  stp = zero(Int)
#  zer = zero(Int)
#  itc = zero(Int)
#
#  # initialize variables
#  str = 1
#  stp = l-1
#
#  # main loop
#  for ii = 1:20*l
#
#    # exit if finished
#    if stp == 0
#      break
#    end
#
#    # check for deflation
#    zer = deflate(h)
#
#    # zero at bottom
#    if zer == stp
#      its[stp] = itcnt
#      itcnt = 0
#      stp -= 1
#      zer = 0
#    # zero anywhere else
#    elseif zer > 0
#      its[zer] = itcnt
#      itcnt = 0
#      str = zer + 1
#
#      # compute shift
#
#      # perform singlestep
#
#    end
#
#  # end main loop
#  end
#
#end
#
## one singleshift step
#function singlestep(a,bb,rho)
#
#  # get length
#  n = length(a)
#
#  # get type
#  T = typeof(bb[1])
#  z = zero(Complex{T})
#  w = one(Complex{T})
#  zz = zero(T)
#  nn = zero(T)
#  cc = one(T)
#  ss = zero(T)
#
#  # initialize 
#  w = -rho
#  z = a[1] + w
#  zz = real(z)^2 + imag(z)^2
#  nn = bb[1] + zz
#
#  # loop
#  for j=1:n-1
#    cc = cc*zz/nn
#    ss = bb[j]/nn
#    w = -rho*w'*(z*z)/zz
#    z = a[j+1] + w
#    zz = real(z)^2 + imag(z)^2
#    nn = bb[j+1] + cc*zz
#    a[j] = rho'*(a[j+1] - cc*z)
#    bb[j] = ss*nn
#  end
#
## end of singlestep function
#end

# root free turnover function
function turnover!(w,cc,ss,u,vv,rho)

  # store inputs
  uold = u
  xx = ss

  # z and zz
  z = w + u
  zz = real(z)^2 + imag(z)^2

  # u
  u = -rho'*(cc*z - u)

  # cc and ss
  cc = cc*zz
  if cc == 0
    nn = vv
    ss = 1.
    w = 1.
  else
    nn = vv + cc
    cc = cc/nn
    ss = vv/nn
    w = -rho*(z*(vv/zz) + u)
  end

  # vv
  vv = xx*nn

  # normalize output
  xx = real(u)^2 + imag(u)^2 + vv
  u = .5*(3.-xx)*u
  vv = vv/xx
  xx = real(w)^2 + imag(w)^2
  w = .5*(3.-xx)*w

# end of root free turnover function
end

# end RootFreeUnitary module
end
