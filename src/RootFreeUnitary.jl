module RootFreeUnitary

# unitary root free qr function
function rfuqr!(u,vv,its)

  # fetch dimension
  n = length(u)
  
#  # fetch type
#  T = typeof(h[1])
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
  # allocate memory
  str = one(Int)
  stp = n-1
  zer = zero(Int)
  itc = zero(Int)
  itm = 20*n

  # initialize variables
  str = 1
  stp = n-1

  # main loop
  for kk = 1:itm

    # exit if finished
    if stp <= 0
      for ii = 1:n-1
        u[n+1-ii] *= u[n-ii]'
      end
      break
    end

    # check for deflation
    zer = deflate!(slice(u,str:stp),slice(vv,str:stp))

    # store itc if deflation
    if zer > 0
      its[str+zer-1] += itc
      itc = 0
    end
    
    # if 1x1 block remove and check again 
    if stp == (str+zer-1) 

      # update indices
      stp -= 1
      zer = 0
      str = 1
    
    # if greater than 1x1 chase a bulge
    else

      # check zer
      if (zer > 0) 
        str += zer
      end

      # set nu for top deflations
      if str > 1 
        nu = u[str-1]'
      else
        nu = 0*u[1] + 1
      end

      # perform singleshift iteration
      singlestep!(slice(u,str:stp+1),slice(vv,str:stp+1),nu)
     
      # update indices
      itc = itc + 1
 
    end
    
    # if ITMAX hit
    if (kk == itm) 
      its[str+stp-1] = itc
    end

  # end main loop
  end

# end rfuqr
end

# deflation check
function deflate!(u,vv)

  # get length
  n = length(u)

  # set tolerance
  tol = eps(one(typeof(vv[1])))^2

  # intialize zer
  zer = 0
  
  # check for deflation
  for ii=1:n
  
    # deflate if subdiagonal is small enough
    if vv[n+1-ii] < tol 
        
      # set zer
      zer = n+1-ii

      # set rotation to diagonal
      vv[zer] *= 0
        
      # renormalize U
      xx = real(u[zer])^2 + imag(u[zer])^2
      u[zer] *= .5*(3-xx)

      # exit loop
      break

    end

  end
 
  # return
  zer

# end of deflation check
end

# one singleshift step
function singlestep!(u,vv,nu)

  # get length
  n = length(u)

  # bottom 2x2 block
  A = zeros(typeof(u[1]),2,2)
  A[1,1] = u[n-1]
  A[2,2] = A[1,1]'
  A[2,1] = sqrt(vv[n-1])
  A[1,2] = -A[2,1]
  A[:,2] *= u[n]
  if n > 2
    if abs(u[n-2]) > 0
      A[1,:] *= u[n-2]'/abs(u[n-2])
    end
  end

  # Wilkinson shift
  evs = zeros(typeof(u[1]),2)
  evs[1] = A[1,1] + A[2,2]
  evs[2] = sqrt(evs[1]^2 - 4*(A[1,1]*A[2,2]-A[2,1]*A[1,2]))
  evs[1] = .5*(evs[1]+evs[2])
  evs[2] = A[1,1] + A[2,2] - evs[1]
  if abs(A[2,2]-evs[1]) < abs(A[2,2]-evs[2])
    rho = evs[1]/abs(evs[1])
  else 
    rho = evs[2]/abs(evs[2])
  end

  # initialize 
  w = -rho
  cc = 0*vv[1] + 1
  ss = 0*cc

  # loop
  for ii=0:n-1

    # set ut and vvt
    ut = nu*u[ii+1]
    vvt = 1*vv[ii+1]

    # turnover
    w,cc,ss,ut,vvt = turnover(w,cc,ss,ut,vvt,rho)

    # store ut and vvt
    if ii > 0 
      u[ii] = nu'*ut
      vv[ii] = vvt
    end

  end

# end of singlestep function
end

# root free turnover function
function turnover(w,cc,ss,u,vv,rho)

  # store inputs
  uold = 1*u
  xx = 1*ss

  # z and zz
  z = w + u
  zz = real(z)^2 + imag(z)^2

  # u
  u = -rho'*(cc*z - u)

  # cc and ss
  cc *= zz
  if cc == 0
    nn = vv
    ss = 0*ss + 1
    w = 0*w + 1
  else
    nn = vv + cc
    cc /= nn
    ss = vv/nn
    w = -rho*(z*(vv/zz) + uold)
  end

  # vv
  vv = xx*nn

  # normalize output
  xx = real(u)^2 + imag(u)^2 + vv
  u *= .5*(3.-xx)
  vv /= xx
  xx = real(w)^2 + imag(w)^2
  w *= .5*(3.-xx)

  # return
  w,cc,ss,u,vv

# end of root free turnover function
end

# end RootFreeUnitary module
end
