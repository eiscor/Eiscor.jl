import Eiscor

# get array of supported types
STs = Eiscor.SupportedTypes

# loop through supported types
for T = STs
  
  # print current type
  print(" current type: ",T,"\n")

  # get current tolerance
  tol = eps(T)
  @printf " current tol: %1.4e\n" tol
    
  # initialize test variables
  w = convert(Complex{T},complex(randn(),randn())); w /= abs(w)
  cc = abs(convert(T,randn())); ss = abs(convert(T,randn()))
  nn = ss + cc; cc /= nn; ss /= nn
  u = convert(Complex{T},complex(randn(),randn())); 
  vv = abs(convert(T,randn()))
  nn = abs(u)^2 + vv; u /= sqrt(nn); vv /= nn
  rho = convert(Complex{T},complex(randn(),randn())); rho /= abs(rho)

  # call root free turnover
  Eiscor.RootFreeUnitary.turnover!(w,cc,ss,u,vv,rho)

  # compute errors
  errW = abs(abs(w)-1)
  errCCSS = abs((cc+ss)-1)
  errUVV = abs((abs(u)^2+vv)-1)
  @printf "    errW: %1.4e\n" errW
  @printf " errCCSS: %1.4e\n" errCCSS
  @printf "  errUVV: %1.4e\n" errUVV

  # print end of test
  print(" end of ",T," tests\n\n")

end

## allocate storage for it counts
#its = zeros(Int64,l-1)
#
## allocate linear storage for compressed matrix H
#a = zeros(Complex{T},l)
#a[l] = (-1)^(l+1) # final entry
#bb = ones(T,l)
#bb[l] = 0
##a = randn(l) + 0*im
##Va[l] = 1 # final entry
##bb = randn(l).^2
##bb[l] = 0
##x = abs(a).^2 + bb
##bb ./= x
#V#a ./= sqrt(x)
#
#rho = randn() + randn()*im
#rho /= abs(rho)
