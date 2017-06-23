import Eiscor

# function to test rfuqr routine
function testrfuqr(T,n)

  # allocate storage 
  u = zeros(Complex{T},n)
  u[n] += (-1)^(n-1)
  vv = ones(T,n)
  vv[n] *= 0
  its = zeros(Int,n-1)
  
  # call routine
  Eiscor.RootFreeUnitary.rfuqr!(u,vv,its)

  # print output
  print("\n")
  for ii = 1:n
    @printf "%+1.15e %+1.15e %1.15e\n" real(u[ii]) imag(u[ii]) abs(abs(u[ii])-1)
  end  

end 

# get array of supported types
STs = Eiscor.SupportedTypes

# loop through supported types
for T = STs
  
  # print current type
  print(" current type: ",T,"\n")

  # get current tolerance
  tol = eps(T)
  @printf " current tol: %1.4e\n" tol
    
#  # initialize test variables
#  w = convert(Complex{T},complex(randn(),randn())); w /= abs(w)
#  cc = abs(convert(T,randn())); ss = abs(convert(T,randn()))
#  nn = ss + cc; cc /= nn; ss /= nn
#  u = convert(Complex{T},complex(randn(),randn())); 
#  vv = abs(convert(T,randn()))
#  nn = abs(u)^2 + vv; u /= sqrt(nn); vv /= nn
#  rho = convert(Complex{T},complex(randn(),randn())); rho /= abs(rho)
#
#  # call root free turnover
#  Eiscor.RootFreeUnitary.turnover!(w,cc,ss,u,vv,rho)

  # call test routine
  testrfuqr(T,2^4)

#  # compute errors
#  errW = abs(abs(w)-1)
#  errCCSS = abs((cc+ss)-1)
#  errUVV = abs((abs(u)^2+vv)-1)
#  @printf "    errW: %1.4e\n" errW
#  @printf " errCCSS: %1.4e\n" errCCSS
#  @printf "  errUVV: %1.4e\n" errUVV

  # print end of test
  print(" end of ",T," tests\n\n")

end

