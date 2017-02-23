reload("Eiscor")

n = 2^2

A = randn(n,n)
A = A + A'
for ii=1:n-2
  for jj=1:ii
    A[jj,ii+2] = 0
    A[ii+2,jj] = 0
  end
end
na = norm(A)

# initialize Ttot
Ttot = eye(n)

m = n
for kk=1:10*n

# first rfr
C,S,NRM = Eiscor.RFR.RFR2Vec2Gen(A[1,1]-A[m,m],A[2,1])
T = [C -S;S C]
Tinv = T'/NRM
A[:,1:2] = A[:,1:2]*T
A[1:2,:] = Tinv*A[1:2,:]
Ttot[:,1:2] = Ttot[:,1:2]*T

# chase bulge
for ii=1:m-2

C,S,NRM = Eiscor.RFR.RFR2Vec2Gen(A[ii+1,ii],A[ii+2,ii])
T = [C -S;S C]
Tinv = T'/NRM
@printf("%1.15e\n",sqrt(NRM))
A[:,ii+1:ii+2] = A[:,ii+1:ii+2]*T
A[ii+1:ii+2,:] = Tinv*A[ii+1:ii+2,:]
Ttot[:,ii+1:ii+2] = Ttot[:,ii+1:ii+2]*T

end

if abs(A[m,m-1]) < 1e-16*na
  m -= 1
end

if m == 1
  break
end

end

norm(Ttot)*norm(inv(Ttot))

#D = zeros(T,2*n,1)
#for ii=1:n
#  D[2*ii-1] = 1
#end
#D[end-1] = (-1)^(n-1)
#ITS = zeros(Int,n-1)
#
#@time Eiscor.UnitaryFA(Q,D,ITS)
#print("\n done!\n")
#
#@profile Eiscor.UnitaryFA(Q,D,ITS)
#Profile.print()

