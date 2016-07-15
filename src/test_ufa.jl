reload("Eiscor")

n = 2^7

#set_bigfloat_precision(128)
#T = BigFloat
T = Float64
#T = Float32
#T = Float16

Q = zeros(T,3*(n-1),1)
for ii=1:n-1
  Q[3*ii] = 1
end
D = zeros(T,2*n,1)
for ii=1:n
  D[2*ii-1] = 1
end
D[end-1] = (-1)^(n-1)
ITS = zeros(Int,n-1)

@time Eiscor.UnitaryFA(Q,D,ITS)
print("\n done!\n")

#@profile Eiscor.UnitaryFA(Q,D,ITS)
#Profile.print()

