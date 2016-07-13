module Eiscor

include("Rotation.jl")

function StartChase(Q,D,shift)

  # first transformation
  x = [complex(Q[1],Q[2]); complex(Q[3],0)]
  x = x*complex(D[1],D[2])
  x[1] -= shift
  phase = x[2]/abs(x[2])
  x = x*conj(phase)
  cr,ci,s,nrm = Rotation.Rot3Vec3Gen(real(x[1]),imag(x[1]),real(x[2]))
  B = [cr; ci; s]

  # fuse
  Binv = [cr; -ci; -s]
  Binv,Q = Rotation.Rot3Fusion(false,Binv,Q)

  # swap diag 
  D,B = Rotation.Rot3SwapDiag(D,B)

  # fuse D and Binv
  scl1 = complex(D[1],D[2])*complex(Binv[1],Binv[2])
  scl2 = complex(D[3],D[4])*complex(Binv[1],-Binv[2])
  D = [real(scl1); imag(scl1); real(scl2); imag(scl2)]

  Q,D,B

end

function ChaseDown(Q,D,B)

  # turnover
  B,Q[1:3],Q[4:6] = Rotation.Rot3Turnover(Q[1:3],Q[4:6],B)

  # swap diag
  D,B = Rotation.Rot3SwapDiag(D,B)

  # return
  Q,D,B
  
end

function EndChase(Q,D,B)

  # fuse
  Q,B = Rotation.Rot3Fusion(true,Q,B)

  # fuse D and Binv
  scl1 = complex(D[1],D[2])*complex(B[1],B[2])
  scl2 = complex(D[3],D[4])*complex(B[1],-B[2])
  D = [real(scl1); imag(scl1); real(scl2); imag(scl2)]

  Q,D

end

function SingleStep(flag,Q,D)

  # compute n
  n = round(Int,length(D)/2)

  if (flag)
    # random shift
    shift = complex(randn(),randn())
  else
    # rayleigh quotient shift
    shift = complex(Q[end-2],-Q[end-1])*complex(D[end-1],D[end])
  end
print(" shift = ",shift,"\n")

  # start chase
  Q[1:3],D[1:4],B = StartChase(Q[1:3],D[1:4],shift)

  # loop for chasing
  for ii = 1:(n-2)

    # set indices
    iq1 = 3*ii-2
    iq2 = 3*ii+3
    id1 = 2*ii+1
    id2 = 2*ii+4

    # chase down
    Q[iq1:iq2],D[id1:id2],B = ChaseDown(Q[iq1:iq2],D[id1:id2],B)   

  end

  # end chase
  Q[end-2:end],D[end-3:end] = EndChase(Q[end-2:end],D[end-3:end],B)

  Q,D

end

end
