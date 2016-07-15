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
  Rotation.Rot3Fusion(false,Binv,Q)

  # swap diag 
  Rotation.Rot3SwapDiag(D,B)

  # fuse D and Binv
  scl1 = complex(D[1],D[2])*complex(Binv[1],Binv[2])
  scl2 = complex(D[3],D[4])*complex(Binv[1],-Binv[2])
  D[1] = real(scl1)
  D[2] = imag(scl1)
  D[3] = real(scl2)
  D[4] = imag(scl2)

  B

end

function ChaseDown(Q,D,B)

  # turnover
  Rotation.Rot3Turnover(slice(Q,1:3),slice(Q,4:6),B)

  # swap diag
  Rotation.Rot3SwapDiag(D,B)

end

function EndChase(Q,D,B)

  # fuse
  Rotation.Rot3Fusion(true,Q,B)

  # fuse D and Binv
  scl1 = complex(D[1],D[2])*complex(B[1],B[2])
  scl2 = complex(D[3],D[4])*complex(B[1],-B[2])
  D[1] = real(scl1)
  D[2] = imag(scl1)
  D[3] = real(scl2)
  D[4] = imag(scl2)

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
  if abs(shift) > 0
    shift = shift/abs(shift)
  end

  # start chase
  B = StartChase(slice(Q,1:3),slice(D,1:4),shift)

  # loop for chasing
  for ii = 1:(n-2)

    # set indices
    iq1 = 3*ii-2
    iq2 = 3*ii+3
    id1 = 2*ii+1
    id2 = 2*ii+4

    # chase down
    ChaseDown(slice(Q,iq1:iq2),slice(D,id1:id2),B)   

  end

  # end chase
  EndChase(slice(Q,3*n-5:3*n-3),slice(D,2*n-3:2*n),B)

end

function DeflationCheck(Q,D)

  # compute n
  n = round(Int,length(D)/2)

  # set tolerance based on type
  tol = eps(typeof(D[1]))

  # initialize zero
  ZERO = 0

  # loop for chasing
  for ii = 1:(n-1)

    # get off-diagonal of Q
    nrm = abs(Q[3*(n-ii)])
    if ( nrm < tol )

      # set zero
      ZERO = max(0,n-ii)

      # fetch entries of Q
      qr = Q[3*(n-ii)-2]
      qi = Q[3*(n-ii)-1]

      # update entries of Q
      Q[3*(n-ii)-2] = 1
      Q[3*(n-ii)-1] = 0
      Q[3*(n-ii)] = 0

      # update entries of D
      scl1 = complex(qr,qi)*complex(D[2*ZERO-1],D[2*ZERO])
      scl2 = complex(qr,-qi)*complex(D[2*ZERO+1],D[2*ZERO+2])
      D[2*ZERO-1] = real(scl1)
      D[2*ZERO] = imag(scl1)
      D[2*ZERO+1] = real(scl2)
      D[2*ZERO+2] = imag(scl2)
 
      # break out
      break 

    end

  end

  ZERO

end

function UnitaryFA(Q,D,ITS)

  # get size 
  N = round(Int,length(D)/2)

  # initialize variables
  STR = 1
  STP = N-1
  ZERO = 0
  ITMAX = 20*N
  ITCNT = 0

  # iteration loop
  for kk = 1:ITMAX

    # check for completion
    if(STP <= 0) 
      break
    end 
    
    # check for deflation
    ZERO = DeflationCheck(slice(Q,(3*STR-2):(3*STP)),slice(D,(2*STR-1):(2*STP+2)))

    # if ZERO > 0 deflate
    if ZERO > 0
    
      # set itcnt
      ITS[STR+ZERO-1] = ITCNT
      ITCNT = 0

      # move STP if at bottom
      if ZERO == (STP-STR+1)
        STP -= 1
        ZERO = 0
        STR = 1
      # otherwise move STR
      else
        STR = STR + ZERO 
        ZERO = 0
      end
    
    # if greater than 2x2 chase a bulge
    else

      # perform singleshift iteration
      if ( ITCNT == 10 )
        SingleStep(true,slice(Q,(3*STR-2):(3*STP)),slice(D,(2*STR-1):(2*STP+2)))
      else
        SingleStep(false,slice(Q,(3*STR-2):(3*STP)),slice(D,(2*STR-1):(2*STP+2)))
      end     

      # update indices
      ITCNT += 1
 
    end
    
    # if ITMAX hit
    if (kk == ITMAX) 
    end
    
  end

end

end
