module Rotation

function Rot3Vec3Gen(AR,AI,B)

  NRM = norm([AR;AI;B])
  CR = AR/NRM
  CI = AI/NRM
  S = B/NRM

  CR,CI,S,NRM

end


function Rot3Turnover(G1,G2,G3)

  # set local variables
  c1r = G1[1]
  c1i = G1[2]
  s1 = G1[3]
  c2r = G2[1]
  c2i = G2[2]
  s2 = G2[3]
  c3r = G3[1]
  c3i = G3[2]
  s3 = G3[3]

  if ((s1 == 0) & (s3 == 0)) 
     # the case s1=s3=0 is special
     # using the procedure for the generic case results in a rotation with
     # c=1 and s=0 and in a rotation with a not necessarily real sine
     
     # initialize c4r, c4i and s4
     T = [c1r*c2r + c1i*c2i; -c1i*c2r + c1r*c2i; s2]
     
     # compute first rotation
     c4r,c4i,s4,nrm = Rot3Vec3Gen(T[1],T[2],T[3])
     
     # initialize c5r, c5i and s5
     T = [c1r*c3r - c1i*c3i; c1r*c3i + c1i*c3r; 0]
     
     # compute second rotation
     c5r,c5i,s5,nrm = Rot3Vec3Gen(T[1],T[2],T[3])
     
     # initialize c6r, c6i and s6
     T = [c2r*c4r + c2i*c4i + c1r*s2*s4; c2i*c4r - c2r*c4i + c1i*s2*s4;
          s1*s2*s5 - (c5r*c2r + c5i*c2i)*s4 + s2*(c1r*(c4r*c5r + c4i*c5i) 
          + c1i*(c4r*c5i - c4i*c5r))]
     
     # compute third rotation
     c6r,c6i,s6,nrm = Rot3Vec3Gen(T[1],T[2],T[3])
     
  else
     # initialize c4r, c4i and s4
     T = [s1*c3r + (c1r*c2r + c1i*c2i)*s3; s1*c3i + (-c1i*c2r + c1r*c2i)*s3; s2*s3]
     
     # compute first rotation
     c4r,c4i,s4,nrm = Rot3Vec3Gen(T[1],T[2],T[3])
     
     # initialize c5r, c5i and s5
     T = [c1r*c3r - c1i*c3i - s1*c2r*s3; c1r*c3i + c1i*c3r - s1*c2i*s3; nrm]
     
     # compute second rotation
     c5r,c5i,s5,nrm = Rot3Vec3Gen(T[1],T[2],T[3])
     
     # initialize c6r, c6i and s6
     T = [c2r*c4r + c2i*c4i + c1r*s2*s4; c2i*c4r - c2r*c4i + c1i*s2*s4;
          s1*s2*s5 - (c5r*c2r + c5i*c2i)*s4 + s2*(c1r*(c4r*c5r + c4i*c5i) 
          + c1i*(c4r*c5i - c4i*c5r))]
     
     # compute third rotation
     c6r,c6i,s6,nrm = Rot3Vec3Gen(T[1],T[2],T[3])

  end

  # store output
  G4 = [c4r; c4i; s4]
  G5 = [c5r; c5i; s5]
  G6 = [c6r; c6i; s6]
  
  G4,G5,G6

  end

function Rot3Fusion(FLAG,G1,G2)

    # retrieve G1  
    c1r = G1[1]
    c1i = G1[2]
    s1 = G1[3]
       
    # retrieve G2  
    c2r = G2[1]
    c2i = G2[2]
    s2 = G2[3]
      
    # compute givens product
    c3r = c1r*c2r - c1i*c2i - s1*s2
    c3i = c1r*c2i + c1i*c2r
    s3r = s1*c2r + c1r*s2
    s3i = s1*c2i - c1i*s2
       
    # compute phase
    scl = complex(s3r,s3i)
    nrm = abs(scl)
    if nrm > 0
      phr = real(scl)/nrm
      phi = imag(scl)/nrm
    else
      phr = 1
      phi = 0
    end
  
    # store product in G1 and diagonal in G2
    if (FLAG) 
  
        # update G3
        c2r = c3r*phr + c3i*phi
        c2i = -c3r*phi + c3i*phr
        s2 = s3r*phr + s3i*phi
        c2r,c2i,s2,nrm = Rot3Vec3Gen(c2r,c2i,s2)
        G3 = [c2r; c2i; s2]
  
        # set G4
        G4 = [phr; phi; 0]
    
    # store product in G4 and diagonal in G3
    else
  
        # update G4
        c2r = c3r*phr - c3i*phi
        c2i = c3r*phi + c3i*phr
        s2 = s3r*phr + s3i*phi
        c2r,c2i,s2,nrm = Rot3Vec3Gen(c2r,c2i,s2)
        G4 = [c2r; c2i; s2]
  
        # set G3
        G3 = [phr; -phi; 0]
  
    end

    G3,G4
  
  end 

function Rot3SwapDiag(D,G)

  # set inputs
  c1r = G[1]
  c1i = G[2]
  s1 = G[3]
  
  # retrieve D
  d1r = D[1]
  d1i = D[2]
  d2r = D[3]
  d2i = D[4]  
  
  # pass through diagonal
  nrm = (d1r*d2r + d1i*d2i)*c1r - (-d1r*d2i + d1i*d2r)*c1i
  c1i = (d1r*d2r + d1i*d2i)*c1i + (-d1r*d2i + d1i*d2r)*c1r
  c1r = nrm

  # renormalize
  c1r,c1i,s1,nrm = Rot3Vec3Gen(c1r,c1i,s1) 
  Gout = [c1r; c1i; s1]
    
  # set D
  Dout = [d2r; d2i; d1r; d1i]

  Dout,Gout
  
end 

end
