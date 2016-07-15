module Rotation

function Rot3Vec3Gen(AR,AI,B)

  NRM = sqrt(AR*AR+AI*AI+B*B)
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
     
     # compute first rotation
     c4r,c4i,s4,nrm = Rot3Vec3Gen(c1r*c2r + c1i*c2i, -c1i*c2r + c1r*c2i, s2)
     
     # compute second rotation
     c5r,c5i,s5,nrm = Rot3Vec3Gen(c1r*c3r - c1i*c3i, c1r*c3i + c1i*c3r, 0)
     
     # compute third rotation
     c6r,c6i,s6,nrm = Rot3Vec3Gen(c2r*c4r + c2i*c4i + c1r*s2*s4, 
                                  c2i*c4r - c2r*c4i + c1i*s2*s4,
          s1*s2*s5 - (c5r*c2r + c5i*c2i)*s4 + s2*(c1r*(c4r*c5r + c4i*c5i) 
          + c1i*(c4r*c5i - c4i*c5r)))
     
  else
     # compute first rotation
     c4r,c4i,s4,nrm = Rot3Vec3Gen(s1*c3r + (c1r*c2r + c1i*c2i)*s3, 
                                  s1*c3i + (-c1i*c2r + c1r*c2i)*s3, s2*s3)
     
     # compute second rotation
     c5r,c5i,s5,nrm = Rot3Vec3Gen(c1r*c3r - c1i*c3i - s1*c2r*s3,
                                  c1r*c3i + c1i*c3r - s1*c2i*s3, nrm)
     
     # compute third rotation
     c6r,c6i,s6,nrm = Rot3Vec3Gen(c2r*c4r + c2i*c4i + c1r*s2*s4,
                                  c2i*c4r - c2r*c4i + c1i*s2*s4,
          s1*s2*s5 - (c5r*c2r + c5i*c2i)*s4 + s2*(c1r*(c4r*c5r + c4i*c5i) 
          + c1i*(c4r*c5i - c4i*c5r)))

  end

  # store output
  G1[1] = c5r
  G1[2] = c5i
  G1[3] = s5
  G2[1] = c6r
  G2[2] = c6i
  G2[3] = s6
  G3[1] = c4r
  G3[2] = c4i
  G3[3] = s4
  
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
  
        # update G1
        c2r = c3r*phr + c3i*phi
        c2i = -c3r*phi + c3i*phr
        s2 = s3r*phr + s3i*phi
        c2r,c2i,s2,nrm = Rot3Vec3Gen(c2r,c2i,s2)
        G1[1] = c2r
        G1[2] = c2i
        G1[3] = s2
  
        # set G2
        G2[1] = phr
        G2[2] = phi
        G2[3] = 0
    
    # store product in G1 and diagonal in G2
    else
  
        # update G2
        c2r = c3r*phr - c3i*phi
        c2i = c3r*phi + c3i*phr
        s2 = s3r*phr + s3i*phi
        c2r,c2i,s2,nrm = Rot3Vec3Gen(c2r,c2i,s2)
        G2[1] = c2r
        G2[2] = c2i
        G2[3] = s2
  
        # set G1
        G1[1] = phr
        G1[2] = -phi
        G1[3] = 0
  
    end

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
  G[1] = c1r
  G[2] = c1i
  G[3] = s1
    
  # set D
  D[1] = d2r
  D[2] = d2i
  D[3] = d1r
  D[4] = d1i

end 

end
