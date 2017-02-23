module RFR

function RFR2Vec2Gen(A,B)

  NRM = max(abs(A),abs(B))
  C = A/NRM
  S = B/NRM
  NRM = 1-C*C+S*S
  NRM = 1 + NRM*(1 - NRM/4)/2
  C = C/NRM
  S = S/NRM
  NRM = C*C + S*S

  C,S,NRM

end

end
