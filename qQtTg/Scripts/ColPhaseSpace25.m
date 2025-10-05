(* ::Package:: *)

Momenta[beta_, costheta_, phi_,z_,lambda_] := 
Module[{p1,p, p2, p3, p4, K, n1,n2, p34, p34sq,
  sqp34sq, sintheta, boostp, kt,ktsq,pp1 },
  p =  (1/2)*{1, 1, 0,0};
  p1 = (1/2)*{1,-1, 0, 0};
  kt = {0,0,1,0};
  pp1 = p[[1]] p1[[1]]-p[[2]] p1[[2]]-p[[3]] p1[[3]]-p[[4]] p1[[4]];
  ktsq = -1;
  K = -((1-z)p - lambda*kt -lambda^2 ktsq/pp1/2/(1-z)*p1) ;
  p2 = (z)p + lambda*kt -lambda^2 ktsq/pp1/2/(z)*p1; 
  p34 = p1 + p2 - K;
  p34sq = p34[[1]]^2 - p34[[2]]^2 - p34[[3]]^2 - p34[[4]]^2;
  sqp34sq = Sqrt[p34sq];

  sintheta = Sqrt[1 - costheta^2]; 
  p3 = (1/2)*Flatten[{sqp34sq,Sqrt[p34sq -beta^2]*{Sin[phi]*sintheta,
       Cos[phi]*sintheta, costheta}}];
  p4 = Flatten[{p3[[1]],-p3[[{2, 3, 4}]]}];

  boostp = p34[[2]]*p3[[2]] + p34[[3]]*p3[[3]] + p34[[4]]*p3[[4]]; 

  p3[[{2, 3, 4}]] = (p3[[1]] + boostp/(p34[[1]] + sqp34sq))*p34[[{2, 3, 4}]] + 
                    sqp34sq*p3[[{2, 3, 4}]];
  p3[[1]] = p34[[1]]*p3[[1]] + boostp;
  p3 = p3/sqp34sq; 

  p4[[{2, 3, 4}]] = (p4[[1]] - boostp/(p34[[1]] + sqp34sq))*p34[[{2, 3, 4}]] + 
                    sqp34sq*p4[[{2, 3, 4}]];
  p4[[1]] = p34[[1]]*p4[[1]] - boostp;
  p4 = p4/sqp34sq; 

  Return[{p1, p2, p3, p4, K,p}]];


S[p_, q_] := p[[1]]*q[[1]] - p[[2]]*q[[2]] - p[[3]]*q[[3]] - p[[4]]*q[[4]]
S[p_]:= S[p,p];


SetAttributes[S, Orderless];

