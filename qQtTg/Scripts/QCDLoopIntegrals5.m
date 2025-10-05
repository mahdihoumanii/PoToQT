(* ::Package:: *)

EpsExpand[expr_, order_] := (Series[expr, {Epsilon, 0, order}]// Normal)//. Log[ScaleMu]-> Log[ScaleMu^2]/2/.Log[ScaleMu^2]-> Log[ScaleMu^2/mt^2]+ Log[mt^2]/. Log[a_]:> Log[a// Together// Factor]// Collect[#,{Epsilon, Log[__], PolyLog[2,__]}, Simplify[#/.eps-> 0]&]&


NumEps[expr_, num_] := N[EpsExpand[expr, 0]/. num, 10]


rGamma = Gamma[1-Epsilon]^2Gamma[1+Epsilon]/Gamma[1-2 Epsilon];
betar[s_]:= Sqrt[1- 4 mt^2/s]
Fac = Exp[-EulerGamma Epsilon]/rGamma;
xr[s_]:= ( betar[s]-1 )/(betar[s]+1)
(*MathCalLi2[x_, y_]:= PolyLog[2, 1- x y] + Log[1- x y ] (Log[x y] - Log[x] - Log[y])*)
ContinuedDiLogReplace[expr_]:= expr /. MathCalLi2[x_,y_]:>  PolyLog[2, 1- x y] + Log[1- x y ] (Log[x y] - Log[x] - Log[y])
KallenExpand[expr_]:= expr /. Kallen[a_,b_,c_] :>  a^2-2 a b+b^2-2 a c-2 b c+c^2


repA0 = {A0[mt^2]-> - ScaleMu^(2 Epsilon) Gamma[-1+ Epsilon]/rGamma (mt^2 )^(1- Epsilon)};


CzakonB0[s_,mt^2,mt^2]:= Module[{x},
x = (Sqrt[1- 4 mt^2/(s+ I eps)]-1)/(Sqrt[1- 4 mt^2/(s+ I eps)]+1); 

(1/ Epsilon + 
2 - Log[mt^2/ScaleMu^2] + (1+x)/(1-x) Log[x] + 
Epsilon*( 
4- 2 Log[mt^2/ScaleMu^2]+ 1/2 Log[mt^2/ScaleMu^2]^2 + Zeta[2]/2 + (1+x)/(1-x) ((2- Log[mt^2/ScaleMu^2] + 1/2 Log[x] - 2 Log[1+x])Log[x] - 2 PolyLog[2,-x] - Zeta[2])
)
)
]


repB0 = {
B0[s_, 0,0]:> (ScaleMu^2/(-s- I eps))^Epsilon(1/Epsilon/(1-2 Epsilon)),
B0[mt^2, 0, mt^2]:> Gamma[1+ Epsilon]/rGamma (ScaleMu^2/mt^2)^Epsilon (1/Epsilon/(1-2 Epsilon)), 
B0[s_, 0, mt^2]:> 
                (ScaleMu^2/mt^2)^Epsilon (
                1/ Epsilon 	
				+ 2+ (mt^2-s)/s Log[(mt^2-(s+ I eps))/mt^2]
				+ Epsilon (Pi^2/6 + 4 +(mt^2-s)/2/s ( 4 Log[(mt^2-(s+I eps))/mt^2]- Log[(mt^2-(s+ I eps))/mt^2]^2 + 2 PolyLog[2, -(s+ I eps)/(mt^2-(s+ I eps))]))
											), 
B0[s_, mt^2, mt^2]:> Fac * CzakonB0[s,mt^2,mt^2]

};


Symmetries = {C0[0 , s_, 0, 0,0,0]:> C0[0, 0, s , 0,0,0],D0[0, mt^2, mt^2, p3sq_, s12_, s23_, 0,0,mt^2, 0]:>  D0[0, p3sq, mt^2, mt^2, s23, s12, 0,0,0,mt^2], 
D0[0, mt^2, mt^2, psq_, s12_, s23_, mt^2, mt^2, 0, mt^2]:>  D0[mt^2, 0, psq, mt^2,s12,s23, 0, mt^2, mt^2, mt^2],  
D0[0, mt^2, 0, s_, s12_, s23_, mt^2, mt^2, 0,0]:> D0[0, mt^2, 0, s, s23, s12, 0, 0, mt^2, mt^2], C0[mt^2,s_, 0, 0, mt^2, 0] :> C0[0, mt^2, s, 0,0,mt^2],C0[0,s_, mt^2,0,0,mt^2] :>  C0[0,mt^2,s,0,0,mt^2] , 
C0[mt^2,mt^2,s_,mt^2,0,mt^2] :> C0[mt^2,s, mt^2, 0, mt^2, mt^2]
}; 


repD0 = { 
D0[0,0,0,0, s12_, s23_, 0,0,0,0]:> (*rGamma*) ScaleMu^(2 Epsilon) /s12/s23 ( 2/ Epsilon^2 (( -s12 - I eps)^(- Epsilon) + (-s23 - I eps)^(- Epsilon))
																		- Log[(-s12 - I eps)/(-s23 - I eps)]^2 - Pi^2 ), 




D0[0,0,0, p4sq_, s12_, s23_, 0,0,0,0]:> ReplaceAll[(*rGamma*) ScaleMu^(2 Epsilon)/s12/s23  * ( 
2/Epsilon^2 ((-s12)^(- Epsilon) + (-s23)^(- Epsilon) - (-p4sq )^(- Epsilon))
- 2 PolyLog[2, 1- p4sq/s12] - 2 PolyLog[2, 1- p4sq/s23] - Log[(-s12)/(-s23)]^2- Pi^2/3 
), {s12-> s12+ I eps, s23-> s23+ I eps, p4sq -> p4sq+ I eps}],




D0[0,0,mt^2, p4sq_, s12_, s23_, 0,0,0, mt^2]:> (ScaleMu^2/mt^2)^Epsilon * (*rGamma*)1 /s12/(s23- mt^2) *(

3/2 1/Epsilon^2 - 
1/Epsilon ( 2 Log[1- (s23+I eps)/mt^2] + Log[-(s12+ I eps)/ mt^2] - Log[1- (p4sq+ I eps)/ mt^2])- 
2 PolyLog[2, (s23 - p4sq)/ (s23+ I eps- mt^2)] + 2 Log[-(s12+ I eps)/mt^2] Log[1- (s23+ I eps)/mt^2] - Log[1- (p4sq+ I eps)/mt^2]^2 - 5 Pi^2/12

)
, 




D0[0, mt^2, p3sq_, mt^2, s12_, s23_, 0,0, mt^2, mt^2]:> ReplaceAll[
(ScaleMu^2/mt^2)^Epsilon /(s12- mt^2)/(s23- mt^2)* 
(
1/Epsilon^2 - 

1/Epsilon ( Log[1- s23/mt^2] + Log[1-s12/mt^2])- Log[xr[p3sq]]^2 + 

2 Log[1- s12/mt^2] Log[1- s23/mt^2] - Pi^2/2
)
,{p3sq-> p3sq+ I eps, s12-> s12 + I eps, s23 -> s23 + I eps}]





,
D0[0, p2sq_, p3sq_, mt^2, s12_, s23_, 0, 0, 0, mt^2] :> 
  (*rGamma*)1/(s12 (s23 - mt^2)) * (
    1/(2 Epsilon^2)
    - (1/Epsilon) * Log[((s12+ I eps)/(p2sq+ I eps)) * ((mt^2 - (s23+ I eps))/(ScaleMu mt))]
    + PolyLog[2, 1 + ((mt^2 - (p3sq+ I eps)) (mt^2 - (s23+ I eps)))/(mt^2 (p2sq+ I eps))]
    + 2 PolyLog[2, 1 - (s12+ I eps)/(p2sq+ I eps)]
    + Pi^2/12
    + Log[((s12+ I eps)/(p2sq+ I eps)) * ((mt^2 - s23- I eps)/(ScaleMu mt))]^2  
  ), 
 
 
  D0[0, mt^2, 0, p4sq_, s12_, s23_, 0, 0, mt^2, mt^2] :> 

  (ScaleMu^2/mt^2)^Epsilon * (1/((s12 - mt^2) (s23 - mt^2))) * (
    1/(2 Epsilon^2)
    - (1/Epsilon) * (
      Log[1 - (s12+ I eps)/mt^2]
      + Log[1 - (s23+ I eps)/mt^2]
      - Log[1 - (p4sq+I eps)/mt^2]
    )
    - 2 PolyLog[2, (s23 - p4sq)/(s23+ I eps - mt^2)]
    - 2 PolyLog[2, (s12 - p4sq)/(s12 + I eps- mt^2)]
    + 2 Log[1 - (s12+ I eps)/mt^2] * Log[1 - (s23+ I eps)/mt^2]
    - Log[1 - (p4sq+ I eps)/mt^2]^2
    - Pi^2/12
  )
 , 
 
 
 
 D0[mt^2, 0, p3sq_, mt^2, s12_, s23_, 0, mt^2, mt^2, mt^2] :> 

  (ScaleMu^2/mt^2)^Epsilon/ (s12- mt^2)/s23 /betar[s23+ I eps] * 
  
  ( 1/Epsilon Log[xr[s23+ I eps]] 
  
  - 2 Sum[MathCalLi2[xr[s23+ I eps], xr[p3sq+ I eps]^rho], {rho, {-1,1}}]
  
  - PolyLog[2, xr[s23+ I eps]^2] -  2 Log[xr[s23+ I eps]] Log[1- xr[s23+ I eps]^2]
  
  - 2 Log[xr[s23+ I eps]]Log[1- (s12+ I eps)/mt^2] - Log[xr[p3sq+ I eps]]^2 + Pi^2/6
  
  
  )/. eps-> 0

 
 

};


repC0 = {
C0[0,0, s_, 0,0,0]:> ScaleMu^(2 Epsilon)/Epsilon^2 ((-s-I eps)^(- Epsilon)/s)

,
 C0[mt^2,s_, mt^2, 0, mt^2, mt^2]:> ReplaceAll[-((Sqrt[-((4*mt^2 - s)*s)]*Log[(2*mt^2 - s + Sqrt[s*(-4*mt^2 + s)])/(2*mt^2)])/(Epsilon*(4*mt^2 - s)*s)) + 
 (Sqrt[-((4*mt^2 - s)*s)]*(Pi^2 - 3*Log[(2*mt^2 - s + Sqrt[s*(-4*mt^2 + s)])/(2*mt^2)]*Log[-1/2*(mt^2*(-2*mt^2 + s + Sqrt[s*(-4*mt^2 + s)]))/(-4*mt^2 + s)^2] - 
    6*Log[(2*mt^2 - s + Sqrt[s*(-4*mt^2 + s)])/(2*mt^2)]*Log[ScaleMu^2/mt^2] + 12*PolyLog[2, (-2*mt^2 + s - Sqrt[s*(-4*mt^2 + s)])/(2*mt^2)]))/(6*(4*mt^2 - s)*s), s-> s + I eps]
,
C0[0, mt^2, p2sq_, 0,0,mt^2]:> ReplaceAll[
(ScaleMu^2/mt^2)^Epsilon /(p2sq - mt^2) * (
1/2/Epsilon^2 + 

1/Epsilon Log[mt^2/(mt^2-p2sq)] + 

Pi^2/12 + 1/2 Log[mt^2/(mt^2 -p2sq)]^2 - PolyLog[2, -p2sq/(mt^2-p2sq)]

)

, {p2sq ->  p2sq + I eps}]


};


repScalarIntegrals= {

C0[mt^2, mt^2, s12_, 0, mt^2, 0]:> ReplaceAll[(4*Pi^2 + 3*Log[1 + ((-1 + Sqrt[1 - (4*mt^2)/s12])*s12)/(2*mt^2)]^2 + 12*PolyLog[2, 1 + ((-1 + Sqrt[1 - (4*mt^2)/s12])*s12)/(2*mt^2)])/(6*s12*Sqrt[(-4*mt^2 + s12)/s12]), {s12-> s12+ I eps}],

C0[0, s12_,s23_, mt^2, mt^2, 0] :> ReplaceAll[((Log[-(mt^2/s12)]^2 - Log[-(mt^2/s23)]^2 + 2*PolyLog[2, mt^2/s12] - 2*PolyLog[2, mt^2/s23])/(2*(s12 - s23))), {s12-> s12+ I eps, s23 -> s23 + I eps }], 

 C0[0, s12_, s23_, mt^2, mt^2, mt^2]:> ReplaceAll[((Log[(2*mt^2 - s12 + Sqrt[s12*(-4*mt^2 + s12)])/(2*mt^2)] - Log[(2*mt^2 - s23 + Sqrt[s23*(-4*mt^2 + s23)])/(2*mt^2)])*
  (Log[(2*mt^2 - s12 + Sqrt[s12*(-4*mt^2 + s12)])/(2*mt^2)] + Log[(2*mt^2 - s23 + Sqrt[s23*(-4*mt^2 + s23)])/(2*mt^2)]))/(2*(s12 - s23)), {s12-> s12+ I eps, s23 -> s23 + I eps }],

  D0[0, mt^2, p3sq_, mt^2, s12_, s23_, mt^2, mt^2, 0, 0]:> 1/2 (-((-2 Log[(-mt^2+s12)/p3sq] Log[(mt^4+mt^2 (2 p3sq-s12-s23)+s12 s23+Sqrt[mt^2-s12] Sqrt[mt^2-s23] Sqrt[mt^4+mt^2 (4 p3sq-s12-s23)+s12 s23])/(2 mt^2 p3sq)]+2 Log[mt^2/(mt^2-s23)] Log[(mt^4+mt^2 (2 p3sq-s12-s23)+s12 s23+Sqrt[mt^2-s12] Sqrt[mt^2-s23] Sqrt[mt^4+mt^2 (4 p3sq-s12-s23)+s12 s23])/(2 mt^2 p3sq)]-Log[(mt^4+mt^2 (2 p3sq-s12-s23)+s12 s23+Sqrt[mt^2-s12] Sqrt[mt^2-s23] Sqrt[mt^4+mt^2 (4 p3sq-s12-s23)+s12 s23])/(2 mt^2 p3sq)]^2-4 PolyLog[2,-((mt^4+s12 s23-mt^2 (s12+s23)+Sqrt[mt^2-s12] Sqrt[mt^2-s23] Sqrt[mt^4+mt^2 (4 p3sq-s12-s23)+s12 s23])/(2 mt^2 p3sq))])/(Sqrt[mt^2-s12] Sqrt[mt^2-s23] Sqrt[mt^4+4 mt^2 p3sq-mt^2 s12-mt^2 s23+s12 s23]))+(-2 Log[(-mt^2+s12)/p3sq] Log[(mt^4+mt^2 (2 p3sq-s12-s23)+s12 s23-Sqrt[mt^2-s12] Sqrt[mt^2-s23] Sqrt[mt^4+mt^2 (4 p3sq-s12-s23)+s12 s23])/(2 mt^2 p3sq)]+2 Log[mt^2/(mt^2-s23)] Log[(mt^4+mt^2 (2 p3sq-s12-s23)+s12 s23-Sqrt[mt^2-s12] Sqrt[mt^2-s23] Sqrt[mt^4+mt^2 (4 p3sq-s12-s23)+s12 s23])/(2 mt^2 p3sq)]-Log[(mt^4+mt^2 (2 p3sq-s12-s23)+s12 s23-Sqrt[mt^2-s12] Sqrt[mt^2-s23] Sqrt[mt^4+mt^2 (4 p3sq-s12-s23)+s12 s23])/(2 mt^2 p3sq)]^2-4 PolyLog[2,(-mt^4-s12 s23+mt^2 (s12+s23)+Sqrt[mt^2-s12] Sqrt[mt^2-s23] Sqrt[mt^4+mt^2 (4 p3sq-s12-s23)+s12 s23])/(2 mt^2 p3sq)])/(Sqrt[mt^2-s12] Sqrt[mt^2-s23] Sqrt[mt^4+4 mt^2 p3sq-mt^2 s12-mt^2 s23+s12 s23]))
 ,
C0[mt^2, s12_, s23_, 0, mt^2, mt^2]  :> -(\[Pi]^2/(6 Sqrt[Kallen[mt^2,s12,s23]]))-Log[(mt^2-s23)/Sqrt[Kallen[mt^2,s12,s23]]]^2/(2 Sqrt[Kallen[mt^2,s12,s23]])+PolyLog[2,-((-2 mt^2 (mt^2-s23)+2 mt^2 Sqrt[Kallen[mt^2,s12,s23]])/(2 mt^2 (mt^2-s23)))]/Sqrt[Kallen[mt^2,s12,s23]]-PolyLog[2,(-s12 (-3 mt^2+s12-s23)-s12 Sqrt[Kallen[mt^2,s12,s23]])/(s12 (3 mt^2-s12+s23)-Sqrt[s12 (-4 mt^2+s12)] Sqrt[Kallen[mt^2,s12,s23]])]/Sqrt[Kallen[mt^2,s12,s23]]+PolyLog[2,(-s12 (-3 mt^2+s12-s23)+s12 Sqrt[Kallen[mt^2,s12,s23]])/(s12 (3 mt^2-s12+s23)-Sqrt[s12 (-4 mt^2+s12)] Sqrt[Kallen[mt^2,s12,s23]])]/Sqrt[Kallen[mt^2,s12,s23]]-PolyLog[2,(-s12 (-3 mt^2+s12-s23)-s12 Sqrt[Kallen[mt^2,s12,s23]])/(s12 (3 mt^2-s12+s23)+Sqrt[s12 (-4 mt^2+s12)] Sqrt[Kallen[mt^2,s12,s23]])]/Sqrt[Kallen[mt^2,s12,s23]]+PolyLog[2,(-s12 (-3 mt^2+s12-s23)+s12 Sqrt[Kallen[mt^2,s12,s23]])/(s12 (3 mt^2-s12+s23)+Sqrt[s12 (-4 mt^2+s12)] Sqrt[Kallen[mt^2,s12,s23]])]/Sqrt[Kallen[mt^2,s12,s23]]-PolyLog[2,((mt^2-s23) (mt^2-s12+s23)+(-mt^2-s23) Sqrt[Kallen[mt^2,s12,s23]])/((mt^2-s23) (mt^2-s12+s23)+(mt^2-s23) Sqrt[Kallen[mt^2,s12,s23]])]/Sqrt[Kallen[mt^2,s12,s23]]-PolyLog[2,((mt^2-s23) (mt^2-s12+s23)+(-mt^2-s23) Sqrt[Kallen[mt^2,s12,s23]])/((mt^2-s23) (mt^2-s12+s23)+(-mt^2+s23) Sqrt[Kallen[mt^2,s12,s23]])]/Sqrt[Kallen[mt^2,s12,s23]]+PolyLog[2,((mt^2-s23) (mt^2-s12+s23)+(-mt^2+s23) Sqrt[Kallen[mt^2,s12,s23]])/((mt^2-s23) (mt^2-s12+s23)+(mt^2-s23) Sqrt[Kallen[mt^2,s12,s23]])]/Sqrt[Kallen[mt^2,s12,s23]]-PolyLog[2,(-mt^2+s23+Sqrt[Kallen[mt^2,s12,s23]])/Sqrt[Kallen[mt^2,s12,s23]]]/Sqrt[Kallen[mt^2,s12,s23]],
  C0[mt^2,s12_,s45_,mt^2,0,0]:> (\[Pi]^2/(6 Sqrt[Kallen[mt^2,s12,s45]])+Log[s12/Sqrt[Kallen[mt^2,s12,s45]]]^2/(2 Sqrt[Kallen[mt^2,s12,s45]])-PolyLog[2,(2 mt^2 s12-2 mt^2 Sqrt[Kallen[mt^2,s12,s45]])/(2 mt^2 s12)]/Sqrt[Kallen[mt^2,s12,s45]]-PolyLog[2,(-s12 (mt^2+s12-s45)-s12 Sqrt[Kallen[mt^2,s12,s45]])/(-s12 (mt^2+s12-s45)+s12 Sqrt[Kallen[mt^2,s12,s45]])]/Sqrt[Kallen[mt^2,s12,s45]]+PolyLog[2,(-s12 (mt^2+s12-s45)+s12 Sqrt[Kallen[mt^2,s12,s45]])/(-s12 (mt^2+s12-s45)-s12 Sqrt[Kallen[mt^2,s12,s45]])]/Sqrt[Kallen[mt^2,s12,s45]]-PolyLog[2,(-mt^4+mt^2 s12+2 mt^2 s45+s12 s45-s45^2+(mt^2-s45) Sqrt[Kallen[mt^2,s12,s45]])/(-mt^4+mt^2 s12+2 mt^2 s45+s12 s45-s45^2+(-mt^2+s45) Sqrt[Kallen[mt^2,s12,s45]])]/Sqrt[Kallen[mt^2,s12,s45]]+PolyLog[2,(-mt^4+mt^2 s12+2 mt^2 s45+s12 s45-s45^2+(mt^2+s45) Sqrt[Kallen[mt^2,s12,s45]])/(-mt^4+mt^2 s12+2 mt^2 s45+s12 s45-s45^2+(mt^2-s45) Sqrt[Kallen[mt^2,s12,s45]])]/Sqrt[Kallen[mt^2,s12,s45]]+PolyLog[2,(-mt^4+mt^2 s12+2 mt^2 s45+s12 s45-s45^2+(mt^2+s45) Sqrt[Kallen[mt^2,s12,s45]])/(-mt^4+mt^2 s12+2 mt^2 s45+s12 s45-s45^2+(-mt^2+s45) Sqrt[Kallen[mt^2,s12,s45]])]/Sqrt[Kallen[mt^2,s12,s45]]+PolyLog[2,(-s12+Sqrt[Kallen[mt^2,s12,s45]])/Sqrt[Kallen[mt^2,s12,s45]]]/Sqrt[Kallen[mt^2,s12,s45]])
 
  };


Note


QCDLoopIntegrals[expr_] := expr /. Symmetries/. repA0/. repB0/. repC0/. repD0
