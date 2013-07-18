#define MAX "6";
autodeclare I a,b,A,B,C,D,E;
Set a:a1,a2,a3,a4,a5;
Set C:C1,C2,C3,C4,C5;
Set D:D1,D2,D3,D4,D5;
Set E:E1,E2,E3,E4,E5;
T t,f,G1,G2,G3,G4,psi1,psi2,psi3,psi4;
S N,i;

*
*tensor d_(a1,a2)
*
L delta=d_(a1,a2)*(2*t(a1,A1,B1))*(2*t(a2,A2,B2));
*
*tensor d_(a1,a2)d_(a3,a4)
*
L ddelta=d_(a1,a2)*d_(a3,a4)*(2*t(a1,A1,B1))*(2*t(a2,A2,B2))*(2*t(a3,A3,B3))*(2*t(a4,A4,B4));
*
*tensor d_(a1,a2)d_(a3,a4)+d_(a1,a4)d_(a3,a2)
*
L ddeltaplus=(d_(a1,a2)*d_(a3,a4)+d_(a1,a4)*d_(a3,a2))*(2*t(a1,A1,B1))*(2*t(a2,A2,B2))*(2*t(a3,A3,B3))*(2*t(a4,A4,B4));
*
*tensor d_(a1,a2)d_(a3,a4)-d_(a1,a4)d_(a3,a2)
*
L ddeltamin=(d_(a1,a2)*d_(a3,a4)-d_(a1,a4)*d_(a3,a2))*(2*t(a1,A1,B1))*(2*t(a2,A2,B2))*(2*t(a3,A3,B3))*(2*t(a4,A4,B4));
*
*ttensor
*
L ttensor=t(a1,A2,A3)*(2*t(a1,A1,B1));
*
*ttensor contraction
*
L ttcontr=t(a1,A1,A2)*t(a1,A3,A4);
*
*ttensor multiplication
*
L ttmult=t(a1,A3,A5)*t(a2,A5,A4)*(2*t(a1,A1,B1))*(2*t(a2,A2,B2));
*
*ttensor commutator
*
L ttcomm=(t(a1,A3,A5)*t(a2,A5,A4)-t(a2,A3,A5)*t(a1,A5,A4))*(2*t(a1,A1,B1))*(2*t(a2,A2,B2));
*
*ttensor anti-commutator
*
L ttanticomm=(t(a1,A3,A5)*t(a2,A5,A4)+t(a2,A3,A5)*t(a1,A5,A4))*(2*t(a1,A1,B1))*(2*t(a2,A2,B2));
*
*ftensor
*
L structconst=f(a1,a2,a3)*(2*t(a1,A1,B1))*(2*t(a2,A2,B2))*(2*t(a3,A3,B3));
*
*ftensor contraction 
*
L structconstcontr= f(a1,a2,a5)*f(a3,a4,a5)*(2*t(a1,A1,B1))*(2*t(a2,A2,B2))*(2*t(a3,A3,B3))*(2*t(a4,A4,B4));

id f(b1?a[i],b2?,b3?)=2*i_*(t(b1,C[i],D[i])*t(b2,D[i],E[i])-t(b2,C[i],D[i])*t(b1,D[i],E[i]))*t(b3,E[i],C[i]);

id t(a1?,A1?,B1?)*t(a1?,A2?,B2?)=1/2*(d_(A1,B2)*d_(A2,B1)-1/N*d_(A1,B1)*d_(A2,B2));

.sort
print delta;
print ddelta;
print ddeltaplus;
print ddeltamin;
print ttensor;
print ttcontr;
print ttmult;
print ttcomm;
print ttanticomm;
print structconst;
print structconstcontr;
.end

