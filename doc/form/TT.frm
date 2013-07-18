Off Statistics;

*
* Declarations:
*

autodeclare I a,b,A,B,C,D,E;
Set a:a1,a2,a3,a4,a5;
Set A:A1,A2,A3,A4,A5;
Set B:B1,B2,B3,B4,B5;
Set C:C1,C2,C3,C4,C5;
Set D:D1,D2,D3,D4,D5;
Set E:E1,E2,E3,E4,E5;
T ADvertex,t,f,phi1,phi2,phi3,phi4,G,G1,G2,G3,G4;
S N,i;

*
* Recursive relation definitions in the adjoint representation:
*

L ADglu1=ADvertex(a1,a2,A3,B4)*G2(a2)*phi3(A3)*phi4(B4);
L ADglu2=ADvertex(a2,a1,A3,B4)*G1(a2)*phi3(A3)*phi4(B4);
L ADphi3=ADvertex(a2,a3,A1,B4)*G1(a2)*G2(a3)*phi4(B4);
L ADphi4=ADvertex(a2,a3,A4,B1)*G1(a2)*G2(a3)*phi3(A4);

*
* Vertex tensor definition in the adjoint representation:
*

L CFvertex=ADvertex(a1,a2,A3,B4)*(2*t(a1,A1,B1))*(2*t(a2,A2,B2));

*
* Replacing adjoint octets with colour-flow nonets:
*

id G?(b1?a[i])=2*t(b1,B[i],A[i])*G(A[i],B[i]);

.sort

*
* Defining the recursive relations in the colour-flow representation:
*

L CFglu1=t(a1,A1,B1)*ADglu1;
L CFglu2=t(a1,A1,B1)*ADglu2;
L CFphi3=ADphi3;
L CFphi4=ADphi4;

*
* Insertion of the vertex tensor in the adjoint representation:
*

id ADvertex(a1?,a2?,A3?,B4?)=t(a1,A3,A5)*t(a2,A5,B4);

*
* Replacing structure constants and generator matrices by their colour-flow
* expressions:
*

id f(b1?a[i],b2?,b3?)=2*i_*(t(b1,C[i],D[i])*t(b2,D[i],E[i])-t(b2,C[i],D[i])*t(b1,D[i],E[i]))*t(b3,E[i],C[i]);

id t(a1?,A1?,B1?)*t(a1?,A2?,B2?)=1/2*(d_(A1,B2)*d_(A2,B1)-1/N*d_(A1,B1)*d_(A2,B2));

.sort
Bracket N;

*
* Colour-flow vertex tensor:
*

print CFvertex;

*
* First (nonet) recursive relation:
*

print CFglu1;

*
* Second (nonet) recursive relation:
*

print CFglu2;

*
* Third (triplet) recursive relation:
*

print CFphi3;

*
* Fourth (triplet) recursive relation:
*

print CFphi4;

.end
