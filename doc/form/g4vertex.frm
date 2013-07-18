Off Statistics;

*
* Declarations:
*

autodeclare I a,b,A,B,C,D,E,mu;
Set a:a1,a2,a3,a4,a5;
Set A:A1,A2,A3,A4,A5;
Set B:B1,B2,B3,B4,B5;
Set C:C1,C2,C3,C4,C5;
Set D:D1,D2,D3,D4,D5;
Set E:E1,E2,E3,E4,E5;
Set mu:mu1,mu2,mu3,mu4,mu5;
T ADvertex,t,f,G,G1,G2,G3,G4,ddmin,YM;
S N,i,j,k,l;
F ym;

*
* Recursive relation definitions in the adjoint representation:
*

L ADglu1=ADvertex(a1,mu1,a2,mu2,a3,mu3,a4,mu4)*G2(a2,mu2)*G3(a3,mu3)*G4(a4,mu4);
L ADglu2=ADvertex(a2,mu2,a1,mu1,a3,mu3,a4,mu4)*G1(a2,mu2)*G3(a3,mu3)*G4(a4,mu4);
L ADglu3=ADvertex(a3,mu3,a2,mu2,a1,mu1,a4,mu4)*G1(a3,mu3)*G2(a2,mu2)*G4(a4,mu4);
L ADglu4=ADvertex(a4,mu4,a2,mu2,a3,mu3,a1,mu1)*G1(a4,mu4)*G2(a2,mu2)*G3(a3,mu3);

*
* Vertex tensor definition in the adjoint representation:
*

L CFvertex=ADvertex(a1,mu1,a2,mu2,a3,mu3,a4,mu4)*(2*t(a1,A1,B1))*(2*t(a2,A2,B2))*(2*t(a3,A3,B3))*(2*t(a4,A4,B4));

*
* Replacing adjoint octets with colour-flow nonets:
*

id G?(b1?a[i],mu1?)=2*t(b1,B[i],A[i])*G(A[i],B[i],mu1);

.sort

*
* Defining the recursive relations in the colour-flow representation:
*

L CFglu1=t(a1,A1,B1)*ADglu1;
L CFglu2=t(a1,A1,B1)*ADglu2;
L CFglu3=t(a1,A1,B1)*ADglu3;
L CFglu4=t(a1,A1,B1)*ADglu4;

*
* Insertion of the vertex tensor in the adjoint representation:
*

id
ADvertex(a1?,mu1?,a2?,mu2?,a3?,mu3?,a4?,mu4?) =
f(a1,a2,a5)*f(a3,a4,a5)*ddmin(mu1,mu3,mu2,mu4) +
f(a1,a3,a5)*f(a4,a2,a5)*ddmin(mu1,mu4,mu2,mu3) +
f(a1,a4,a5)*f(a2,a3,a5)*ddmin(mu1,mu2,mu3,mu4);

AntiBracket ddmin;
*
* Replacing structure constants and generator matrices by their colour-flow
* expressions:
*

id f(b1?a[i],b2?,b3?)=2*i_*(t(b1,C[i],D[i])*t(b2,D[i],E[i])-t(b2,C[i],D[i])*t(b1,D[i],E[i]))*t(b3,E[i],C[i]);

id t(a1?,A1?,B1?)*t(a1?,A2?,B2?)=1/2*(d_(A1,B2)*d_(A2,B1)-1/N*d_(A1,B1)*d_(A2,B2));


*id dd1(mu1?,mu2?,mu3?,mu4?)=1/2*YM(mu1,mu2,mu3,mu4)-1/2*YM(mu1,mu3,mu2,mu4)+dd3(mu1,mu2,mu3,mu4);

*id dd2(mu1?,mu2?,mu3?,mu4?)=YM(mu1,mu3,mu2,mu4)-dd3(mu1,mu2,mu3,mu4);
.sort
Collect ym;
.sort
*
* Colour-flow vertex tensor:
*

print CFvertex;

*
* First nonet recursive relation:
*

print CFglu1;

*
* Second nonet recursive relation:
*

print CFglu2;

*
* Third nonet recursive relation:
*

print CFglu3;

*
* Fourth nonet recursive relation:
*

print CFglu4;

.end
