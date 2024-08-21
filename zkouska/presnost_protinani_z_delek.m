clc; clear; format long G
RAD=pi/200;
A=[563894.483 1056937.030];
B=[562162.433 1056937.030];

P=[563028.457 1057437.030];


el=[5,5];          %přesnost měření v pořadí s_délkaAP[mm], s_délkaBP[mm]
ex=[30,30,30,30];                %přenost podkladu v pořadní sxA[mm],syA[mm],sxB[mm],syB[mm]
%% nastavení
X=2;Y=1;
SAP=sqrt((A(1)-P(1))^2+(A(2)-P(2))^2);
SBP=sqrt((B(1)-P(1))^2+(B(2)-P(2))^2);
SAB=sqrt((A(1)-B(1))^2+(A(2)-B(2))^2);

EL=eye(length(el)).*el.^2;
EX=eye(length(ex)).*ex.^2;

A1=[-(A(X)-P(X))/(SAP),-((A(Y)-P(Y))/(SAP));-(B(X)-P(X))/(SBP),-((B(Y)-P(Y))/(SBP))];

A2=[(A(X)-P(X))/(SAP),((A(Y)-P(Y))/(SAP)),0,0
    0,0,(B(X)-P(X))/(SBP),((B(Y)-P(Y))/(SBP))];

K=A1^(-1);
L=A1^(-1)*A2;

EX1=K*EL*K';
c1=sqrt((EX1(1)-EX1(4))^2+4*EX1(2)^2);
a1=sqrt((EX1(1)+EX1(4)+c1)/2);
b1=sqrt((EX1(1)+EX1(4)-c1)/2);
alf1=atan2(2*EX1(2),EX1(1)-EX1(4))/2/RAD;
if alf1<0
    alf1=(alf1+400);
end
sx1=sqrt(EX1(1));sy1=sqrt(EX1(4));sxy1=sqrt((EX1(1)+EX1(4))/2);

EX2=K*EL*K'+L*EX*L';
c2=sqrt((EX2(1)-EX2(4))^2+4*EX2(2)^2);
a2=sqrt((EX2(1)+EX2(4)+c2)/2);
b2=sqrt((EX2(1)+EX2(4)-c2)/2);
alf2=atan2(2*EX2(2),EX2(1)-EX2(4))/2/RAD;
if alf2<0
    alf2=(alf2+400);
end
sx2=sqrt(EX2(1));sy2=sqrt(EX2(4));sxy2=sqrt((EX2(4)+EX2(1))/2);


presnost=[sx1,sy1,sxy1,alf1,a1,b1;sx2,sy2,sxy2,alf2,a2,b2]