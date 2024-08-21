clc; clear; format long G
RAD=pi/200;
A=[563894.483 1056937.030];
B=[562162.433 1056937.030];
C=[563028.458 1058437.030];

P=[563028.457  1057437.030];


el=[0.7*RAD,0.7*RAD,0.7*RAD];          %přesnost měření v pořadí s_směrPA[mgon], s_směrPB[mgon], s_směrPC[mgon], 
ex=[30,30,30,30,30,30];                %přenost podkladu v pořadní sxA[mm],syA[mm],sxB[mm],syB[mm],sxC[mm],syC[mm]
%% nastavení
X=2;Y=1;
SAP=sqrt((A(1)-P(1))^2+(A(2)-P(2))^2);
SBP=sqrt((B(1)-P(1))^2+(B(2)-P(2))^2);
SCP=sqrt((C(1)-P(1))^2+(C(2)-P(2))^2);

EL=eye(length(el)).*el.^2;
EX=eye(length(ex)).*ex.^2;

A1=[((B(Y)-P(Y))/(SBP^2))-((A(Y)-P(Y))/(SAP^2)),-(B(X)-P(X))/(SBP^2)+(A(X)-P(X))/(SAP^2)
    ((C(Y)-P(Y))/(SCP^2))-((B(Y)-P(Y))/(SBP^2)),-(C(X)-P(X))/(SCP^2)+(B(X)-P(X))/(SBP^2)];


A2=[((A(Y)-P(Y))/(SAP^2)),-(A(X)-P(X))/(SAP^2),-((B(Y)-P(Y))/(SBP^2)),(B(X)-P(X))/(SBP^2),0,0
    0,0,((B(Y)-P(Y))/(SBP^2)),-(B(X)-P(X))/(SBP^2),-((C(Y)-P(Y))/(SCP^2)),(C(X)-P(X))/(SCP^2)];

D=[-1,1,0;0,-1,1];
K=A1^(-1)*D;
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


presnost=[sx1,sy1,sxy1,alf1,a1,b1;sx2,sy2,sxy2,alf2,a2,b2]     %úhel stočení elipsy chyb je jedno protože se jedná o kružnici
                                                               % má vyjít
                                                               % 0g