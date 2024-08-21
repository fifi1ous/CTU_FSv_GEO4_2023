clc; clear; format long G
%% nastavení
IB=[501 740468.64 1057661.18
    502 740403.72 1057576.40
    503 739854.96 1057931.68];
SS=[501 1156.54 0966.75
    502 1105.89 0873.58
    503 0510.32 1149.57
    1 1067.88 1193.32
    2 0996.48 1138.00
    3 0993.44 1140.13
    4 0969.60 1175.76
    5 0996.98 1202.78
    6 1041.84 1233.86];
%%
pom=[1,1,1]; pom1=IB(:,2:end); pom=[SS(1:3,3),SS(1:3,2),pom'];
D=det(pom);

dyAC=SS(3,2)-SS(1,2); dxAC=SS(3,3)-SS(1,3);
dYAC=IB(3,2)-IB(1,2); dXAC=IB(3,3)-IB(1,3);

dyAB=SS(2,2)-SS(1,2); dxAB=SS(2,3)-SS(1,3);
dYAB=IB(2,2)-IB(1,2); dXAB=IB(2,3)-IB(1,3);

a1=(dyAC*dXAB-dyAB*dXAC)/D;
a2=(dxAB*dYAC-dxAC*dYAB)/D;
b1=(dxAC*dXAB-dxAB*dXAC)/D;
b2=(dyAC*dYAB-dyAB*dYAC)/D;

X0=IB(:,3)-a1*SS(1:3,3)+b1*SS(1:3,2);
Y0=IB(:,2)-b2*SS(1:3,3)-a2*SS(1:3,2);
X0=mean(X0);    Y0=mean(Y0);

X=X0+a1*SS(:,3)-b1*SS(:,2);
Y=Y0+b2*SS(:,3)+a2*SS(:,2);