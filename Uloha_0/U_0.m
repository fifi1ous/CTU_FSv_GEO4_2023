clc; clear; format long G

a=17; b=11;

A=[a,2,3;2,4,5;3,5,b];

%%
fid=fopen('Vysledky_Roucka.txt','w');
%% Determinant
D_MAT=det(A);
D_VYP=(A(1)*A(5)*A(9)+A(4)*A(8)*A(3)+A(7)*A(2)*A(6))-(A(7)*A(5)*A(3)+A(1)*A(8)*A(6)+A(9)*A(2)*A(4));
fprintf(fid,'Determinant matice vypočítaný z Matlabu: %.0f',D_MAT);
fprintf(fid,'\nDeterminant matice vypočítaný ručně: %.0f',D_VYP);
%% Stopa matice
T_MAT=trace(A);
T_VYP=A(1)+A(5)+A(9);
fprintf(fid,'\nStopa matice vypočítaný z Matlabu: %.0f',T_MAT);
fprintf(fid,'\nStopa matice vypočítaný ručně: %.0f',T_VYP);
%% Inverzní matice 
I_MAT=A^(-1);
if D_VYP==0
    error('Inverzní matice nejde vypočítat determinant matice je roven 0');
end
for n=1:size(A,1)
    for m=1:size(A,2)
        M = A;
        M(n,:) = [];
        M(:,m) = [];
        if rem(n+m,2)==1
            A1(n,m)=det(M)*(-1);
        else
            A1(n,m)=det(M)*(1);
        end
    end
end
A1=A1';
I_VYP=A1/D_VYP;
fprintf(fid,'\n\nInverzní matice vypočítaný z Matlabu:\n');
fprintf(fid,'%7.4f %7.4f %7.4f\n',I_MAT);
fprintf(fid,'\nInverzní matice vypočítaný ručně:\n');
fprintf(fid,'%7.4f %7.4f %7.4f\n',I_VYP);
%% Vlastní císla matice
VL_MAT=eig(A);
syms y
KOEF=det(A-y*eye(3));
subs(KOEF,y,1);
KOEF_MAT=poly(A);
KOEF_VYP=[-1,+32,-261,303];
koef=KOEF_VYP;
m=1000;
x1=1.3;
f1=koef(1)*x1^3 + koef(2)*x1^2 + koef(3)*x1 + koef(4);
f1d=koef(1)*3*x1^2 + koef(2)*2*x1 + koef(3);
x2=11.3;
f2=koef(1)*x2^3 + koef(2)*x2^2 + koef(3)*x2 + koef(4);
f2d=koef(1)*3*x2^2 + koef(2)*2*x2 + koef(3);
x3=19.3;
f3=koef(1)*x3^3 + koef(2)*x3^2 + koef(3)*x3 + koef(4);
f3d=koef(1)*3*x3^2 + koef(2)*2*x3 + koef(3);

for n=1:m
    VL_VYP(1)=x1-f1/f1d;
    VL_VYP(2)=x2-f2/f2d;
    VL_VYP(3)=x3-f3/f3d;
end

fprintf(fid,'\nVlastní čísla matice vypočítaný z Matlabu: %8.3f %8.3f %8.3f',VL_MAT);
fprintf(fid,'\nVlastní čísla matice vypočítaný ručně: %8.3f %8.3f %8.3f\n',VL_VYP);
%%
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
%%Podobnostní transformace
ss=SS(1:3,2:3); ib=IB(:,2:3);

xT=mean(ss(:,2)); yT=mean(ss(:,1));
XT=mean(ib(:,2)); YT=mean(ib(:,1));

xr=ss(:,2)-xT; yr=ss(:,1)-yT;
Xr=ib(:,2)-XT; Yr=ib(:,1)-YT;

lam1=(xr'*Xr+yr'*Yr)/(xr'*xr+yr'*yr);
lam2=(xr'*Yr-yr'*Xr)/(xr'*xr+yr'*yr);

q=sqrt(lam1^2+lam2^2);  w=atan2(lam2,lam1);

X0=XT-lam1*xT+lam2*yT;
Y0=YT-lam1*yT-lam2*xT;

x=SS(:,3); y=SS(:,2);

XX=X0+lam1*x-lam2*y;
YY=Y0+lam1*y+lam2*x;

Souradnice_transformace=[SS(:,1),YY,XX];
Plocha_h=polyarea(XX(4:end),YY(4:end));

%% Afinní transformace
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

Plocha_A=polyarea(Y(4:end),X(4:end));

delt=Plocha_A-Plocha_h;
fprintf(fid,'\nVýměra vypočítaná pomocí podobnostní transformace: %d m^2',round(Plocha_h));
fprintf(fid,'\nVýměra vypočítaná pomocí afinní transformace: %d m^2',round(Plocha_A));
fprintf(fid,'\nRozdíl vypočítaných výměr: %d m^2', round(delt));
fclose(fid);
