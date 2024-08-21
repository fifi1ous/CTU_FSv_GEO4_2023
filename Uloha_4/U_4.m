clc; clear; format long G
RAD=pi/180;
%% Import dat
fid=fopen('Body_ETRS89.txt','r');
IB_ETRS=fscanf(fid,'%f %f %f %f %f %f %f %f',[8, inf])';IB_ETRS=sortrows(IB_ETRS,1);
fclose(fid);
fid=fopen('mereni_1.txt','r');
mer1=fscanf(fid,'%f %f %f %f %f %f %f %f',[8, inf])';
fclose(fid);
fid=fopen('mereni_stroj.txt','r');
mer_kontrola=fscanf(fid,'%f %f %f %f',[4, inf])';
fclose(fid);
fid=fopen('Body_S_JTSK.txt','r');
IB_S_JTSK=fscanf(fid,'%f %f %f %f',[4, inf])';IB_S_JTSK=sortrows(IB_S_JTSK,1);
fclose(fid);
%% a - Převod na geocetnrické souřadnice
BE_IB=(IB_ETRS(:,2)+IB_ETRS(:,3)/60+IB_ETRS(:,4)/3600)*RAD;
LE_IB=(IB_ETRS(:,5)+IB_ETRS(:,6)/60+IB_ETRS(:,7)/3600)*RAD;
HE_IB=IB_ETRS(:,8);

a_ETRS=6378137;  b_ETRS=6356752.314;
e_ETRS=0.006694380022901;
W_ETRS=sqrt(1-e_ETRS*(sin(BE_IB).^2));
N=a_ETRS./W_ETRS;

X_ETRS_GEO=cos(BE_IB).*cos(LE_IB).*(N+HE_IB);
Y_ETRS_GEO=cos(BE_IB).*sin(LE_IB).*(N+HE_IB);
Z_ETRS_GEO=(N*(1-e_ETRS)+HE_IB).*sin(BE_IB);
P_ETRS_GEO=[X_ETRS_GEO';Y_ETRS_GEO';Z_ETRS_GEO'];
%% b - Převod z S-JTSK na geocentrické
a_BESSEL=6377397.15508;
e_BESSEL=0.006674372230622;
e_bessel=sqrt(e_BESSEL);
fi0=(49+30/60)*RAD; S0=(78+30/60)*RAD; UQ=(59+42/60+42.6969/3600)*RAD; R_C=0.9999;
R_K=(a_BESSEL*sqrt(1-e_BESSEL))/(1-e_BESSEL*sin(fi0)^2);
r0=R_C*R_K*cot(S0);
alfa=sqrt(1+((e_BESSEL*cos(fi0)^4)/(1-e_BESSEL)));
U0=asin(sin(fi0)/alfa); n=sin(S0);n1=n;
k=((((1-e_bessel*sin(fi0))/(1+e_bessel*sin(fi0)))^(alfa*e_bessel/2))*(tan(fi0/2+pi/4)^alfa))/tan(U0/2+pi/4);

ro=(IB_S_JTSK(:,2).^2+IB_S_JTSK(:,3).^2).^(1/2);
eps=atan2(IB_S_JTSK(:,2),IB_S_JTSK(:,3));

S=2*atan(((r0./ro).^(1/n))*tan(S0/2+pi/4))-pi/2; D=eps/sin(S0);
U=asin(sin(UQ).*sin(S)-cos(UQ).*cos(S).*cos(D));
DV=asin((sin(D).*cos(S))./cos(U));

L=(24+50/60)*RAD-DV/alfa;
B=U;

for n=1:10
    B=2*atan((k.^(1./alfa)).*((1-e_bessel.*sin(B))./(1+e_bessel.*sin(B))).^((-e_bessel)./2).*tan(U./2+pi/4).^(1./alfa))-pi/2;   
end
B_BESSEL=B; L_BESSEL=L; H_BESSEL=IB_S_JTSK(:,4);

W_BESSEL=sqrt(1-e_BESSEL*(sin(B_BESSEL).^2));
N_BESSEL=a_BESSEL./W_BESSEL;

X_BESS_GEO=cos(B_BESSEL).*cos(L_BESSEL).*(N_BESSEL+H_BESSEL);
Y_BESS_GEO=cos(B_BESSEL).*sin(L_BESSEL).*(N_BESSEL+H_BESSEL);
Z_BESS_GEO=(N_BESSEL*(1-e_BESSEL)+H_BESSEL).*sin(B_BESSEL);
P_BESS_GEO=[X_BESS_GEO';Y_BESS_GEO';Z_BESS_GEO'];
%% c - tranformační klíč
X_ETRS_T=mean(X_ETRS_GEO);Y_ETRS_T=mean(Y_ETRS_GEO);Z_ETRS_T=mean(Z_ETRS_GEO);
X_BESS_T=mean(X_BESS_GEO);Y_BESS_T=mean(Y_BESS_GEO);Z_BESS_T=mean(Z_BESS_GEO);

X_ETRS_R=X_ETRS_GEO-X_ETRS_T;   Y_ETRS_R=Y_ETRS_GEO-Y_ETRS_T;   Z_ETRS_R=Z_ETRS_GEO-Z_ETRS_T;
X_BESS_R=X_BESS_GEO-X_BESS_T;   Y_BESS_R=Y_BESS_GEO-Y_BESS_T;   Z_BESS_R=Z_BESS_GEO-Z_BESS_T;

A=[];X=[];
for n=1:size(IB_ETRS,1)
    Ai=[0,-Z_ETRS_R(n),Y_ETRS_R(n),1,0,0,X_ETRS_R(n)
        Z_ETRS_R(n),0,-X_ETRS_R(n),0,1,0,Y_ETRS_R(n)
        -Y_ETRS_R(n),X_ETRS_R(n),0,0,0,1,Z_ETRS_R(n)];
    A=[A;Ai];
    X=[X;X_BESS_R(n);Y_BESS_R(n);Z_BESS_R(n)];
end

h=(A'*A)^(-1)*A'*X;
q=h(end);
R=[1,h(3),-h(2) 
  -h(3),1,h(1)    
   h(2),-h(1),1];

v=A*h-X;
vx=sum(v(1:3:end)); vy=sum(v(2:3:end)); vz=sum(v(3:3:end));

X0=[X_BESS_T;Y_BESS_T;Z_BESS_T]-q*R*[X_ETRS_T;Y_ETRS_T;Z_ETRS_T];
%% d - transformace
B_PB=(mer1(:,2)+mer1(:,3)/60+mer1(:,4)/3600)*RAD;
L_PB=(mer1(:,5)+mer1(:,6)/60+mer1(:,7)/3600)*RAD;
H_PB=mer1(:,8);

W_ETRS1=sqrt(1-e_ETRS*(sin(B_PB).^2));
N_ETRS1=a_ETRS./W_ETRS1;

X_ETRS_GEO_PB=cos(B_PB).*cos(L_PB).*(N_ETRS1+H_PB);
Y_ETRS_GEO_PB=cos(B_PB).*sin(L_PB).*(N_ETRS1+H_PB);
Z_ETRS_GEO_PB=(N_ETRS1*(1-e_ETRS)+H_PB).*sin(B_PB);
P_ETRS_GEO_PB=[X_ETRS_GEO_PB';Y_ETRS_GEO_PB';Z_ETRS_GEO_PB']';
P_ETRS_GEO_B=[P_ETRS_GEO';P_ETRS_GEO_PB];

B_BESS=(X0+q*R*P_ETRS_GEO_B')';

%% e převod z geocentrických souřadnic na zeměpisné
B_BESS_X=B_BESS(:,1); B_BESS_Y=B_BESS(:,2); B_BESS_Z=B_BESS(:,3);

L_PB=atan2(B_BESS_Y,B_BESS_X);
B_PB=atan(B_BESS_Z./sqrt(B_BESS_X.^2+B_BESS_Y.^2)).*(1/(1-e_BESSEL));

for n=1:10
    N_BESS=a_BESSEL./sqrt(1-e_BESSEL.*(sin(B_PB).^2));
    H_BESS=(sqrt(B_BESS_X.^2+B_BESS_Y.^2)./cos(B_PB))-N_BESS;
    B_PB=atan((B_BESS_Z./sqrt(B_BESS_X.^2+B_BESS_Y.^2)).*(1-(N_BESS.*e_BESSEL)./(N_BESS+H_BESS)).^(-1));
end
%% f - Převod zeměpiných souřadnic do S-JTSK
U_PB=2*atan((1/k)*(((1-e_bessel*sin(B_PB))./(1+e_bessel*sin(B_PB))).^(alfa*e_bessel*(1/2))).*(tan(B_PB/2+pi/4).^alfa))-pi/2;
DV_PB=alfa*(((24+50/60)*RAD)-L_PB);

S_PB=asin(sin(UQ)*sin(U_PB)+cos(UQ)*cos(U_PB).*cos(DV_PB));
D_PB=asin((sin(DV_PB).*cos(U_PB))./cos(S_PB));

ro_PB=r0.*(tan(S0/2+pi/4)./tan(S_PB/2+pi/4)).^n1;
eps_PB=n1*D_PB;

X_PB=ro_PB.*cos(eps_PB);
Y_PB=ro_PB.*sin(eps_PB);

%% 
PB_kontrola=[IB_S_JTSK;mer_kontrola]; PB=[PB_kontrola(:,1),Y_PB,X_PB,H_BESS];
Vymera_katastr=10076+2705;
y1=PB(7:end,2);x1=PB(7:end,3);
Vym_podrobne_body=polyarea(y1,x1);

IBY=PB_kontrola(1:6,2)-PB(1:6,2); IBX=PB_kontrola(1:6,3)-PB(1:6,3);
p_odch=sqrt((PB_kontrola(:,2)-PB(:,2)).^2+(PB_kontrola(:,3)-PB(:,3)).^2);
sig_v=sqrt((IBY'*IBY+IBX'*IBX)/size(IB_ETRS,1));
delta=PB_kontrola-PB;

T_x=mean(x1);T_y=mean(y1); 
X_min=min(x1);X_max=max(x1);Y_min=min(y1); Y_max=max(y1);
dy_1=delta(7:end,2);dx_1=delta(7:end,3);

figure (1)
plot([-y1;-y1(1)],[-x1;-x1(1)]);
hold on
plot(-y1,-x1,'LineStyle','none','Marker','x','Color','k','MarkerSize',10);
quiver(-y1,-x1,dy_1,dx_1,0.2,'r')
legend('Hranice parcely','Lomové body','Vektory odchylek zvětšna 20x','Location','southeast')
xlabel('Y[m]')
ylabel('X[m]')
text(-T_y,-T_x+20,'Parcela č. 656');
ylim([-X_max-10,-X_min+10]);
xlim([-Y_max-10,-Y_min+10]);
title('Grafické znázornění transformovaných bodů a vektorových odchylek')
text(-y1+2,-x1+2,num2str(PB(7:end,1)));
axis equal
hold off
