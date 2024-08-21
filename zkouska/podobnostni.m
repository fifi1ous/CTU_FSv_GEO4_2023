clc; clear; format long G
%% zadání hodnot
%  cb    Y        X        Vynechání bodu
IB=[1 499.997   500.003      1
    2 499.997  1000.003      1
    3 1000.003  999.997      1
    4 1000.904  500.004      0   ];

ib=[1 357.701   626.169
    2 134.094  1073.383
    2 581.448  1296.849
    4 804.914   849.636];

ss=[1   100      500
    2   125      400
    3   400      600
    4   200      300];
%% nastavení
IB=sortrows(IB); ib=sortrows(ib);
[i]=find(IB(:,4)==1);
IB=IB(i,1:3);   ib=ib(i,:);
n=size(IB,1);
cb=[ib(:,1);ss(:,1)];
cb_ib=ib(:,1);
SS=[ib(:,2:3);ss(:,2:3)];
%% Souřadnice těžiště
XT=mean(IB(:,3));   YT=mean(IB(:,2));
xT=mean(ib(:,3));   yT=mean(ib(:,2));
%% Redukce k těžišti
Xr=IB(:,3)-XT;  Yr=IB(:,2)-YT;
xr=ib(:,3)-xT;  yr=ib(:,2)-yT;
%% transformační klíč
lam1=(xr'*Xr+yr'*Yr)/(xr'*xr+yr'*yr);
lam2=(xr'*Yr-yr'*Xr)/(xr'*xr+yr'*yr);

q=sqrt(lam1^2+lam2^2);
w=atan2(lam2,lam1);

X0=XT-lam1*xT+lam2*yT;
Y0=YT-lam1*yT-lam2*xT;
%% transformace souřadnic
X=X0+lam1*SS(:,2)-lam2*SS(:,1);
Y=Y0+lam1*SS(:,1)+lam2*SS(:,2);
SEZ=[cb,Y,X];
%% kontrola
vxi=lam1*xr-lam2*yr-(IB(:,3)-XT);   vyi=lam1*yr+lam2*xr-(IB(:,2)-YT);
Evxi=sum(vxi)  
Evyi=sum(vyi)
%% identifikace chybného bodu
r=yr.^2+xr.^2;
f=(1-((1/n)+(r/sum(r)))).^(-1);
delvv=f.*(vxi.^2+vyi.^2);
delvv=[cb_ib,delvv]
%% odchylky
mira_identity=sqrt((vxi'*vxi+vyi'*vyi)/n)
smerodatna_odchylka_danych_bodu_klice=sqrt((vxi'*vxi+vyi'*vyi)/(2*n-4))