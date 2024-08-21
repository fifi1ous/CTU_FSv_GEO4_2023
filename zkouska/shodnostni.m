clc; clear; format long G
%  cb    Y        X       
IB=[501 499.997   500.003 
    502 499.997  1000.003 
    503 1000.003  999.997 
    504 1000.004  500.004 ];

ib=[501 357.701   626.169
    502 134.094  1073.383
    502 581.448  1296.849
    504 804.914   849.636];

ss=[1   100      500
    2   125      400
    3   400      600
    4   200      300];
%% nastavení
IB=sortrows(IB); ib=sortrows(ib);
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
w=atan2(xr'*Yr-yr'*Xr,xr'*Xr+yr'*Yr);

X0=XT-xT*cos(w)+yT*sin(w);
Y0=YT-yT*cos(w)-xT*sin(w);
%% Transformační rovnice
X=X0+SS(:,2)*cos(w)-SS(:,1)*sin(w);
Y=Y0+SS(:,1)*cos(w)+SS(:,2)*sin(w);

SEZ=[cb,Y,X];
%% rozdíly souřadnic
wxi=xr*cos(w)-yr*sin(w)-(IB(:,3)-XT);
wyi=yr*cos(w)+xr*sin(w)-(IB(:,2)-YT);
Ewxi=sum(wxi)
Ewyi=sum(wyi)
%% kritéria přesnosti
mira_identity=sqrt((wxi'*wxi+wyi'*wyi)/n)
smerodatna_odchylka_danych_bodu_klice=sqrt((wxi'*wxi+wyi'*wyi)/(2*n-3))