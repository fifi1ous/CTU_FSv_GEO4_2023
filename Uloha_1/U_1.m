clc; clear; format long G
%% Úloha 0
RAD=pi/200;
m=0.99987;

nazev='mer.txt';
[mereni,s,fi] = ul_0(nazev,m);
mereni = sortrows(mereni,1);
%% Import souřadnic
fid=fopen('BODY.txt','r');
IB=fscanf(fid,'%d %f %f %d %d',[5,inf])';
fclose(fid);
IB1=sortrows(IB,5);
%%
KAT=IB1(:,1);
IB = [IB1(:,end),IB1(:,2:3)];

r1=size(IB1,1);
r=size(IB,1);
%% Výpočet souřadnic
[x,y]=pol2cart(mereni(:,3)*RAD,s);
BODY=[mereni(:,1),y+1000,x+5000];
%% Transformace Shodnostní
ss=BODY(1:r,:);
ss(8:9,:)=[];
ss(2:3,:)=[];
IB(8:9,:)=[];
IB(2:3,:)=[];
r=size(IB,1);

Yt=mean(IB(:,2)); Xt=mean(IB(:,3));
yt=mean(ss(:,2)); xt=mean(ss(:,3));

Yr=IB(:,2)-Yt; Xr=IB(:,3)-Xt;
yr=ss(:,2)-yt; xr=ss(:,3)-xt;

%transformační klíč
w=atan2((xr'*Yr-yr'*Xr),(xr'*Xr+yr'*Yr));
X0=Xt-xt*cos(w)+yt*sin(w);
Y0=Yt-yt*cos(w)-xt*sin(w);

%transformační rovnice
XV=X0+BODY(:,3)*cos(w)-BODY(:,2)*sin(w);
YV=Y0+BODY(:,2)*cos(w)+BODY(:,3)*sin(w);
XVi=XV; XVi(8:9,:)=[];  XVi(2:3,:)=[];
YVi=YV; YVi(8:9,:)=[];  YVi(2:3,:)=[];
BODY1=[BODY(:,1),YV,XV];

%kontrola
wx1=xr*cos(w)-yr*sin(w)-(IB(:,3)-Xt);
wy1=yr*cos(w)+xr*sin(w)-(IB(:,2)-Yt);

if round(sum(wy1),6)~=0 || round(sum(wx1),6)~=0
    error('Kontrola není blízká 0')
end

sigmav_shod=((wx1'*wx1+wy1'*wy1)/r)^(1/2);
sigmao_shod=((wx1'*wx1+wy1'*wy1)/(2*r-3))^(1/2);
%% Transformace Podobnostní
lam1=(xr'*Xr+yr'*Yr)/(xr'*xr+yr'*yr);
lam2=(xr'*Yr-yr'*Xr)/(xr'*xr+yr'*yr);

XO=Xt-lam1*xt+lam2*yt;
YO=Yt-lam1*yt-lam2*xt;

q=sqrt(lam1^2+lam2^2);

w1=atan2(lam2,lam1);
if w1~=w
    error('Omega není stejná')
end

XW=XO+lam1*ss(:,3)-lam2*ss(:,2);
YW=YO+lam1*ss(:,2)+lam2*ss(:,3);

vx1=lam1*xr-lam2*yr-(IB(:,3)-Xt);
vy1=lam1*yr+lam2*xr-(IB(:,2)-Yt);

if round(sum(vy1),6)~=0 || round(sum(vx1),6)~=0
    error('Kontrola není blízká 0')
end

BODY3=[-IB(:,2:3),-YVi(1:r,:),-XVi(1:r,:),-YW,-XW];
%% Identifikace chybného bodu
r_2=xr.^2+yr.^2;
f=(1-(1/r+r_2/sum(r_2))).^(-1);
Dvv=f.*(vx1.^2+vy1.^2);

sigmav_pod=((vx1'*vx1+vy1'*vy1)/r)^(1/2);
sigmao_pod=((vx1'*vx1+vy1'*vy1)/(2*r-4))^(1/2);
%% Grafický výstup
[SST,SSK] = ne_vzd_kruznice(IB);

XV(8:9)=[]; XV(2:3)=[];
YV(8:9)=[]; YV(2:3)=[];
figure(10)
plot(-YW,-XW,"Marker","x","LineStyle","none","Color",'blue')
hold on
plot(-YV,-XV,"Marker",".","LineStyle","none","Color",'red')
plot(-IB(:,2),-IB(:,3),"Marker","o","LineStyle","none","Color",'g')
plot(-SSK(:,1),-SSK(:,2),'Color','k')
legend('Body spočítané podobnostní transformací','Body spočítané shodnostní transformací','Identické body','Kružnice','Location','southeast')
hold off

plochaMistni=polyarea(BODY(r1+1:end,3),BODY(r1+1:end,2));
plochaSJ=polyarea(BODY1(r1+1:end,3),BODY1(r1+1:end,2));

for n=1:r
    A=[BODY3(n,1:2);BODY3(n,3:4);BODY3(n,5:6)];
    [SST1,SSK1] = ne_vzd_kruznice(A);
    [del,odch,sor] = delky(BODY3(n,:));
    figure(n)
    plot(BODY3(n,1),BODY3(n,2),"Marker","o","LineStyle","none","Color",'g');
    hold on
    plot(BODY3(n,3),BODY3(n,4),"Marker",".","LineStyle","none","Color",'red');
    plot(BODY3(n,5),BODY3(n,6),"Marker","x","LineStyle","none","Color",'blue');
    plot(sor(:,1),sor(:,2),"Color",'k','LineStyle',':')
    plot(SSK1(:,1),SSK1(:,2),"LineStyle","none");
    TEXT=[num2str(odch(1)),' mm'];
    xlabel(TEXT)
    TEXT=[num2str(odch(2)),' mm'];
    ylabel(TEXT)
    TEXT=['Bod č. ',num2str(IB(n))];
    title (TEXT)
    legend('Identický bod','IB spočítaný shodnostní transformací','IB spočítaný podobnostní transformací','Odchylka','Location','northeast')
    hold off
end

%% funkce
function [mereni,dh,fi] = ul_0(nazev,m)

RAD=pi/200;
fid = fopen(nazev,'r');
mereni=fscanf(fid,'%d %f %f %f',[4,inf])';
fclose(fid);

fi = 0.00998*mereni(:,2)/1000.*sin(mereni(:,4)*RAD)*RAD;
dh=mereni(:,2).*sin(mereni(:,4)*RAD-fi)*m;
end

function [SST,SSK] = ne_vzd_kruznice(ss)
c=size(ss,1);
d=size(ss,2);
D=0;
if d~=2
    D=1;
end
max_vzd=0;
for n=1:c
    vzd1=sqrt((ss(n,1+D)-ss(:,1+D)).^2+(ss(n,2+D)-ss(:,2+D)).^2);
    vzd=max(vzd1);
    r1=find(vzd1==vzd);
    if vzd>max_vzd
        max_vzd=vzd;
        m=r1;
        n1=n;
    end
end
max_vzd1=max_vzd/4*3;
SST=[(ss(m,1+D)+ss(n1,1+D))/2,(ss(m,2+D)+ss(n1,2+D))/2];
th=0:pi/1000:2*pi;
SSK(:,1)=max_vzd1*sin(th)+SST(1);
SSK(:,2)=max_vzd1*cos(th)+SST(2);
end

function [del,odch,sor] = delky(BODY)
for n=1:size(BODY,1)
    A=sqrt((BODY(n,1)-BODY(n,3))^2+(BODY(n,2)-BODY(n,4))^2);
    B=sqrt((BODY(n,1)-BODY(n,5))^2+(BODY(n,2)-BODY(n,6))^2);
    if A>B
        del(n)=A;
        Y(n)=sqrt((BODY(n,1)-BODY(n,3))^2);
        X(n)=sqrt((BODY(n,2)-BODY(n,4))^2);
        odch=[Y,X];
        sor=[BODY(n,5),BODY(n,6);BODY(n,1),BODY(n,2);BODY(n,3),BODY(n,4);];
    else
        del(n)=B;
        Y(n)=sqrt((BODY(n,1)-BODY(n,5))^2);
        X(n)=sqrt((BODY(n,2)-BODY(n,6))^2);
        odch=[Y,X];
        sor=[BODY(n,3),BODY(n,4);BODY(n,1),BODY(n,2);BODY(n,5),BODY(n,6);];
    end   
    odch=odch*100*10;
    odch=round(odch);
end
end



