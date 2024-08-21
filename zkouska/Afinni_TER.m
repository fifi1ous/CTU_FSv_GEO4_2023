clc; clear; format long g
% Tereza Cernohousova, 27.05. 2023
%zadani: a = 17, b = 1

%% PRACE S MATICEMI
a = 17; b = 1;
A = [a 2 3
     2 4 5
     3 5 b];
stopa = trace(A);
determinant = det(A);
inverzni = inv(A);
vlc = eig(A);
fprintf('Stopa matice A = %.0f\n', stopa)
fprintf('Determinant matice A = %.3f\n', determinant)
fprintf('Inverzni matice A = \n')
fprintf('  %11.6f %11.6f %11.6f \n', inverzni')
fprintf('Vlastní čísla matice A = %.5f\n %32.5f\n %32.5f\n', vlc)

%% NACTENI DAT NA TRANSFORMACE
fid = fopen('ss.txt','r');
id = fscanf(fid,'%f', [5 3])';
body = fscanf(fid,'%f', [3 inf])';

%% HELMERTOVA TRANSFORMACE (podobnostni) - podle U1
% REDUKCE SOURADNIC K TEZISTI
% MISTNI
xt = mean(id(:,end));
yt = mean(id(:,end-1));
% SJTSK
Xt = mean(id(:,3));
Yt = mean(id(:,2));
% REDUKOVANE SOURADNICE
xr = id(:,end) - xt;
yr = id(:,end-1) - yt;
Xr = id(:,3) - Xt;
Yr = id(:,2) - Yt;
% TRANSFORMACNI KLIC
lam1 = sum(xr.*Xr + yr.*Yr)/sum(xr.^2 + yr.^2);
lam2 = sum(xr.*Yr - yr.*Xr)/sum(xr.^2 + yr.^2);
X0 = Xt - lam1*xt + lam2*yt;
Y0 = Yt - lam1*yt -lam2*xt;
q = sqrt(lam1^2 + lam2^2);
omega = atan2(lam2,lam1);
%TRANSFORMACNI ROVNICE (identicke a podrobne body)
% 1:3 jsou identicke body, 4:end jsou podrobne body 1-6
X_Hel = X0 + lam1.*[id(:,end); body(:,3)] - lam2.*[id(:,end-1); body(:,2)];
Y_Hel = Y0 + lam1.*[id(:,end-1); body(:,2)] + lam2.*[id(:,end); body(:,3)];
% ROZDILY SOURADNIC (opravy, plati pro identicke body)
v_Xi = lam1.*xr - lam2.*yr - (id(:,3) - Xt);
v_Yi = lam1.*yr + lam2.*xr - (id(:,2) - Yt);
% kontrola musi byt nulova nebo temer nula
kontrola_vx = sum(v_Xi);
kontrola_vy = sum(v_Yi);
% vx = X_(1:3) - id(1:3,3)
% [Y1,X1,q1,omega1,lam11,lam22] = helmertova_transformace(id(:,2:3),[id(:,end-1:end);body(:,2:3)])

%% AFINNI TRANSFORMACE
xa = id(1,end); ya = id(1,end-1);
xb = id(2,end); yb = id(2,end-1);
xc = id(3,end); yc = id(3,end-1);
XA = id(1,3); YA = id(1,2);
XB = id(2,3); YB = id(2,2);
XC = id(3,3); YC = id(3,2);
% DETERMINANT SOUSTAVY
% D = xa*(yb-yc)+xb*(yc-ya)+xc*(ya-yb)
AA = [xa, ya, 1
      xb, yb, 1
      xc, yc, 1];
D = det(AA);
%KONSTANTY a1,b1,a2,b2
a1 = ((yc-ya)*(XB-XA)-(yb-ya)*(XC-XA))/D;
b1 = ((xc-xa)*(XB-XA)-(xb-xa)*(XC-XA))/D;
a2 = ((xb-xa)*(YC-YA)-(xc-xa)*(YB-YA))/D;
b2 = ((yc-ya)*(YB-YA)-(yb-ya)*(YC-YA))/D;
% ↑ tohle jsou ty jedine parametry, ktere musim spocitat ze strany 17.
% I kdyz pise, ze nemusime pocitat X0 a Y0, tak i tak je musime spocitat.
% Vyjadrime je ze vzorce 1.9 na strane 16
X0_501 = id(1,3) - a1*id(1,end) + b1*id(1,end-1);
Y0_501 = id(1,2) - b2*id(1,end) - a2*id(1,end-1);

X0_502 = id(2,3) - a1*id(2,end) + b1*id(2,end-1);
Y0_502 = id(2,2) - b2*id(2,end) - a2*id(2,end-1);

X0_503 = id(3,3) - a1*id(3,end) + b1*id(3,end-1);
Y0_503 = id(3,2) - b2*id(3,end) - a2*id(3,end-1);
% ↑ tohle potom smazu
% ↑ Tady jsem si overila, ze X0 a Y0 jsou pro vsechny 3 identicky body
% stejne, takze jsou spocitany dobre a ja muzu spocitat identicke body

X0 = id(1,3) - a1*id(1,end) + b1*id(1,end-1);
Y0 = id(1,2) - b2*id(1,end) - a2*id(1,end-1);

% TRANSFORMACE PODROBNYCH BODU
X_afin = X0 + a1.*body(:,end) - b1.*body(:,end-1);
Y_afin = Y0 + b2.*body(:,end) + a2.*body(:,end-1);

% pro kontrolu ted i identicke body
X_afin_id = X0 + a1.*id(:,end) - b1.*id(:,end-1);
Y_afin_id = Y0 + b2.*id(:,end) + a2.*id(:,end-1);

%% VYPOCET VYMERY
vymera_Hel = polyarea(X_Hel(4:end), Y_Hel(4:end))
vymera_afin = polyarea(X_afin, Y_afin)

%%
fclose(fid);