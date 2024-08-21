clc; clear; format long G
%nastavení 
global RAD

RAD=pi/200;
%  cb:     Y           X       Byl měřen směr:    Byla měřena délka:
SS=[1  877363.098  1035205.687       1                     1 
   19  877557.207  1034622.298       1                     1 
   29  877160.751  1034593.248       1                     1 
   35  876977.558  1034767.056       1                     1
   48  877238.727  1035213.431       1                     1 ];
 
%  cb:     Y           X         Orientační posun:
ST=[600 877321.655  1034878.252    390.0021];

% byl měřen podrobný bod ANO=1 NE=0
byl_podrobny_bod=1;


if byl_podrobny_bod==1

      %  cb:     Y           X         
    PB= [514 877291.742  1034880.401];

else
    PB=[];
end
if sum(SS(:,4))==0
   ST=[ST(1),ST(2),ST(3)];
end


%      cb:    směr:  délka:   přesnost měření směrů:    přesnost měření délek:
mereni=[1    18.0142 330.057   ((90/1000)*RAD)/18.0142      5/1000
        19  162.6415 347.833   (1/1000)*RAD                 2/1000
        29  242.7221 327.283   (1/1000)*RAD                 5/1000
        35  290.1032 361.626   (1/1000)*RAD                 2/1000
        48  394.5568 345.297   (1/1000)*RAD                 2/1000    ];

merP=[ 514  314.5650  30.001   (1/1000)*RAD                 2/1000    ];

sm0=1;

%% nastavení
SS=sortrows(SS,1);
mereni=sortrows(mereni,1);
[i_sm]=find(SS(:,4)==1);
[i_d]=find(SS(:,5)==1);
SS_smer=SS(i_sm,2:3);
SS_del=SS(i_d,2:3);
X=[ST(:,2:end),PB(:,2:end)];     %X(3)=X(3)*RAD;         %neznámé jsou seřazeny Yst,Xst,OPst,YPB,XPB
if sum(SS(:,4))~=0
   X(3)=X(3)*RAD;
end
mereni_d=mereni(i_d,3);
mereni(:,2)=mereni(:,2)*RAD;
mereni_sm=mereni(i_sm,2);
if byl_podrobny_bod==1
p=[mereni(i_d,5);merP(5);mereni(i_sm,4);merP(4)];                       %sestavování jsou délky, potom směry
else
p=[mereni(i_d,5);mereni(i_sm,4)];   
end
mer_d=length(mereni_d); mer_sm=length(mereni_sm);
nez=length(X);
if byl_podrobny_bod==1
    l=[mereni_d;merP(3);mereni_sm;merP(2)*RAD];
else
    l=[mereni_d;mereni_sm];
end
%% Matice vah                               
P=eye(size(p,1)).*((sm0^2)./(p.^2)); 

%% Vyrovnani 
for kv=1:10                     %Počet chtěných iterací
% Matice A
A1Y=1./(((X(1) - SS_smer(:,1)).^2./(X(2) - SS_smer(:,2)).^2 + 1).*(X(2) - SS_smer(:,2)));
A1X=-(X(1) - SS_smer(:,1))./(((X(1) - SS_smer(:,1)).^2./(X(2) - SS_smer(:,2)).^2 + 1).*(X(2) - SS_smer(:,2)).^2);
AORSM=[A1Y,A1X,ones(mer_sm,1),zeros(mer_sm,2)];

if byl_podrobny_bod==1
A1YPB=1./(((X(1) - X(4)).^2./(X(2) - X(5)).^2 + 1).*(X(2) - X(5)));
A1XPB=-(X(1) - X(4))./(((X(1) - X(4)).^2./(X(2) - X(5)).^2 + 1).*(X(2) - X(5)).^2);
AYsmPB=-1./(((X(1) - X(4)).^2./(X(2) - X(5)).^2 + 1).*(X(2) - X(5)));
AXsmPB=(X(1) - X(4))./(((X(1) - X(4)).^2./(X(2) - X(5)).^2 + 1).*(X(2) - X(5)).^2);
A0PBSM=[A1YPB,A1XPB,1,AYsmPB,AXsmPB];
else
A0PBSM=[];
end

A2Y=(2*X(1) - 2* SS_del(:,1))./(2*sqrt((X(2) -  SS_del(:,2)).^2 + (X(1) -  SS_del(:,1)).^2));
A2X=(2*X(2) - 2* SS_del(:,2))./(2*sqrt((X(2) -  SS_del(:,2)).^2 + (X(1) -  SS_del(:,1)).^2));
AORD=[A2Y,A2X,zeros(mer_d,3)];

if byl_podrobny_bod==1
A2YPB=(2*X(1) - 2* X(4))./(2*sqrt((X(2) -  X(5)).^2 + (X(1) -  X(4)).^2));
A2XPB=(2*X(2) - 2* X(5))./(2*sqrt((X(2) -  X(5)).^2 + (X(1) -  X(4)).^2));
AYdPB=-(2*X(1) - 2*X(4))./(2*sqrt((X(2) -  X(5)).^2 + (X(1) -  X(4)).^2));
AXdPB=-(2*X(2) - 2*X(5))./(2*sqrt((X(2) -  X(5)).^2 + (X(1) -  X(4)).^2));
A0PBD=[A2YPB,A2XPB,0,AYdPB,AXdPB];
else
A0PBD=[];
end

A=[AORD;A0PBD;AORSM;A0PBSM];
if byl_podrobny_bod==0
    A(:,4:5)=[];
end
if sum(SS(:,4))==0
   A(:,3)=[];
end

%matice l
sml=atan2(SS_smer(:,1)-X(1),SS_smer(:,2)-X(2));
for i=1:length(sml)
    if sml(i)<0
        sml(i)=sml(i)+2*pi;
    end
    sml(i)=sml(i)-X(3);
    if sml(i)<0
       sml(i)=sml(i)+2*pi;
    end
end

if byl_podrobny_bod==1
    smlPB=atan2(X(4)-X(1),X(5)-X(2));
    if smlPB<0
        smlPB=smlPB+2*pi;
    end
    smlPB=smlPB-X(3);
    if smlPB<0
       smlPB=smlPB+2*pi;
    end
else
    smlPB=[];
end

dl=sqrt((SS_del(:,1)-X(1)).^2 + (SS_del(:,2)-X(2)).^2);

if byl_podrobny_bod==1
    dlPB=sqrt((X(4)-X(1)).^2 + (X(5)-X(2)).^2);
else
    dlPB=[];
end

l0=[dl;dlPB;sml;smlPB];

lc=l0-l;

N=A'*P*A;
dx=-N^(-1)*A'*P*lc;
v=A*dx+lc;
X=X+dx';
end
%% Kejdy okolo co prostě skori potřebuje, aby byl happy jak dva Grepy

Sig0=sqrt((v'*P*v)/(mer_d+mer_sm-nez));                    %empirická hodnota jednotkové variance

Qx=N^(-1);                  % matice váhových koeficientů neznámých
Qv=P^(-1)-A*(N^(-1))*A';    % matice váhových koeficietnů oprav měření
QLL=A*(N^(-1))*A';          % matice váhových koeficientl vyrovaných měření

Ex=sm0^2*(N^(-1));                  % kovarianční matice neznámých
Ev=sm0^2*(P^(-1)-A*(N^(-1))*A');    % kovarianční matice oprav
ELL=sm0^2*A*(N^(-1))*A';            % kovarianční matice vyrovnaných měření

%Správnost odhadu přesnosti měření
alf=0.05;
kvantil=chi2inv(1-alf,(mer_d+mer_sm-nez));
kvantil=sqrt(kvantil/(mer_d+mer_sm-nez));
nas_odhad=Sig0/sm0;                 % kvantil by měl být větší než nas_odhad, pokud ano tak Odhad přesnosti měření je adekvátní reálné skutečnosti 


if sum(SS(:,4))~=0
   smsmernik=sqrt(Ex(3,3));
end
smYST=Ex(1,1);smXST=Ex(2,2);
if byl_podrobny_bod==1
smYPB=Ex(4,4);smXPB=Ex(5,5);
end

smerodatna_odchylky_stanoviska=sqrt((smYST+smXST)/2)
if byl_podrobny_bod==1
smerodatna_odchylky_podrobneho_bodu=sqrt((smYPB+smXPB)/2)
end
%% elipsa chyb
cST=sqrt((smXST-smYST)^2+4*(Ex(1,2)^2));
aST=sqrt((smXST+smYST+cST)/2)
bST=sqrt((smXST+smYST-cST)/2)
alfaST=(atan2(2*Ex(1,2),smXST-smYST)/2)/RAD;
if alfaST<0
    alfaST=alfaST+400;
end
alfaST

if byl_podrobny_bod==1
cPB=sqrt((smXPB-smYPB)^2+4*(Ex(4,5)^2));
aPB=sqrt((smXPB+smYPB+cPB)/2)
bPB=sqrt((smXPB+smYPB-cPB)/2)
alfaPB=(atan2(2*Ex(4,5),smXPB-smYPB)/2)/RAD;
if alfaPB<0
    alfaPB=alfaPB+400;
end
alfaPB
end
%% Skořiho analýza oprav protože co kdyby to tam dal ten ....
if sum(SS(:,4))~=0
kov_v = 2.1^2*(inv(P) - A*inv(A'*P*A)*A'); % kovariační matice oprav
kov_v = 0.5 * (kov_v+kov_v'); % zaručení symetrie matice, v matici S jsou pak reálná čísla 
[U,S] = eig(kov_v); % rozklad matice
for i = 1:4
S(i,i) = 0; % na místě praktických nulových prvků se vloží přesně nula
end
for i = 5:7
S(i,i) = 1 / sqrt(S(i,i));
end
kov_v_odm = U*S*U';
v_hom = kov_v_odm*v;                % homogenizované rovnice oprav, měli by být menší než 1.96, jinak přesahují interval spolehlivosti
end