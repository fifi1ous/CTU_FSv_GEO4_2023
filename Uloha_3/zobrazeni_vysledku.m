function [X,Y,SSA,SSB,SSC] = zobrazeni_vysledku(P,SS,B)
%Protínání a jeho přesnosti protínání z délek
%Input:         P   -   Výsledná matice
%               SS  -   Souřadnice měřické sítě
%               B   -   Určoavný bod
%
%Output:        X   -   Souřadnice X elipsy chyb
%               Y   -   Souřadnice Y elipsy chyb
%               SSA -   Souřadnice
%               SSB -   Souřadnice
%               SSC -   Souřadnice
%
%%
RAD=pi/200;
t=linspace(0,2*pi,1000);
alfa=P(:,end-2)*RAD; a=P(:,end-1); b=P(:,end);

X=B(2)+a.*cos(t).*cos(alfa)-b.*sin(t).*sin(alfa);
Y=B(1)+b.*sin(t).*cos(alfa)+a.*cos(t).*sin(alfa);
SS=sortrows(SS,1);

if  rem(size(P,2),2)==0
    P2=P(1,1:2);
    PRV=sortrows(P2(1,:)',1);
    SS2=SS(PRV,:);
    SSA=[SS2(1,2:3);B];
    SSB=[SS2(2,2:3);B];
    SSC=[SS2(1,2:3);SS2(2,2:3)];
else
    P2=P(1,1:3);
    PRV=sortrows(P2(1,:)',1);
    SS2=SS(PRV,:);
    SSA=[SS2(1,2:3);B];
    SSB=[SS2(2,2:3);B];
    SSC=[SS2(3,2:3);B];
end
end