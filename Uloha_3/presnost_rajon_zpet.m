function [P1,EX1_kov] = presnost_rajon_zpet(pu,pd,B,SS)
%Protínání a jeho přesnosti rajón zpět
%Input:         pu  -   přesnost měření úhlů
%               pd  -   přesnost měření délek
%               B   -   souřadnice určovaných bodů
%               SS  -   souřadnice měřické sítě
%
%Output:        P1  -   matice všech kombinací bez vlivu podkladu   1.sl měření délky i úhlu
%                                                                   2.sl měření pouze úhlu
%               EX1_kov - kovariančí matice (složená matice 2x2, pořadí
%                                            řádek      -   z bodu 1.sl v matici P1, měření délky i úhlu
%                                            sloupeček  -   z bodu 2.sl v matici P1, měření pouze úhlu   )
%
%% Výpočet
RAD=pi/200;
D=[-1,1,0;0,0,1];
Ex1_kov=[];EX1_kov=[];
UP=[0,0;0,0];
for n=1:size(SS,1)
    for m=1:size(SS,1)
        if n~=m
        [sm, Spa] = cart2pol(B(1)-SS(n,2),B(2)-SS(n,3));
        [sm, Spb] = cart2pol(B(1)-SS(m,2),B(2)-SS(m,3));
        %Matice A1 měření
        A1=[(SS(m,2)-B(1))/(Spb^2)-(SS(n,2)-B(1))/(Spa^2),-(SS(m,3)-B(2))/(Spb^2)+(SS(n,3)-B(2))/(Spa^2)
            -(SS(n,3)-B(2))/(Spa),-(SS(n,2)-B(1))/(Spa)];

        el=[pu;pu;pd];                                   %Přesnost měření
        K=A1^(-1)*D;
        EL=diag(el.^2);                         
        Ex1=K*EL*K';                                     %Kovarianční matice souřadnic bez podkladu
        [a11,b11,alf11] =par_el_chyb(Ex1);               %Parametry elipsy chyb
        a1(m,n)=a11; b1(m,n)=b11; alf1(m,n)=alf11;
        sigx1(m,n)=sqrt(Ex1(1));sigy1(m,n)=sqrt(Ex1(4));sigxy1(m,n)=sqrt((Ex1(1)+Ex1(4))/2);
        KAR(m,n)=m;
        Ex1_kov=[Ex1_kov,Ex1];
        else
            Ex1_kov=[Ex1_kov,UP];
        end
        KAR1(m,n)=n;
    end
    EX1_kov=[EX1_kov;Ex1_kov];
    Ex1_kov=[];
end

KAR1=reshape(KAR1,size(SS,1)*size(SS,1),1);
KAR=reshape(KAR,size(SS,1)*size(SS,1),1);
a1=reshape(a1,size(SS,1)*size(SS,1),1); b1=reshape(b1,size(SS,1)*size(SS,1),1); alf1=reshape(alf1,size(SS,1)*size(SS,1),1);
sigx1=reshape(sigx1,size(SS,1)*size(SS,1),1);sigy1=reshape(sigy1,size(SS,1)*size(SS,1),1);sigxy1=reshape(sigxy1,size(SS,1)*size(SS,1),1);

P1=[KAR1,KAR,sigx1,sigy1,sigxy1,alf1/RAD,a1,b1];
[r,s]=find(P1(:,2)==0);
P1(r,:)=[];
end