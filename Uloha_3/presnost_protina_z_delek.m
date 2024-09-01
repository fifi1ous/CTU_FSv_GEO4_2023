function [P1,EX1_kov] = presnost_protina_z_delek(pd,B,SS)
%Protínání a jeho přesnosti protínání z délek
%Input:         pd  -   přesnost měření délek
%               B   -   souřadnice určovaných bodů
%               SS  -   souřadnice měřické sítě
%
%Output:        P1  -   matice všech kombinací bez vlivu podkladu   1.sl měření délky
%                                                                   2.sl měření délky
%               EX1_kov - kovariančí matice (složená matice 2x2, pořadí
%                                            řádek      -   měření délky z
%                                                           bodu 1.sl v matici P1
%                                            sloupeček  -   měření délky z
%                                                           bodu 2.sl v matici P1)
%
%% Výpočet
RAD=pi/200;
Ex1_kov=[];EX1_kov=[];
UP=[0,0;0,0];
for n=1:size(SS,1)
    for m=1:size(SS,1)
        if n~=m
        [sm, Spa] = cart2pol(B(1)-SS(n,2),B(2)-SS(n,3));
        [sm, Spb] = cart2pol(B(1)-SS(m,2),B(2)-SS(m,3));
        A1=[-((SS(n,3)-B(2))/Spa),-((SS(n,2)-B(1))/Spa)
            -((SS(m,3)-B(2))/Spb),-((SS(m,2)-B(1))/Spb)];
        el=[pd,pd];
        EL=diag(el.^2);
        K=A1^(-1);
        Ex1=K*EL*K';                                    
        [a11,b11,alf11] =par_el_chyb(Ex1);
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