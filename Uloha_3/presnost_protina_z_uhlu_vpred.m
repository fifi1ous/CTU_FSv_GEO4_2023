function [P1,EX1_kov] = presnost_protina_z_uhlu_vpred(pu,B,SS)
%Protínání a jeho přesnosti protinani z úhlů vpřed
%Input:         pd  -   přesnost měření délek
%               B   -   souřadnice určovaných bodů
%               SS  -   souřadnice měřické sítě
%
%Output:        P1  -   matice všech kombinací bez vlivu podkladu   1.sl měření úhlu
%                                                                   2.sl měření úhlu
%               EX1_kov - kovariančí matice (složená matice 2x2, pořadí
%                                            řádek      -   měření úhlu z
%                                                           bodu 1.sl v matici P1
%                                            sloupeček  -   měření úhlu z
%                                                           bodu 2.sl v matici P1)
%
%% Výpočet
RAD=pi/200;
Ex1_kov=[];EX1_kov=[];
UP=[0,0;0,0];
D=[-1,1,0,0;0,0,-1,1];
for n=1:size(SS,1)
    for m=1:size(SS,1)
        if n~=m
        [sm, Spa] = cart2pol(B(1)-SS(n,2),B(2)-SS(n,3));
        [sm, Spb] = cart2pol(B(1)-SS(m,2),B(2)-SS(m,3));
        A1=[-((B(1)-SS(n,2))/(Spa^2)),((B(2)-SS(n,3))/(Spa^2))
            ((B(1)-SS(m,2))/(Spb^2)),-((B(2)-SS(m,3))/(Spb^2))];
        el=[pu,pu,pu,pu];
        EL=diag(el.^2);
        K=A1^(-1)*D;
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