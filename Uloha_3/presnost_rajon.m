function [P1,P2,EX1_kov,EX2_kov] = presnost_rajon(pu,pd,Sxy,B,SS)
%Protínání a jeho přesnosti rajón
%Input:         pu  -   přesnost měření úhlů
%               pd  -   přesnost měření délek
%               Sxy -   přesnost souřadnic podkladu
%               B   -   souřadnice určovaných bodů
%               SS  -   souřadnice měřické sítě
%
%Output:        P1  -   matice všech kombinací bez vlivu podkladu   1.sl stanovisko
%                                                                   2.sl orientace
%               P2  -   matice všech kombinací s vlivem podkladu    1.sl stanovisko
%                                                                   2.sl orientace
%               EX1_kov - kovariančí matice s měřením bez podkladu (složená matice 2x2, pořadí
%                                                                   řádek      -   stanovisko
%                                                                   sloupeček  -   orientace  )
%               EX2_kov - kovariančí matice s měřením s podkladem  (složená matice 2x2, pořadí
%                                                                   řádek      -   stanovisko
%                                                                   sloupeček  -   orientace  )
%%
RAD=pi/200;
D=[-1,1,0;0,0,1];
Ex1_kov=[];EX1_kov=[];
Ex2_kov=[];EX2_kov=[];
UP=[0,0;0,0];
for n=1:size(SS,1)
    for m=1:size(SS,1)
        if n~=m
        [sm, Sap] = cart2pol(B(1)-SS(n,2),B(2)-SS(n,3));
        [sm, Sab] = cart2pol(SS(n,2)-SS(m,2),SS(n,3)-SS(m,3));
        %Matice A1 měření
        A1=[-(B(1)-SS(n,2))/(Sap^2),(B(2)-SS(n,3))/(Sap^2);
                 ((B(2)-SS(n,3))/Sap),(B(1)-SS(n,2))/Sap];          
        %Matice A2 podkladu
        A2=[((B(1)-SS(n,2))/Sap^2)-((SS(m,2)-SS(n,2))/Sab^2),-((B(2)-SS(n,3))/Sap^2)+((SS(m,3)-SS(n,3))/Sab^2),((SS(m,2)-SS(n,2))/Sab^2),-((SS(m,3)-SS(n,3))/Sab^2);
            -((B(2)-SS(n,3))/Sap),-((B(1)-SS(n,2))/Sap),0,0];

        el=[pu;pu;pd];                                   %Přesnost měření
        ex=[Sxy;Sxy;Sxy;Sxy];
        K=A1^(-1)*D;
        L=A1^(-1)*A2;
        EL=diag(el.^2);
        EX=diag(ex.^2);
        Ex2=K*EL*K'+L*EX*L';                             %Kovarianční matice souřadnic s podkladem
        Ex1=K*EL*K';                                     %Kovarianční matice souřadnic bez podkladu
        [a11,b11,alf11] =par_el_chyb(Ex1);               %Parametry elipsy chyb
        a1(m,n)=a11; b1(m,n)=b11; alf1(m,n)=alf11;
        [a22,b22,alf22] =par_el_chyb(Ex2);               %Parametry elipsy chyb
        a2(m,n)=a22; b2(m,n)=b22; alf2(m,n)=alf22;
        sigx1(m,n)=sqrt(Ex1(1));sigy1(m,n)=sqrt(Ex1(4));sigxy1(m,n)=sqrt((Ex1(1)+Ex1(4))/2);
        sigx2(m,n)=sqrt(Ex2(1));sigy2(m,n)=sqrt(Ex2(4));sigxy2(m,n)=sqrt((Ex2(1)+Ex2(4))/2);
        KAR(m,n)=m;
        Ex1_kov=[Ex1_kov,Ex1];
        Ex2_kov=[Ex2_kov,Ex2];
        else
            Ex1_kov=[Ex1_kov,UP];
            Ex2_kov=[Ex2_kov,UP];
        end
        KAR1(m,n)=n;
    end
    EX1_kov=[EX1_kov;Ex1_kov];
    Ex1_kov=[];
    EX2_kov=[EX2_kov;Ex2_kov];
    Ex2_kov=[];
end

KAR1=reshape(KAR1,size(SS,1)*size(SS,1),1);
KAR=reshape(KAR,size(SS,1)*size(SS,1),1);
a1=reshape(a1,size(SS,1)*size(SS,1),1); b1=reshape(b1,size(SS,1)*size(SS,1),1); alf1=reshape(alf1,size(SS,1)*size(SS,1),1);
a2=reshape(a2,size(SS,1)*size(SS,1),1); b2=reshape(b2,size(SS,1)*size(SS,1),1); alf2=reshape(alf2,size(SS,1)*size(SS,1),1);
sigx1=reshape(sigx1,size(SS,1)*size(SS,1),1);sigy1=reshape(sigy1,size(SS,1)*size(SS,1),1);sigxy1=reshape(sigxy1,size(SS,1)*size(SS,1),1);
sigx2=reshape(sigx2,size(SS,1)*size(SS,1),1);sigy2=reshape(sigy2,size(SS,1)*size(SS,1),1);sigxy2=reshape(sigxy2,size(SS,1)*size(SS,1),1);

P1=[KAR1,KAR,sigx1,sigy1,sigxy1,alf1/RAD,a1,b1];
[r,s]=find(P1(:,2)==0);
P1(r,:)=[];
P2=[KAR1,KAR,sigx2,sigy2,sigxy2,alf2/RAD,a2,b2];
[r,s]=find(P2(:,2)==0);
P2(r,:)=[];
end