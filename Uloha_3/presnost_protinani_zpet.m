function [P1,EX1] = presnost_protinani_zpet(pu,B,SS)
%Protínání a jeho přesnosti protinani zpet
%Input:         pu  -   přesnost měření úhlů
%               B   -   souřadnice určovaných bodů
%               SS  -   souřadnice měřické sítě
%
%Output:        P1  -   matice všech kombinací bez vlivu podkladu   1.sl měření úhlu 
%                                                                   2.sl měření úhlu
%                                                                   3.sl měření úhlu
%               EX1_kov - kovariančí matice (složená matice 2x2, pořadí
%                                            řádek      -   z bodu 2.sl v matici P1, měření pouze úhlu
%                                            sloupeček  -   z bodu 3.sl v matici P1, měření pouze úhlu   )
%
%% Výpočet
    ID=[SS(:,1),SS(:,1),SS(:,1)];
    ind=size(SS,1);
    P1=[]; a=[]; b=[]; alf=[]; D=[-1,1,0;0,-1,1];EX1=[]; RAD=pi/200;
    for n=1:ind
        for m=1:ind
            if ID(n,1)~=ID(m,2)
                for o=1:ind
                    if ID(m,2)~=ID(o,3) && ID(n,1)~=ID(o,3)
                        i=find(ID(n,1)==SS(:,1));
                        j=find(ID(m,2)==SS(:,1));
                        k=find(ID(o,3)==SS(:,1));
    
                        [sm, Spa] = cart2pol(B(1)-SS(i,2),B(2)-SS(i,3));
                        [sm, Spb] = cart2pol(B(1)-SS(j,2),B(2)-SS(j,3));
                        [sm, Spc] = cart2pol(B(1)-SS(k,2),B(2)-SS(k,3));
                        A1=[((SS(j,2)-B(1))/(Spb^2))-((SS(i,2)-B(1))/(Spa^2)),-((SS(j,3)-B(2))/(Spb^2))+((SS(i,3)-B(2))/(Spa^2))
                            ((SS(k,2)-B(1))/(Spc^2))-((SS(j,2)-B(1))/(Spb^2)),-((SS(k,3)-B(2))/(Spc^2))+((SS(j,3)-B(2))/(Spb^2))];
                        el=[pu,pu,pu];
                        EL=diag(el.^2);
                        K=A1^(-1)*D;
                        Ex1=K*EL*K';                      
                        [a,b,alf] =par_el_chyb(Ex1);
                        sigx1=sqrt(Ex1(1)); sigy1=sqrt(Ex1(4));sigxy1=sqrt((Ex1(1)+Ex1(4))/2);
                        if alf<0
                            alf=alf+2*pi;
                        end
                        alf=alf-pi;
                        if alf<0
                            alf=alf+2*pi;
                        end
                        P1=[P1;[ID(n,1),ID(m,2),ID(o,3),sigx1,sigy1,sigxy1,alf/RAD,a,b]];
                        EX1=[EX1;Ex1];
                    end
                end
            end
        end
    end
end