function [a1,b1,alf1] = par_el_chyb(EX1)
%Funkce na vytvoření elipsy chyb
%Input:        EX1      -   kovarianční matice vypočteného bodu
%
%Output:        a1      -   poloosa a elipsy chyb
%               b1      -   poloosa b elipsy chyb
%               alf1    -   úhel stočení alfa
%
%%
c1=sqrt((EX1(1)-EX1(4))^2+4*EX1(2)^2);
a1=sqrt((EX1(1)+EX1(4)+c1)/2);
b1=sqrt((EX1(1)+EX1(4)-c1)/2);
alf1=atan2(2*EX1(2),EX1(1)-EX1(4))/2;
if alf1<0
    alf1=alf1+2*pi;
end
end