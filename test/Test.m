clc; clear; format long G
RAD=pi/200;

A=[1100,5000];
B=[1000,5000];

% plot([A(1),B(1)],[A(2),B(2)])

dx=B(2)-A(2);
dy=B(1)-A(1);

[smernik, delka] = cart2pol(dx, dy);

s=delka;
smer=300*RAD;

if smernik <0
    sm=smernik+2*pi;
end

sm=sm+smer;
if sm >2*pi
    sm=sm-2*pi;
end
P=[A(1)+sin(sm)*s,A(2)+cos(sm)*s];


sw=sqrt(1/(444444))*RAD;
ss=sqrt(1/(40000));

sw*1000/RAD
ss*1000

smxy1=sqrt((sw^2*s^2+ss^2)/2)*1000;

sw=(sw*1000)^2;
ss=(ss*1000)^2;

A1=[-((P(1)-A(1))/(s^2)),((P(2)-A(2))/(s^2))
    ((P(2)-A(2))/(s)),((P(1)-A(1))/(s))];

D=[-1,1,0;0,0,1];

A1=A1^(-1);
el=[sw,sw,ss];
EL=diag(el);
K=A1*D;

EX=K*EL*K';

c=sqrt((EX(1)-EX(4))^2);
a=sqrt((EX(1)+EX(4)+c)/2);
b=sqrt((EX(1)+EX(4)-c)/2);
alfa=atan2(0,EX(1)-EX(4))/RAD;

