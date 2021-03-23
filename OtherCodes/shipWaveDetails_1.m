clc
clear all
close all

g=9.81;

h = 5;
c = 0.85*sqrt(g*h);
amp = 1;
syms L
L=solve(c*c-(g/2/pi*L*tanh(2*pi/L*h)));
L=abs(double(L));

T = L/c;

hByL=h/L;
fprintf('f = %f\n',1/T);
fprintf('T = %f\n',T);
fprintf('Amp = %f\n',amp);
fprintf('L = %f\n',L);
fprintf('h = %f\n',h);
fprintf('h/L = %f\n',hByL);
fprintf('h/L0 = %f\n',h/(g/2/pi*T*T));
fprintf('kh = %f\n',2*pi*hByL);
fprintf('Fr = %f\n',c/sqrt(9.81*h));
fprintf('UMax x T = %f\n',amp*g*T/L*T);

