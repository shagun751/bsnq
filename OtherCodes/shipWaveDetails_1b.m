clc
clear all
close all

g=9.81;

h = 5;
c_g = 0.85*sqrt(g*h);
amp = 1;
syms L c kh
kh = 2*pi/L*h;
c = sqrt(g*L/2/pi*tanh(kh));
L=solve(c_g - c/2*(1 + 2*kh/sinh(kh)));
L=abs(double(L));

kh = 2*pi/L*h;
c = sqrt(g*L/2/pi*tanh(kh));
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

