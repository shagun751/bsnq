clc
clear all
close all

g=9.81;

h = 5;
cs = 0.85*sqrt(g*h);
p = g*h/cs/cs;
syms m n kh
m = tanh(kh)/kh;
n = 2*kh/sinh(2*kh);
soln = solve(2/p - m*(3-n));
kh = abs(double(soln));

L=2*pi*h/kh;

kh = 2*pi/L*h;
c = sqrt(g*L/2/pi*tanh(kh));
T = L/c;

hByL=h/L;
fprintf('f = %f\n',1/T);
fprintf('T = %f\n',T);
fprintf('L = %f\n',L);
fprintf('h = %f\n',h);
fprintf('h/L = %f\n',hByL);
fprintf('h/L0 = %f\n',h/(g/2/pi*T*T));
fprintf('kh = %f\n',2*pi*hByL);
fprintf('Fr = %f\n',c/sqrt(9.81*h));

