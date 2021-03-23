clc
clear all
%close all

g=9.81;

h = 5;
c = 0.85*sqrt(g*h);
amp = 1;
syms kh

kh=solve( tanh(kh)/kh - 2*c*c/g/h/(3 - 2*kh/sinh(2*kh)) );
kh=abs(double(kh));

n = 2*kh/sinh(2*kh);

th = rad2deg(acos(sqrt(8*(1-n))/(3-n)));

fprintf('kh = %f\n',kh);
fprintf('theta = %f\n',th);

