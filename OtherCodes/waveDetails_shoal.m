
clc
clear all 
close all

g=9.81;

T=8;
h=8;%0.4572; %0.1524
H=1;
amp=H/2;
syms L
L=solve(L-(g/2/pi*T*T*tanh(2*pi/L*h)));
L=abs(double(L));
k = 2*pi/L;
n = 0.5*(1 + 2*k*h/sinh(2*k*h));

hByL=h/L;
fprintf('f = %f\n',1/T);
fprintf('T = %f\n',T);
fprintf('H = %f\n',H);
fprintf('Amp = %f\n',amp);
fprintf('L = %f\n',L);
fprintf('h = %f\n',h);
fprintf('h/L = %f\n',hByL);
fprintf('h/L0 = %f\n',h/(g/2/pi*T*T));
fprintf('kh = %f\n',2*pi*hByL);
fprintf('w = %f\n',2*pi/T);
fprintf('UMax = %f\n',amp*g*T/L);
fprintf('UMax x T = %f\n',amp*g*T/L*T);
fprintf('h/g/T/T = %f\n',h/g/T/T);


%% Shoal
h2 = 1;
syms L2
L2=solve(L2-(g/2/pi*T*T*tanh(2*pi/L2*h2)));
L2=abs(double(L2));
k2 = 2*pi/L2;
n2 = 0.5*(1 + 2*k2*h2/sinh(2*k2*h2));
H2 = H*sqrt(L/L2)*n/n2;
fprintf('\n\nh2 = %f\n',h2);
fprintf('H2 = %f\n',H2);
fprintf('L2 = %f\n',L2);
