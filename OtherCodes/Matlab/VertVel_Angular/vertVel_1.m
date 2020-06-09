clc
clear all

g=9.81;
d=0.7;
H=0.1;
T=2;
syms L;
L=solve(L-(g/2/pi*T*T*tanh(2*pi/L*d)));
L=abs(double(L));

dByL=d/L;
fprintf('f = %f\n',1/T);
fprintf('T = %f\n',T);
fprintf('Height = %f\n',H);
fprintf('L = %f\n',L);
fprintf('d = %f\n',d);
fprintf('d/L = %f\n',dByL);
fprintf('d/L0 = %f\n',d/(g/2/pi*T*T));
fprintf('kd = %f\n',2*pi*dByL);
fprintf('w = %f\n',2*pi/T);
fprintf('UMax = %f\n',H/2*g*T/L);
fprintf('UMax x T = %f\n',H*g*T/L*T);

%%
%%Wavemaker reference at 0, 0
x0=1;
y0=1;
t=0;
nx=1/sqrt(2);
ny=sqrt(1-nx*nx);
k=2*pi/L;
w=2*pi/T;

[eta, p, q] = getVars(H, w, k, nx, ny, t, x0, y0);


%%
function [eta, p, q] = getVars(H, w, k, nx, ny, t, x, y)

    % Wavemaker reference at 0, 0
    kx=k*nx;
    ky=k*ny;
    eta = H/2*sin(kx*x + ky*y - w*t);
    p = H/2 * w/k * nx * sin(kx*x + ky*y - w*t);
    q = H/2 * w/k * ny * sin(kx*x + ky*y - w*t);

end


