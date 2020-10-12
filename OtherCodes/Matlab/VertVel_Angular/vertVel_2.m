clc
clear all

g=9.81;
d=0.7;
H=0.1;
T=2;
syms L;
L=solve(L-(g/2/pi*T*T*tanh(2*pi/L*d)));
L=abs(double(L));

k=2*pi/L;
w=2*pi/T;
fprintf('f = %f\n',1/T);
fprintf('T = %f\n',T);
fprintf('Height = %f\n',H);
fprintf('L = %f\n',L);
fprintf('d = %f\n',d);
fprintf('d/L = %f\n',d/L);
fprintf('d/L0 = %f\n',d/(g/2/pi*T*T));
fprintf('kd = %f\n',k*d);
fprintf('w = %f\n',w);
fprintf('C = %f\n\n',L/T);

t=1.25;
nx=1;
ny=sqrt(1 - nx*nx);
z0 = -0.35;

%%
% Wavemaker reference at 0, 0
% 'd' is still water depth
% 'z' is 0 at mean sea level

%[eta, p, q] = getDepInt(H, w, k, nx, ny, t, x0, y0);

dx=0.2;
[px, py]=meshgrid([-1:dx:1],[-1:dx:1]);
indX=7;
indY=6;
x0 = px(indY, indX);
y0 = py(indY, indX);
fprintf('( x0, y0 ) = ( %10.4f, %10.4f )\n', x0, y0);
fprintf('z = %10.4f \n', z0);

etM=zeros(size(px));
pM=zeros(size(px));
qM=zeros(size(px));

for i=1:size(px,1)
    for j = 1:size(px,2)
        x = px(i,j);
        y = py(i,j);
        [eta, p, q] = getDepInt(H, w, k, nx, ny, t, x, y);
        etM(i,j) = eta;
        pM(i,j) = p;
        qM(i,j) = q;
    end
end

umM = pM ./ (d);
vmM = qM ./ (d);
umdM = umM * d;
vmdM = vmM * d;

[~, p0, q0] = getDepInt(H, w, k, nx, ny, t, x0, y0);
[et0, u0, v0, w0, pr0] = getAiry(H, d, w, k, nx, ny, t, x0, y0, z0);

um0 = umM(indY, indX);
um0_x = ( umM(indY, indX+1) - umM(indY, indX-1) ) / 2/dx;
umd0_x = ( umdM(indY, indX+1) - umdM(indY, indX-1) ) / 2/dx;
um0_xx = ( umM(indY, indX+1) - 2*umM(indY, indX) + umM(indY, indX-1) ) / dx/dx;
umd0_xx = ( umdM(indY, indX+1) - 2*umdM(indY, indX) + umdM(indY, indX-1) ) / dx/dx;
um0_xxx = ( 0.5*umM(indY, indX+2) - umM(indY, indX+1) + umM(indY, indX-1) -0.5*umM(indY, indX-2) ) / dx^3;
umd0_xxx = ( 0.5*umdM(indY, indX+2) - umdM(indY, indX+1) + umdM(indY, indX-1) -0.5*umdM(indY, indX-2) ) / dx^3;

u0c = um0 - 0.5*d * ( umd0_xx ) ...
       + d*d/6 * ( um0_xx ) ...
       - z0 * ( umd0_xx ) ...
       - z0*z0/2 * ( um0_xx );

w0c = - ( umd0_x ) ...
      - z0 * ( um0_x ) ...
      + z0*d/2 * ( umd0_xxx ) ...
      - z0*d*d/6 * ( um0_xxx ) ...
      + z0*z0/2 * ( umd0_xxx ) ...
      + z0*z0*z0/6 * ( um0_xxx );

fprintf('( p, q ) = ( %10.4f, %10.4f ) \n', p0, q0 );
fprintf('u = %10.4f, %10.4f, %10.4f \n', u0, u0c, (u0c - u0)/u0*100 );
fprintf('w = %10.4f, %10.4f, %10.4f \n', w0, w0c, (w0c - w0)/w0*100 );



%%
function [eta, p, q] = getDepInt(H, w, k, nx, ny, t, x, y)

    % Wavemaker reference at 0, 0
    kx=k*nx;
    ky=k*ny;
    eta = H/2*sin(kx*x + ky*y - w*t);
    p = H/2 * w/k * nx * sin(kx*x + ky*y - w*t);
    q = H/2 * w/k * ny * sin(kx*x + ky*y - w*t);

end

function [eta, u, v, w, pr] = getAiry(H, d, w, k, nx, ny, t, x, y, z)

    % Wavemaker reference at 0, 0
    kx=k*nx;
    ky=k*ny;
    eta = H/2*sin(kx*x + ky*y - w*t);
    u = H*w/2 * cosh(k*(d+z)) / sinh(k*d) * nx * sin(kx*x + ky*y - w*t);
    v = H*w/2 * cosh(k*(d+z)) / sinh(k*d) * ny * sin(kx*x + ky*y - w*t);
    w = -H*w/2 * sinh(k*(d+z)) / sinh(k*d) * cos(kx*x + ky*y - w*t);
    pr = cosh(k*(d+z)) / cosh(k*d) * eta - z;

end


