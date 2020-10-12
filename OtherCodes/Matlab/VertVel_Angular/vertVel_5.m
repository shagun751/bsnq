clc
clear all
close all

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
wvAngDeg = 30;
if(wvAngDeg==90)
    nx = 0;
else    
    nx=cos(deg2rad(wvAngDeg));
end
ny=sqrt(1 - nx*nx);

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
fprintf('( x0, y0 ) = ( %10.4f, %10.4f )\n\n', x0, y0);

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

data = zeros(15,10);
data(:,1) = [-0.70:0.05:0]';

for i = 1:size(data,1)
    z0 = data(i,1);
    [et0, u0, v0, w0, pr0] = getAiry(H, d, w, k, nx, ny, t, x0, y0, z0);
    [u0c, v0c, w0c, uErr, vErr, wErr] = calcAll(d, indX, indY, et0, u0, v0, w0, pr0, p0, q0, umM, vmM, umdM, vmdM, dx, z0);
    data(i,2:10)=[u0, u0c, uErr, v0, v0c, vErr, w0, w0c, wErr];
end

figure(1)
subplot(1,4,1)
hold on
plot(data(:,2), data(:,1), 'k', 'LineWidth', 2)
plot(data(:,3), data(:,1), 'r', 'LineWidth', 2)
legend('Airy','Calc', 'Location', 'southeast')
grid on
set(gca, 'GridAlpha', 1, 'GridLineStyle','--')
title('u Vel (m/s)')
xlabel('u (m/s)')
ylabel('z (m)')


figure(1)
subplot(1,4,2)
hold on
plot(data(:,5), data(:,1), 'k', 'LineWidth', 2)
plot(data(:,6), data(:,1), 'r', 'LineWidth', 2)
legend('Airy','Calc', 'Location', 'southeast')
grid on
set(gca, 'GridAlpha', 1, 'GridLineStyle','--')
title('v Vel (m/s)')
xlabel('v (m/s)')
ylabel('z (m)')

figure(1)
subplot(1,4,3)
hold on
plot(data(:,8), data(:,1), 'k', 'LineWidth', 2)
plot(data(:,9), data(:,1), 'r', 'LineWidth', 2)
legend('Airy','Calc', 'Location', 'southeast')
grid on
set(gca, 'GridAlpha', 1, 'GridLineStyle','--')
title('w Vel (m/s)')
xlabel('w (m/s)')
ylabel('z (m)')


figure(1)
subplot(1,4,4)
hold on
plot(data(:,4), data(:,1), 'k', 'LineWidth', 2)
plot(data(:,7), data(:,1), 'r', 'LineWidth', 2)
plot(data(:,10), data(:,1), 'b', 'LineWidth', 2)
legend('u', 'v', 'w', 'Location', 'southeast')
grid on
set(gca, 'GridAlpha', 1, 'GridLineStyle','--')
title('% Errors')
xlabel('%')
ylabel('z (m)')


sgtitle(strcat('T = 2s, d = 0.7m, H = 0.1m, t=1.25s, wvAngDeg = ', num2str(wvAngDeg) ))

%%

function [u0c, v0c, w0c, uErr, vErr, wErr] = calcAll(d, indX, indY, et0, u0, v0, w0, pr0, p0, q0, umM, vmM, umdM, vmdM, dx, z0)
    fprintf('z = %10.4f \n', z0);

    um0 = umM(indY, indX);
    um0_x = ( umM(indY, indX+1) - umM(indY, indX-1) ) / 2/dx;
    umd0_x = ( umdM(indY, indX+1) - umdM(indY, indX-1) ) / 2/dx;
    um0_xx = ( umM(indY, indX+1) - 2*umM(indY, indX) + umM(indY, indX-1) ) / dx/dx;
    umd0_xx = ( umdM(indY, indX+1) - 2*umdM(indY, indX) + umdM(indY, indX-1) ) / dx/dx;
    um0_xxx = ( 0.5*umM(indY, indX+2) - umM(indY, indX+1) + umM(indY, indX-1) -0.5*umM(indY, indX-2) ) / dx^3;
    umd0_xxx = ( 0.5*umdM(indY, indX+2) - umdM(indY, indX+1) + umdM(indY, indX-1) -0.5*umdM(indY, indX-2) ) / dx^3;
       
    um0_y = ( umM(indY+1, indX) - umM(indY-1, indX) ) / 2/dx;
    umd0_y = ( umdM(indY+1, indX) - umdM(indY-1, indX) ) / 2/dx;
    um0_yy = ( umM(indY+1, indX) - 2*umM(indY, indX) + umM(indY-1, indX) ) / dx/dx;
    umd0_yy = ( umdM(indY+1, indX) - 2*umdM(indY, indX) + umdM(indY-1, indX) ) / dx/dx;
    um0_yyy = ( 0.5*umM(indY+2, indX) - umM(indY+1, indX) + umM(indY-1, indX) -0.5*umM(indY-2, indX) ) / dx^3;
    umd0_yyy = ( 0.5*umdM(indY+2, indX) - umdM(indY+1, indX) + umdM(indY-1, indX) -0.5*umdM(indY-2, indX) ) / dx^3;

    vm0 = vmM(indY, indX);
    vm0_y = ( vmM(indY+1, indX) - vmM(indY-1, indX) ) / 2/dx;
    vmd0_y = ( vmdM(indY+1, indX) - vmdM(indY-1, indX) ) / 2/dx;
    vm0_yy = ( vmM(indY+1, indX) - 2*vmM(indY, indX) + vmM(indY-1, indX) ) / dx/dx;
    vmd0_yy = ( vmdM(indY+1, indX) - 2*vmdM(indY, indX) + vmdM(indY-1, indX) ) / dx/dx;
    vm0_yyy = ( 0.5*vmM(indY+2, indX) - vmM(indY+1, indX) + vmM(indY-1, indX) -0.5*vmM(indY-2, indX) ) / dx^3;
    vmd0_yyy = ( 0.5*vmdM(indY+2, indX) - vmdM(indY+1, indX) + vmdM(indY-1, indX) -0.5*vmdM(indY-2, indX) ) / dx^3;
    
    vm0_x = ( vmM(indY, indX+1) - vmM(indY, indX-1) ) / 2/dx;
    vmd0_x = ( vmdM(indY, indX+1) - vmdM(indY, indX-1) ) / 2/dx;
    vm0_xx = ( vmM(indY, indX+1) - 2*vmM(indY, indX) + vmM(indY, indX-1) ) / dx/dx;
    vmd0_xx = ( vmdM(indY, indX+1) - 2*vmdM(indY, indX) + vmdM(indY, indX-1) ) / dx/dx;
    vm0_xxx = ( 0.5*vmM(indY, indX+2) - vmM(indY, indX+1) + vmM(indY, indX-1) -0.5*vmM(indY, indX-2) ) / dx^3;
    vmd0_xxx = ( 0.5*vmdM(indY, indX+2) - vmdM(indY, indX+1) + vmdM(indY, indX-1) -0.5*vmdM(indY, indX-2) ) / dx^3;
    
    um0_xy = ( umM(indY, indX+1) - umM(indY, indX-1) ) * ( umM(indY+1, indX) - umM(indY-1, indX) ) / 4/dx/dx;
    umd0_xy = ( umdM(indY, indX+1) - umdM(indY, indX-1) ) * ( umdM(indY+1, indX) - umdM(indY-1, indX) ) / 4/dx/dx;
    vm0_xy = ( vmM(indY, indX+1) - vmM(indY, indX-1) ) * ( vmM(indY+1, indX) - vmM(indY-1, indX) ) / 4/dx/dx;
    vmd0_xy = ( vmdM(indY, indX+1) - vmdM(indY, indX-1) ) * ( vmdM(indY+1, indX) - vmdM(indY-1, indX) ) / 4/dx/dx;

    k = 1;
    u0c = um0 - 0.5*d * ( umd0_xx + vmd0_xx + k*umd0_xy - k*vmd0_xy ) ...
           + d*d/6 * ( um0_xx + vm0_xx + k*um0_xy - k*vm0_xy ) ...
           - z0 * ( umd0_xx + vmd0_xx + k*umd0_xy - k*vmd0_xy ) ...
           - z0*z0/2 * ( um0_xx + vm0_xx + k*um0_xy - k*vm0_xy );

    v0c = vm0 - 0.5*d * ( vmd0_yy + umd0_yy - k*umd0_xy + k*vmd0_xy ) ...
           + d*d/6 * ( vm0_yy + um0_yy - k*um0_xy + k*vm0_xy ) ...
           - z0 * ( vmd0_yy + umd0_yy - k*umd0_xy + k*vmd0_xy ) ...
           - z0*z0/2 * ( vm0_yy + um0_yy - k*um0_xy + k*vm0_xy );

    w0c = - ( umd0_x ) ...
          - z0 * ( um0_x ) ...
          + z0*d/2 * ( umd0_xxx ) ...
          - z0*d*d/6 * ( um0_xxx ) ...
          + z0*z0/2 * ( umd0_xxx ) ...
          + z0*z0*z0/6 * ( um0_xxx );
    w0c = w0c - ( vmd0_y ) ...
          - z0 * ( vm0_y ) ...
          + z0*d/2 * ( vmd0_yyy ) ...
          - z0*d*d/6 * ( vm0_yyy ) ...
          + z0*z0/2 * ( vmd0_yyy ) ...
          + z0*z0*z0/6 * ( vm0_yyy );
     
    uErr = (u0c - u0)/u0*100;
    vErr = (v0c - v0)/v0*100;
    wErr = (w0c - w0)/w0*100;
    
    fprintf('( p, q ) = ( %10.4f, %10.4f ) \n', p0, q0 );
    fprintf('u = %10.4f, %10.4f, %10.4f \n', u0, u0c, uErr );
    fprintf('v = %10.4f, %10.4f, %10.4f \n', v0, v0c, vErr );
    fprintf('w = %10.4f, %10.4f, %10.4f \n\n', w0, w0c, wErr );

end 

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



