clc
clear all
close all

% dy/dt = f(t)
% y = sin(2*pi/a t)

a = 10;
dt = 0.1;
t0 = 6.0;
nt = 10;
dtf = dt/nt;
t = t0 + dtf*[0:nt];
ya = analytical_y(a,t);
fa = analytical_f(a,t);

yn=zeros(size(ya));
yn(1) = ya(1);

kC = getK(a,dt,t0);
rC = rk4_coeffs(dt);

for i = 1:nt
    rF = rk4_coeffs(t(i+1) - t(1));
    Tr = rF.C * rF.P * rC.Pinv * rC.Cinv
    ktil = Tr * kC;    
    yn(i+1) = yn(1) + getdy(t(i+1) - t(1), ktil);
end

yn2 = zeros(size(ya));
yn2(end) = yn(end);

kC = getK(a,dt,t(end));
rC = rk4_coeffs(dt);

for i = 0:nt-1    
    rF = rk4_coeffs(t(i+1) - t(end));
    ktil = rF.C * rF.P * rC.Pinv * rC.Cinv * kC;    
    yn2(i+1) = yn2(end) + getdy(t(i+1) - t(end), ktil);
end


%%
figure(1)
subplot(2,1,1)
hold on
plot(t, ya, 'k-', 'LineWidth', 3);
plot(t, yn, 'r--', 'LineWidth', 3);
%plot(t, yn2, 'b--', 'LineWidth', 3);
grid on
xticks(t)

subplot(2,1,2)
hold on
plot(t, (yn-ya)./ya, 'r-', 'LineWidth', 3);
%plot(t, (yn2-ya)./ya, 'b-', 'LineWidth', 3);
grid on
xticks(t)

%%
function [y] = analytical_y(a,t)

    y = sin(2*pi/a*t);

end

function [f] = analytical_f(a,t)

    f = 2*pi/a*cos(2*pi/a*t);

end

function [k] = getK(a,dt,t)
    k1 = analytical_f(a,t);
    k2 = analytical_f(a,t+dt/2);
    k3 = analytical_f(a,t+dt/2);
    k4 = analytical_f(a,t+dt);
    k=[k1;k2;k3;k4];
end

function [dy] = getdy(dt, k)
    dy = dt/6 * (k(1) + 2*k(2) + 2*k(3) + k(4));
end

function [rk4] = rk4_coeffs(dt)

    rk4.a21 = 0.5;
    rk4.a31 = 0;
    rk4.a32 = 0.5;
    rk4.a41 = 0;
    rk4.a42 = 0;
    rk4.a43 = 1;
    
    rk4.c22 = rk4.a21;
    rk4.c32 = rk4.a31 + rk4.a32;
    rk4.c33 = rk4.a32 * rk4.a21;
    rk4.c42 = rk4.a41 + rk4.a42 + rk4.a43;
    rk4.c43 = rk4.a43 * rk4.a32 + rk4.a43 * rk4.a31 + rk4.a42 * rk4.a21;
    rk4.c44 = rk4.a43 * rk4.a32 * rk4.a21;
    
    rk4.C = [   1  0   0   0;
                1   rk4.c22    0    0;
                1   rk4.c32     rk4.c33     0;
                1   rk4.c42     rk4.c43     rk4.c44];
                
    rk4.P = [   1  0   0   0;
                0   dt  0   0;
                0   0   dt*dt   0;
                0   0   0   dt^3];        
    
    rk4.Cinv = inv(rk4.C);
    rk4.Pinv = inv(rk4.P);
end

function [B] = getBCoef(dt)
    
    B = [   1   dt  dt*dt/2     dt*dt*dt/6;
            0   1   dt          dt*dt/2;
            0   0   1           dt;
            0   0   0           1];

end