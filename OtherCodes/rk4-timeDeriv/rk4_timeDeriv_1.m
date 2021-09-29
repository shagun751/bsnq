clc
clear all
close all

% dy/dt = f(t)
% y = sin(2*pi/a t)

a = 10;
dt = 0.1;

tn = 6;
yn = analytical_y(a,tn);

tn1 = tn+dt;

k1 = analytical_f(a,tn);
k2 = analytical_f(a,tn+dt/2);
k3 = analytical_f(a,tn+dt/2);
k4 = analytical_f(a,tn+dt);
kMat = [k1;k2;k3;k4];

yn1 = yn + dt/6*(k1 + 2*k2 + 2*k3 + k4);
yn1a = analytical_y(a,tn1);

fprintf('t    %15s, %15s, %15s \n','y num', 'y ana', 'dydt ana');
fprintf('tn   %15.6f, %15.6f, %15.6f\n', yn, yn, analytical_f(a,tn));
fprintf('tn+1 %15.6f, %15.6f, %15.6f, %15.6f\n', yn1, yn1a, analytical_f(a,tn1), (yn1-yn1a)/yn1a*100);

rk4 = rk4_coeffs(dt);

figure(1)
t = [0:dt:10]';
y = analytical_y(a,t);
plot(t,y,'LineWidth',3)

rk4.Pinv * rk4.Cinv * kMat


function [y] = analytical_y(a,t)

    y = sin(2*pi/a*t);

end

function [f] = analytical_f(a,t)

    f = 2*pi/a*cos(2*pi/a*t);

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