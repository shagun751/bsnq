%% %%
% Reference
% Orszaghova, J., Borthwick, A.G.L., Taylor, P.H., 2012. 
% From the paddle to the beach - A Boussinesq shallow water 
% numerical wave tank based on Madsen and Sorensen equations. 
% Journal of Computational Physics 231, 328–344. 
% https://doi.org/10.1016/j.jcp.2011.08.028
%%%%%

%%
clc
clear all
close all

B = 1/15;
g=9.81;
h=1;
A=0.2;

C = getC(g,A,h);

qmax = C*A;

%odefnc(0,qmax)

opts = odeset('Reltol',1e-13,'AbsTol',1e-14,'Stats','on');
[x,q] = ode113(@(x,q) odefnc(x,q), [0 100], qmax);

plot(x,q)


%%
function C = getC(g,A,h)
    
    tmp1 = g*h*A*A*(A+3*h);
    tmp2 = 6*h*h*(A - h*log(1+A/h));
    
    C = sqrt(tmp1/tmp2);

end

%%
function dq = odefnc(x,q)

    B = 1/15;
    g = 9.81;
    h = 1;
    A = 0.2;
    C = 3.44066076038234;

    tmp1 = -q*C*C/h;
    tmp2 = g*q*q/(2*C*h);
    tmp3 = g*q*q*q/(6*C*C*h*h);
    tmp4 = C*C*C*log(1+q/C/h);
    
    tmp5 = ( B*g*h/C - C*(B+1/3) )*0.5;
    
    tmp6 = (tmp1 + tmp2 + tmp3 + tmp4) / tmp5;
   
    dq = sqrt(tmp6);

end