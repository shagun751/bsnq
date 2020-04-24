clc
clear all
close all

% x=[1,3,2,2,2.5,1.5];
% y=[1,1,3,1,2,2];

x=[1,3,-5,2,-1,-2];
y=[1,3,7,2,5,4];

ea=[];
na=[];
pha=[];
xm=[];
ym=[];

ne=11;
nn=11;
de=1/(ne-1);

k=0;
for i=1:ne
    e=(i-1)*de;        
    dn=(1-e)/(nn-1);
    for j=1:nn        
        n=(j-1)*dn;
        k=k+1;
        ea(k)=e;
        na(k)=n;
        phi=[1 - 3*e - 3*n + 2*e^2 + 2*n^2 + 4*e*n,...
            2*e^2 - e,...
            2*n^2 - n,...
            4*e - 4*e^2 - 4*e*n,...
            4*e*n,...
            4*n - 4*n^2 - 4*e*n];
        pha(k,:)=phi;
        xm(k)=phi*x';
        ym(k)=phi*y';
    end
end

figure(1)
hold on
plot(x,y,'+')
plot([x(1),x(2)],[y(1),y(2)])
plot([x(2),x(3)],[y(2),y(3)])
plot([x(3),x(1)],[y(3),y(1)])
plot(xm,ym,'.','MarkerSize',20)

figure(2)
plot(ea,na,'.')

for i=1:6
    figure(10+i)
    plot3(ea,na,pha(:,i),'.')
end