clc
close all

for i=1:3
    figure(i)
    hold on
    plot(x,uderv(:,i),'LineWidth',2)
    plot(x,udervA(:,i),'LineWidth',2)
    legend('MLS','Analytical')
    title(strcat(num2str(i),'rd Derivative of sin(x). Mesh 0.2032'))
    grid on
    set(gca,'GridAlpha',1,'GridLineStyle','--')
end