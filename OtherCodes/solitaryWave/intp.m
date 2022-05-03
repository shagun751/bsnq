clc
clear all
close all

load('solitaryWave.mat')

da = importdata('nodes.dat');

for i = 1:max(size(da))
    x = da(i,1);
    ind = find(abs(sx-x)<dx/10);   
    ssz = size(ind);
    if(ssz(2)==1)
        da(i,3) = se(ind);
        da(i,4) = sp(ind);
    end
end

plot(da(:,1),da(:,3))