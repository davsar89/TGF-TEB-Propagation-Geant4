clear all
close all
clc

yy=importdata('/Home/siv29/dsa030/Desktop/GEANT4/G4_TGF_TEB_PROPA_agile_test_sample/build/lala.txt');

[n,xout] = hist(yy,100) 
% rings = pi.* xout(2:end).^2 - pi.* xout(1:end-1).^2 ;
% xout(end)=[];
% n(end)=[];
% flux = n./rings;
plot(xout, n)

