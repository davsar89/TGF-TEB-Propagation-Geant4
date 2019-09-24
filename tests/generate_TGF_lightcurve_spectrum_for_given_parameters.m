clear all
close all
clc

E = wgs84Ellipsoid('meter');

TGFsource.lat = -13;
TGFsource.lon = 32;

%%

yy = importdata('/Home/siv29/dsa030/Desktop/GEANT4/GEANT4-TGF-Propagation/build/output/detParticles_478_400_15_40_Uniform_0.out');

ecef(:,1) = yy(:,11)*1000; % km to m
ecef(:,2) = yy(:,12)*1000;
ecef(:,3) = yy(:,13)*1000;

lla = ecef2lla(ecef);

[RD,az] = distance(lla(:,1),lla(:,2),TGFsource.lat,TGFsource.lon,E); % arc-length output is radial distance

RD = RD/1000; % m to km