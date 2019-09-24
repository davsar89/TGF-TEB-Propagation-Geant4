clear all
close all
clc


load('/Home/siv29/dsa030/Desktop/GEANT4/TGF-TEB-Propagation-Geant4/analysis/inside_50km_square.mat');

histogram('BinEdges',EDGES,'BinCounts',spec,'DisplayStyle','stairs','LineWidth',2);
hold on

load('/Home/siv29/dsa030/Desktop/GEANT4/TGF-TEB-Propagation-Geant4/analysis/inside_10km_square.mat');

histogram('BinEdges',EDGES,'BinCounts',spec*18,'DisplayStyle','stairs','LineWidth',2);
hold on

set(gca,'xscale','log')
set(gca,'yscale','log')