clear all
close all
clc


%% RHESSI TGF distribution as function of radial distance (from Kjetil)

% histogram values for figure
% First row is the weak tgfs, second row is the search algorithm TGFs with WWLLN match
% 0-50km, 50-100, 100-150, etc
%
% 13 30 39 82 56 85 81 90 85 83 73 61 58 40 31 22 24 13 17 17
%
% 5 12 16 32 43 33 53 36 41 21 16 9 9 3 3 1 5 0 0 0


TGF_nb = [13 30 39 82 56 85 81 90 85 83 73 61 58 40 31 22 24 13 17 17 ];

rd_grid = 0:50:(length(TGF_nb))*50;

centers_rd_grid = (rd_grid(2:end)+rd_grid(1:end-1))/2.0;

histogram('BinEdges',rd_grid,'BinCounts',TGF_nb,'DisplayStyle','stairs','LineWidth',2);


%%

load('/Home/siv29/dsa030/Desktop/GEANT4/TGF-TEB-Propagation-Geant4/tests/t50_rad_dist_simu.mat');

t50_samples=[];

centers = (bins_rad_dist(2:end)+bins_rad_dist(1:end-1))/2.0;

funct_t50 = @(x) interp1(centers,t50,x);

funct_tgf_rad_dist_density = @(x) interp1(centers_rd_grid,TGF_nb,x);

sampled_rad_dist = sampleDist(funct_tgf_rad_dist_density,max(TGF_nb),100000,[0 1000],0);

t50_samples = funct_t50(sampled_rad_dist);


hist(t50_samples,0:20:500)