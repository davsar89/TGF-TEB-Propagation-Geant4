clear all
close all
clc

%%

load('/Home/siv29/dsa030/Desktop/GEANT4/TGF-TEB-Propagation-Geant4/tests/data_t50_rd_15km.mat');
histogram('BinEdges',bins_rad_dist,'BinCounts',t50,'DisplayStyle','stairs','LineWidth',2);
hold on

load('/Home/siv29/dsa030/Desktop/GEANT4/TGF-TEB-Propagation-Geant4/tests/data_t50_rd_12km.mat');
histogram('BinEdges',bins_rad_dist,'BinCounts',t50,'DisplayStyle','stairs','LineWidth',2);

load('/Home/siv29/dsa030/Desktop/GEANT4/TGF-TEB-Propagation-Geant4/tests/data_t50_rd_10km.mat');
histogram('BinEdges',bins_rad_dist,'BinCounts',smooth(t50),'DisplayStyle','stairs','LineWidth',2);

load('/Home/siv29/dsa030/Desktop/GEANT4/TGF-TEB-Propagation-Geant4/tests/data_t50_t90_rd_15km_1MeV.mat');
histogram('BinEdges',bins_rad_dist,'BinCounts',smooth(t50),'DisplayStyle','stairs','LineWidth',2);

load('/Home/siv29/dsa030/Desktop/GEANT4/TGF-TEB-Propagation-Geant4/tests/data_t50_t90_rd_15km_1MeV.mat');
histogram('BinEdges',bins_rad_dist,'BinCounts',smooth(t90),'DisplayStyle','stairs','LineWidth',2);

legend('source at 15 km, t_{50}, E>300keV','source at 12 km t_{50}, E>300keV','source at 10 km t_{50}, E>300keV', ...
    'source at 15 km, t_{50}, E>1MeV','source at 15 km, t_{90}, E>1MeV')

xlabel('TGF ISS radial distance (km)')
ylabel('T_{50} or T_{90}  duration (micro-second)')
grid on