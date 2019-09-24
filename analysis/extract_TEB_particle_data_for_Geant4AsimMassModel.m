clear all
close all
clc

%%

integration_rad_dist_interval = 50; % km
beaming_wanted = 30;
min_energy = 300; % keV

%%

simu_data=load('/Home/siv29/dsa030/Desktop/GEANT4/TGF-TEB-Propagation-Geant4/build/output_fram/simu_data_uniform_beaming_ratio_electron_photon.mat');
simu_data=simu_data.simu_data;

%%

beaming_tmp = simu_data(:,3);
ener_tmp = simu_data(:,9);

to_keep = beaming_tmp == beaming_wanted & ener_tmp > min_energy;
simu_data = simu_data(to_keep,:);

%%

type_tmp = simu_data(:,7);

simu.type = simu_data(:,7);

simu.energies = simu_data(:,9);

simu.momx = simu_data(:,17);
simu.momy = simu_data(:,18);
simu.momz = simu_data(:,19);

simu.lat = simu_data(:,11);
simu.lon = simu_data(:,12);

%%
leptons_only = simu_data(type_tmp==-11 | type_tmp==11,:);
lat_leptons = leptons_only(:,11);
lon_leptons = leptons_only(:,12);
mid_lon = find_histogram_max(lon_leptons).*ones(size(simu.lat));
mid_lat = find_histogram_max(lat_leptons).*ones(size(simu.lat));

%%

wgs84 = wgs84Ellipsoid('meters');
[dists, ~] = distance(mid_lat, mid_lon, simu.lat, simu.lon, wgs84); % m
dists = dists./1000; % m to km

to_keep2 = dists < integration_rad_dist_interval;

%%

simu.energies = simu.energies(to_keep2);
simu.momx = simu.momx(to_keep2);
simu.momy = simu.momy(to_keep2);
simu.momz = simu.momz(to_keep2);
simu.type = simu.type(to_keep2);

photons.energies = simu.energies(simu.type==22);
photons.momx = simu.momx(simu.type==22);
photons.momy = simu.momy(simu.type==22);
photons.momz = simu.momz(simu.type==22);
photons.type = simu.type(simu.type==22);

leptons.energies = simu.energies(simu.type~=22);
leptons.momx = simu.momx(simu.type~=22);
leptons.momy = simu.momy(simu.type~=22);
leptons.momz = simu.momz(simu.type~=22);
leptons.type = simu.type(simu.type~=22);

%%
M_out = [photons.energies photons.type photons.momx photons.momy photons.momz];
dlmwrite('photons.txt',M_out,'delimiter',' ','precision',6)

%%
% clear M_out;
% M_out = [leptons.energies leptons.type leptons.momx leptons.momy leptons.momz];
% dlmwrite('leptons.txt',M_out,'delimiter',' ','precision',6)













%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function max_of_distribution = find_histogram_max(values)

[N,EDGES] = histcounts(values,128);

[~,I] = max(N);

max_of_distribution = (EDGES(I)+EDGES(I+1))/2.0;

end



