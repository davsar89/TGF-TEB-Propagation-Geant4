clear all
close all
clc

addpath('./export_fig/');

beaming_type = "gaussian";
% beaming_type = "uniform";

%%

real_measurement_fluence = 170/150; % particles per cm2; 150 cm2 effective are, 170 counts repported for the event

%%

% yy = importdata('../build/fused.out');
% yy=load('/Home/siv29/dsa030/Desktop/GEANT4/TGF-TEB-Propagation-Geant4/build/output_fram/simu_data_filtered.mat');
yy=load('/Data/ift/ift_romfys1/dsarria/SIMULATION_DATA/TEB/simu_data_180602/simu_data_contour_fluence_alt_ang_gaussian.mat');
% yy=load('/Data/ift/ift_romfys1/dsarria/TEB/simu_data_180602/simu_data_contour_fluence_alt_ang.mat');
yy=yy.simu_data;

lat_temp = yy(:,11);
yy(lat_temp<2,:)=[];

simu.seed = yy(:,1);
simu.type = yy(:,7);
simu.nb_posi = sum(simu.type==-11);
simu.nb_elec = sum(simu.type==11);
simu.ratio_ee = simu.nb_posi/simu.nb_elec*100;
simu.nb_event = yy(:,5);
simu.lat = yy(:,11);
simu.lon = yy(:,12);

% sanity check
if (sum(simu.nb_event==0)>0)
    error('nb_event is 0 or at least one record. It should be >=1.')
end

source.alt = yy(:,2);
source.opening_angle = yy(:,3);

%%

griddd.alt = unique(source.alt)';
griddd.opening_angle = unique(source.opening_angle)';
griddd.opening_angle(griddd.opening_angle==15)=[];
griddd.opening_angle(griddd.opening_angle==5)=[];
griddd.opening_angle(griddd.opening_angle==60)=[];
griddd

%%
nb_sampled = get_ini_sampled(griddd, yy);

%%
source_fluences = get_required_source_fluence(griddd, source, yy, nb_sampled, real_measurement_fluence);

e_e_ratio = get_positron_electron_ratio(griddd, source, yy);

source_fluences_log10 = log10(source_fluences);
% source_fluences_log10=fillmissing(source_fluences_log10,'linear');

source_fluences_log10_smoothed=source_fluences_log10;

% for ii=1:size(source_fluences_log10,1)
%     source_fluences_log10_smoothed(ii,:) = smoothdata(source_fluences_log10_smoothed(ii,:));
% end
% %
% for ii=1:size(source_fluences_log10,2)
%     source_fluences_log10_smoothed(:,ii) = smoothdata(source_fluences_log10_smoothed(:,ii));
% end

[X,Y] = meshgrid(griddd.opening_angle,griddd.alt);

minnn = round(min(source_fluences_log10_smoothed(:)),1);
maxxx = round(max(source_fluences_log10_smoothed(:)),1);
levels = minnn:0.1:maxxx ;

[c,h] = contour(X,Y,source_fluences_log10_smoothed, levels, 'ShowText', 'on');


ylabel('altitude (km)','interpreter','latex','fontsize',17)

set(gcf,'renderer','Painters')
% print -dpdf -painters -bestfit -r0 fluence_vs_alt_angle.pdf

if strcmp(beaming_type, "gaussian")
    xlabel('standard deviation $\sigma_\theta$ (degrees)','interpreter','latex','fontsize',17)
    title('Gaussian beaming')
    axis([10 40 10 16])
else
    axis([10 50 10 16])
    title('Uniform beaming')
    xlabel('$\theta$ (degrees)','interpreter','latex','fontsize',17)
end
grid on

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function e_e_ratio = get_positron_electron_ratio(gridd, source, yy)

wgs84 = wgs84Ellipsoid('meters');

e_e_ratio = zeros(length(gridd.alt),length(gridd.opening_angle));

for i_a = 1:length(gridd.alt)
    for i_ang = 1:length(gridd.opening_angle)
        
        idx_slc = source.alt == gridd.alt(i_a) & source.opening_angle == gridd.opening_angle(i_ang);
        yy_kept = yy(idx_slc,:);
        
        simu.lat = yy_kept(:,11);
        simu.lon = yy_kept(:,12);
        
        simu.type = yy_kept(:,7);
        
        %         teb.center.lat = median(simu.lat);
        %         teb.center.lon = median(simu.lon);
        q1_lat = quantile(simu.lat,0.25);
        q2_lat = quantile(simu.lat,0.75);
        q1_lon = quantile(simu.lon,0.25);
        q2_lon = quantile(simu.lon,0.75);
        
        qmid_lon = (q2_lon+q1_lon)/2.0;
        qmid_lat = (q2_lat+q1_lat)/2.0;
        
        nb_inside = sum(simu.lat>q1_lat & simu.lat<q2_lat & simu.lon>q1_lon & simu.lon<q2_lon);
        
        [dist_lat, ~] = distance(q1_lat, qmid_lon, q2_lat, qmid_lon,wgs84); % m
        dist_lat = dist_lat/1000; % m to km
        
        [dist_lon, ~] = distance(qmid_lat, q1_lon, qmid_lat, q2_lon,wgs84); % m
        dist_lon = dist_lon/1000; % m to km
        
        disp([num2str(gridd.alt(i_a)) ' ' num2str(gridd.opening_angle(i_ang)) ' ' num2str(nb_inside)]);
        
        nb_elec = sum(simu.type==11);
        nb_posi = sum(simu.type==-11);
        
        e_e_ratio(i_a,i_ang) = nb_posi/nb_elec*100;
        
    end
end

end

%%
function source_fluences = get_required_source_fluence(grid, source, yy, nb_sampled, real_measurement_fluence)

wgs84 = wgs84Ellipsoid('meters');

source_fluences = zeros(length(grid.alt),length(grid.opening_angle));

for i_a = 1:length(grid.alt)
    for i_ang = 1:length(grid.opening_angle)
        
        idx_slc = source.alt == grid.alt(i_a) & source.opening_angle == grid.opening_angle(i_ang);
        yy_kept = yy(idx_slc,:);
        
        simu.lat = yy_kept(:,11);
        simu.lon = yy_kept(:,12);
        
        %         teb.center.lat = median(simu.lat);
        %         teb.center.lon = median(simu.lon);
        q1_lat = quantile(simu.lat,0.25);
        q2_lat = quantile(simu.lat,0.75);
        q1_lon = quantile(simu.lon,0.25);
        q2_lon = quantile(simu.lon,0.75);
        
        qmid_lon = (q2_lon+q1_lon)/2.0;
        qmid_lat = (q2_lat+q1_lat)/2.0;
        
        nb_inside = sum(simu.lat>q1_lat & simu.lat<q2_lat & simu.lon>q1_lon & simu.lon<q2_lon);
        
        [dist_lat, ~] = distance(q1_lat, qmid_lon, q2_lat, qmid_lon,wgs84); % m
        dist_lat = dist_lat/1000; % m to km
        
        [dist_lon, ~] = distance(qmid_lat, q1_lon, qmid_lat, q2_lon,wgs84); % m
        dist_lon = dist_lon/1000; % m to km
        
        disp([num2str(grid.alt(i_a)) ' ' num2str(grid.opening_angle(i_ang)) ' ' num2str(nb_inside)]);
        
        simulated_fluence = nb_inside / (dist_lat*dist_lon); % particles per km2
        simulated_fluence = simulated_fluence ./ 1e10; % particles / cm2
        source_fluences(i_a,i_ang) = real_measurement_fluence * nb_sampled(i_a,i_ang) / simulated_fluence;
        
    end
end

end

%%
function nb_sampled = get_ini_sampled(grid, yy)

nb_sampled=zeros(length(grid.alt),length(grid.opening_angle));

alts = yy(:,2);
opening_angles = yy(:,3);

% getting initial number of photons sampled per parameter set
for i_a = 1:length(grid.alt)
    for i_ang = 1:length(grid.opening_angle)
        
        idx_slc = alts == grid.alt(i_a) & opening_angles == grid.opening_angle(i_ang);
        
        yy_kept = yy(idx_slc,:);
        seed = yy_kept(:,1);
        nb_event = yy(:,5);
        
        list_seeds = unique(seed);
        
        nb_sampled(i_a,i_ang) = 0;
        
        for i_seed = 1:length(list_seeds)
            idx_seed = list_seeds(i_seed)==seed;
            nb_sampled(i_a,i_ang) = nb_sampled(i_a,i_ang) + max(nb_event(idx_seed));
        end
        
    end
end
end

