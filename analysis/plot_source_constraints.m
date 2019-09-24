clear all
close all
clc

addpath('./export_fig/');

%%

real_measurement_fluence = 170/150; % particles per cm2; 150 cm2 effective are, 170 counts repported for the event

%%

% yy = importdata('../build/fused.out');
yy=load('../build/simu_data_contour_fluence_alt_ang.mat');
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

gridd.alt = unique(source.alt)';
gridd.opening_angle = unique(source.opening_angle)';
gridd.opening_angle(gridd.opening_angle==15)=[];
gridd

%%
nb_sampled = get_ini_sampled(gridd,source, yy);

%%
source_fluences = get_required_source_fluence(gridd, source, yy, nb_sampled, real_measurement_fluence);

source_fluences_log10 = log10(source_fluences);
% source_fluences_log10=fillmissing(source_fluences_log10,'linear');

source_fluences_log10_smoothed=source_fluences_log10;

for ii=1:size(source_fluences_log10,1)
    source_fluences_log10_smoothed(ii,:) = smoothdata(source_fluences_log10_smoothed(ii,:));
end
% 
for ii=1:size(source_fluences_log10,2)
    source_fluences_log10_smoothed(:,ii) = smoothdata(source_fluences_log10_smoothed(:,ii));
end

[X,Y] = meshgrid(gridd.opening_angle,gridd.alt);

levels =  17:0.1:18.6;

[c,h] = contour(X,Y,source_fluences_log10_smoothed, levels, 'ShowText', 'on');

xlabel('Opening angle (degrees)','interpreter','latex','fontsize',17)
ylabel('altitude (km)','interpreter','latex','fontsize',17)

set(gcf,'renderer','Painters')
% print -dpdf -painters -bestfit -r0 fluence_vs_alt_angle.pdf


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function source_fluences = get_required_source_fluence(grid, source, yy, nb_sampled, real_measurement_fluence)
tgf.center.lat = 11.01;
tgf.center.lon = -95.40;

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
        
        simulated_fluence = nb_inside / (dist_lat*dist_lon); % particles per km2
        simulated_fluence = simulated_fluence ./ 1e10; % particles / cm2
        source_fluences(i_a,i_ang) = real_measurement_fluence * nb_sampled(i_a,i_ang) / simulated_fluence;
        
    end
end

end
%%
function nb_sampled = get_ini_sampled(grid,source, yy)

nb_sampled=zeros(length(grid.alt),length(grid.opening_angle));

% getting initial number of photons sampled per parameter set
for i_a = 1:length(grid.alt)
    for i_ang = 1:length(grid.opening_angle)
        
        idx_slc = source.alt == grid.alt(i_a) & source.opening_angle == grid.opening_angle(i_ang);
        yy_kept = yy(idx_slc,:);
        simu2.seed = yy_kept(:,1);
        simu2.nb_event = yy(:,5);
        
        list_seeds = unique(simu2.seed);
        
        nb_sampled(i_a,i_ang) = 0;
        
        for i_seed = 1:length(list_seeds)
            idx_seed = list_seeds(i_seed)==simu2.seed;
            nb_sampled(i_a,i_ang) = nb_sampled(i_a,i_ang) + max(simu2.nb_event(idx_seed));
        end
        
    end
end
end

