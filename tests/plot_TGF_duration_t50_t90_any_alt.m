clear all
close all
clc

% simu_data=importdata('../build/fused_output.out');
simu_data = load('/Data/ift/ift_romfys1/dsarria/SIMULATION_DATA/TGF_t50_calc/data_for_t50_par.mat');
simu_data=simu_data.simu_data;

%%

energy = simu_data(:,9);
to_remove = energy<300;
simu_data(to_remove,:) = [];

tgf_alt = simu_data(:,2);
tgf_angle = simu_data(:,3);
rec_alt = simu_data(:,10);

grid.alt=unique(tgf_alt);
grid.angles=unique(tgf_angle);
grid.rec_alt = unique(rec_alt);

source_lat = unique(simu_data(:,21));
source_lon = unique(simu_data(:,22));

simu_data = simu_data(tgf_alt==15,:);

% part to remove becaus due to TEB : lat : -1.5 to 1  ; lon : 54.5 56
photons.lat = simu_data(:,11);
photons.lon = simu_data(:,12);
% to_remove = photons.lat > -2. & photons.lat < 1 & photons.lon>54.5 & photons.lon<56;
% simu_data(to_remove,:) = [];

to_remove = photons.lat > -6.80;
simu_data(to_remove,:) = [];

%%
figure
plot(photons.lon(1:10000),photons.lat(1:10000),'+');
% return;
%%

photons.rad_dist = simu_data(:,13);
photons.times = simu_data(:,8);

bin_size = 0.5; % km
bins_rad_dist = 0:bin_size:1200;

%%

length(photons.times)

for ii=1:length(bins_rad_dist)-1
    
    to_keep = photons.rad_dist>bins_rad_dist(ii) & photons.rad_dist<bins_rad_dist(ii+1);
    
    t50(ii) = get_t50(photons.times(to_keep));
    
    t90(ii) = get_t90(photons.times(to_keep));
    
end

t50 = fillmissing(t50,'linear');
t50(t50>3000)=0;

t90 = fillmissing(t90,'linear');
t90(t90>3000)=0;

%%
figure
histogram('BinEdges',bins_rad_dist,'BinCounts',t50,'DisplayStyle','stairs','LineWidth',2);

xlabel('TGF ISS radial distance (km)')
ylabel('T_{50} duration (micro-second)')

axis tight

%%
figure
histogram('BinEdges',bins_rad_dist,'BinCounts',t90,'DisplayStyle','stairs','LineWidth',2);

xlabel('TGF ISS radial distance (km)')
ylabel('T_{90} duration (micro-second)')

axis tight

%%

[N,EDGES] = histcounts(photons.rad_dist,bins_rad_dist);

centers = (EDGES(1:end-1)+EDGES(2:end))/2.0;
ring_areas = pi.*(EDGES(2:end).^2) - pi.*(EDGES(1:end-1).^2);

simulated_fluence = N ./ (ring_areas); % particles per km2
simulated_fluence = simulated_fluence ./ 1e10; % particles / cm2
nb_sam = 1;
source_fluences = simulated_fluence .* 1 ./ nb_sam;

%%
figure
histogram('BinEdges',EDGES,'BinCounts',source_fluences, ...
    'DisplayStyle','stairs','LineWidth',2);

%%
figure
histogram('BinEdges',EDGES,'BinCounts',source_fluences.*t50, ...
    'DisplayStyle','stairs','LineWidth',2);


%%

function t50 = get_t50(time_list)

% t01 = quantile(time_list,0.005); % to remove outliers
% time_list(time_list<t01) = [];
% t99 = quantile(time_list,0.995); % to remove outliers
% time_list(time_list>t99) = [];

t25 = quantile(time_list,0.25);
t75 = quantile(time_list,0.75);

t50 = t75-t25;

end

function t90 = get_t90(time_list)

% t01 = quantile(time_list,0.005); % to remove outliers
% time_list(time_list<t01) = [];
% t99 = quantile(time_list,0.995); % to remove outliers
% time_list(time_list>t99) = [];

t05 = quantile(time_list,0.05);
t95 = quantile(time_list,0.95);

t90 = t95-t05;

end


%%

% alt_wanted = 12;
% beam_wanted = 30;
% rad_dist_wanted = 330;
% date_wanted = '2018-Jun-21';
%
% fluence_out = get_fluence(rad_dist_wanted,date_wanted,alt_wanted,beam_wanted)

%%

% function rad_dist_out = get_fluence(rad_dist_wanted,date_wanted,alt_wanted,beam_wanted)
%
% if (strcmp(date_wanted,"2018-Jun-21"))
%     simu_data = load('simu_data_2018-Jun-21.mat');
% elseif(strcmp(date_wanted,"2018-Sep-05"))
%     simu_data = load('simu_data_2018-Sep-05.mat');
% elseif(strcmp(date_wanted,"2018-Oct-11"))
%     simu_data = load('simu_data_2018-Jun-21.mat');
% else
%     error('second input argument "date_wanted" is invalid. It should be "2018-Jun-21", "2018-Sep-05" or "2018-Oct-11"');
% end
%
% nVarargs = length(varargin);
%
% if nVarargs==2
%     i_a = 2;
%     i_ang = 3;
% elseif nVarargs==3
%     i_a = find(simu_data.gridd.alt==alt_wanted);
%     i_ang = 3;
% else
%     i_a = find(simu_data.gridd.alt==alt_wanted);
%     i_ang = find(simu_data.gridd.op_angle==beam_wanted);
% end
%
% simu_data = simu_data.simulation_data;
%
% fluences = squeeze(simu_data.source_fluences(i_a,i_ang,:));
%
% rad_dist_out = interp1(simu_data.centers,fluences,rad_dist_wanted);
%
% end

% save('simulation_data.mat','simulation_data','-v7.3');


%%
% exit

% wgs84 = wgs84Ellipsoid('meters');
% [dists, ~] = distance(tgf.lat*ones(size(photons.lat)), tgf.lon*ones(size(photons.lon)), photons.lat, photons.lon, wgs84); % m
% dists = dists./1000; % m to km

function FolderInfo = dir2(varargin)
ignored_folder_names ={'.','..'};
% A custom dir function that does not list . and ..

if nargin == 0
    name = '.';
elseif nargin == 1
    name = varargin{1};
else
    error('Too many input arguments.')
end

FolderInfo = dir(name);

inds = [];
n    = 0;
k    = 1;

while n < length(ignored_folder_names) && k <= length(FolderInfo)
    if any(strcmp(FolderInfo(k).name, ignored_folder_names))
        inds(end + 1) = k;
        n = n + 1;
    end
    k = k + 1;
end

FolderInfo(inds) = [];
end


%%

function data_temp = better_importdata(full_file_path,number_of_columns)
fid = fopen(full_file_path);

text_pattern_to_scan = repmat('%f ',1,number_of_columns);

C = textscan(fid,text_pattern_to_scan,'CommentStyle','#');
% data_temp = single([C{:}]); % single precision to use les memory
data_temp = [C{:}]; % single precision to use les memory

fclose(fid);
end
