clear all
close all
clc

%%

FolderInfo = dir2('./output_ascii/');

nb_files=length(FolderInfo);

simu_data0 = cell(1,nb_files);

seeds=1:length(FolderInfo);

%%
parfor ii=1:length(FolderInfo)
    yy_tmp = better_importdata([FolderInfo(ii).folder '/' FolderInfo(ii).name],20);
    yy_tmp(:,1) = ones(size(yy_tmp,1),1).*seeds(ii);
    simu_data0{ii} = yy_tmp;
end

%%

simu_data = cat(1,simu_data0{:});

% save('simu_data.mat','simu_data','-v7.3');


%% keep data only around TEB

lat_tmp = simu_data(:,11);
simu_data(lat_tmp<1,:)=[]; % removing negative latitudes

%%

lat_tmp = simu_data(:,11);
lon_tmp = simu_data(:,12);

PDG_tmp = simu_data(:,7);

to_keep = PDG_tmp~=22;
lat_tmp_leptons_only = lat_tmp(to_keep);
lon_tmp_leptons_only = lon_tmp(to_keep);

q1_lat = quantile(lat_tmp_leptons_only,0.01);
q2_lat = quantile(lat_tmp_leptons_only,0.99);
q1_lon = quantile(lon_tmp_leptons_only,0.01);
q2_lon = quantile(lon_tmp_leptons_only,0.99);

qmid_lon = (q2_lon+q1_lon)/2.0;
qmid_lat = (q2_lat+q1_lat)/2.0;

%
inside = lat_tmp>q1_lat & lat_tmp<q2_lat & lon_tmp>q1_lon & lon_tmp<q2_lon;

simu_data = simu_data(inside,:);


%%

save('simu_data_filtered.mat','simu_data','-v7.3');


%%
% exit

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
