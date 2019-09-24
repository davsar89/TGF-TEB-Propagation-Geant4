clear all
close all
clc

%%

simu_data = better_importdata('fused_15km_gauss_20deg.txt',24);

save('simu_data_closeToTEB_15km_gauss_20deg.mat','simu_data','-v7.3');

%%

function data_temp = better_importdata(full_file_path,number_of_columns)
fid = fopen(full_file_path);

text_pattern_to_scan = repmat('%f ',1,number_of_columns);

C = textscan(fid,text_pattern_to_scan,'CommentStyle','#');
data_temp = [C{:}];

fclose(fid);
end
