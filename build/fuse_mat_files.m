clear all
close all
clc

%%

FolderInfo = dir2('./output/');

nb_files=length(FolderInfo);

simu_data = cell(1,nb_files);

yy_tmp = load([FolderInfo(1).folder '/' FolderInfo(1).name]);
fields = fieldnames(yy_tmp);

%%
seeds=1:length(FolderInfo);
%%

for ii=1:length(fields)
    yy.(fields{ii}) = [];
end

%

for ii=1:length(FolderInfo)
    
    yy_tmp = load([FolderInfo(ii).folder '/' FolderInfo(ii).name]);
    yy_tmp.RANDOM_SEED = seeds(ii);
    
    for i_field = 1:length(fields)
        if (length(yy_tmp.(fields{i_field}))>1)
            yy.(fields{i_field}) = [yy.(fields{i_field}); double(yy_tmp.(fields{i_field}))];
        elseif (length(yy_tmp.(fields{i_field}))==1)
            array_tmp = ones(1,length(yy_tmp.energy)) .* double(yy_tmp.(fields{i_field}));
            yy.(fields{i_field}) = [yy.(fields{i_field}); array_tmp'];
        end
    end
end

%%
save('simu_data.mat','yy','-v7.3');





%%
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